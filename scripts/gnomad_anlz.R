library(data.table)
library(parallel)
library(ggplot2)

source("gnomad_ultis.R")

## Main Code

gad_files <-
  grep("chr[0-9]{1,2}.annotated.csv", 
       list.files("../gnomad_data", full.names = TRUE), 
       value = TRUE)

#gad1 <- fread("../gnomad_data/gnomad.exomes.r2.1.sites.chr1.annotated.csv")
system.time(
gad <- 
  rbindlist(mclapply(gad_files, 
                     function(x){filter_gad(fread(x))},
                     mc.cores = 11)
           )
)

setkey(gad, Chr, Start, Ref)
system.time(make_delrank(gad))
system.time(rescale_ac(gad))
system.time(normalize_del_score(gad))
system.time(make_delquantile(gad, breaks = c(seq(0, 0.5, 0.1), 1)))

system.time(summarise_scores_norm_quantile(gad))

del_info <-
  grep("_quantile$", names(gad), value = TRUE)

pop_freq <- 
  grep("^AF_[a-z]{3}$", names(gad), value = TRUE)

pop_list <- 
  sort(gsub("(^AF_)([a-z]{3}$)", "\\2", pop_freq))

pop_comb <- combn(pop_list[c(1,2,4)], 2)
pop_comb <- combn(pop_list, 2)

scores_list <- 
  gsub("(^.+)(_quantile$)", "\\1",
    grep("_quantile", names(gad), value = TRUE))

system.time(
ttest_data <- 
  rbindlist(
  lapply(scores_list, function(score){
    rbindlist(
    mclapply(1:dim(pop_comb)[2],
           function(i){
             pop1 <- pop_comb[1,i] 
             pop2 <- pop_comb[2,i]
             ac <- paste0("AC_", c(pop1, pop2))
             an <- paste0("AN_", c(pop1, pop2))
             t_dt <- 
               apply_ttest(rm_fixed_alleles(gad, ac, an),
                           pop1, pop2, score, inplace = TRUE)

             return(t_dt)
           }, mc.cores = length(scores_list)))
         })
  )
)

setkey(ttest_data1, score, pop1, pop2, quantiles)

ttest_data1 <- 
  copy(ttest_data)

ttest_data[pop1 == 'afr' & pop2 == 'nfe' & p.value < 0.05][order(quantiles)]
ttest_data[pop1 == 'afr' & pop2 == 'nfe' & p.value < 0.05][, .N, by = .(score)]

ttest_data[p.value < 0.05][pop1 == 'afr']

ttest_data[pop1 == 'afr' & mean_AF_pop1 < mean_AF_pop2 & p.value < 0.5]

ttest_data[score == 'Polyphen2_HDIV' & p.value < 0.05]

scols <- c(paste0(score, '_quantile'), pop_freq)

scols <- 
    c("Chr", "Start", "Ref","Alt", "rsid", 
			"GERP_RS_score", "GERP_RS_rank", "GERP_RS_quantile",
      "AC_afr", "AN_afr", "AC_amr", "AN_amr", "AC_nfe", "AN_nfe")


score = scores_list[21]

sfs_afr_full <- sfs_1d_from_ac(gad, pop = 'afr', score = score)

system.time(sfs_afr_small <- project_sfs_1d(sfs_afr_full, AN_afr, folded = FALSE))
setkeyv(sfs_afr_small, key(sfs_afr_full))
all.equal(sfs_afr_small, sfs_afr_full)

ac <- paste0("AC_", pop_list)
an <- paste0("AN_", pop_list)

sfs_full <- sfs_from_ac(gad, pop_list, score)

sfs_afr_nfe <- 
  rm_fixed_alleles(
    sfs_from_ac(gad, c("afr", "nfe"), score),
    paste0('AC_', c("afr", "nfe")),
    paste0('AN_', c("afr", "nfe")))

AN_afr <- sfs_afr_nfe[1, AN_afr]
AN_nfe <- sfs_afr_nfe[1, AN_nfe]
delta_AN <- round((AN_nfe - AN_afr) / 4)

system.time(
sfs_small_afr_nfe <- 
  rbindlist(
  mclapply(
    seq(AN_afr, AN_nfe - delta_AN, delta_AN), 
    function(AN_new){
      project_onepop_sfs(sfs_afr_nfe[quantiles == "[0,0.1]"], 
                         "nfe", 
                          AN_new,
                          inplace = TRUE)
    }, 
    mc.cores = 4
  ))
)

sfs_afr_nfe[AC_afr > 0, .(sum((AC_afr / AN_afr) * num_sites) / sum(num_sites), 
                sum((AC_nfe / AN_nfe) * num_sites) / sum(num_sites)),
                by = .(quantiles, AN_nfe)]

sfs_afr_nfe[, 
            .(sum((AC_afr / AN_afr) * num_sites) / sum(num_sites), 
              sum((AC_nfe / AN_nfe) * num_sites) / sum(num_sites)),
                by = quantiles]


sfs_small_afr_nfe[, 
                  .(sum((AC_afr / AN_afr) * num_sites) / sum(num_sites), 
                    sum((AC_nfe / AN_nfe) * num_sites) / sum(num_sites)),
                   by = .(quantiles, AN_nfe)]

sfs_small_afr_nfe[AC_afr > 0, 
                  .(sum((AC_afr / AN_afr) * num_sites) / sum(num_sites), 
                    sum((AC_nfe / AN_nfe) * num_sites) / sum(num_sites)),
                   by = .(quantiles, AN_nfe)]



conditional = 'all'
make_round = F
digits = 0

emv_sfs_afr_nfe  <-
  function(sfs_dt, conditional = 'all', make_round = FALSE, digits = 0L){

    if(conditional == 'all'){
      sfs_dt <- copy(sfs_dt)
    }else{
      sfs_dt <- sfs_dt[eval(parse(text=conditional))]
    }

    if(make_round){
      sfs_dt <- 
        sfs_dt[, num_sites := round(num_sites, digits = digits)
               ][num_sites > 0]
    }else{
      digits = -1L
    }

    sfs_dt[, `:=`(AF_afr = AC_afr / AN_afr, AF_nfe = AC_nfe / AN_nfe)]
    mean_dt <- 
      sfs_dt[,.(mean_AF_afr = sum(AF_afr * num_sites) / sum(num_sites), 
                mean_AF_nfe = sum(AF_nfe * num_sites) / sum(num_sites),
                total_num_sites = sum(num_sites),
                conditional = conditional,
                digits = digits),
            by = .(score, quantiles, AN_nfe)]

    var_dt <- 
      sfs_dt[mean_dt, on = .(score, quantiles, AN_nfe)
            ][, .(var_AF_afr = sum((AF_afr - mean_AF_afr) ^ 2) / sum(num_sites) ^ 2,
                  var_AF_nfe = sum((AF_nfe - mean_AF_nfe) ^ 2) / sum(num_sites) ^ 2),
                by = .(score, quantiles, AN_nfe)]

    out_dt <- var_dt[mean_dt, on = .(score, quantiles, AN_nfe)]

    return(out_dt)
  }

mdt <- 
  emv_sfs_afr_nfe(sfs_small_afr_nfe, 
                  conditional = 'AC_afr > 0', 
                  make_round = TRUE, digits = 4)

mdt
