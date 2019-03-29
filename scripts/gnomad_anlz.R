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
system.time(make_delquantile(gad))
system.time(rescale_ac(gad))
system.time(normalize_del_score(gad))
system.time(summarise_scores_norm_quantile(gad))

scols <- 
    c("Chr", "Start", "Ref","Alt", "rsid", 
			"GERP_RS_score", "GERP_RS_rank", "GERP_RS_quantile",
      "AC_afr", "AN_afr", "AC_amr", "AN_amr", "AC_nfe", "AN_nfe")

del_info <-
  grep("_quantile$", names(gad), value = TRUE)
del_info

AF_pop1 = 'AF_afr'
AF_pop2 = 'AF_nfe'
score_quantile = 'REVEL_quantile'


gad[AF_nfe > 0.0, mean(AF_nfe)]
gad[AF_nfe > 0.0, .N]
gad[AF_afr > 0.0, .N]

combn(1:3, 2)

sort(pop_list)

pop_freq <- 
  grep("^AF_[a-z]{3}$", names(gad), value = TRUE)

pop_list <- 
  sort(gsub("(^AF_)([a-z]{3}$)", "\\2", pop_freq))

pop_comb <- combn(pop_list[c(1,2,4)], 2)
pop_comb <- combn(pop_list, 2)

scores_list <- 
  gsub("(^.+)(_quantile$)", "\\1",
    grep("_quantile", names(gad), value = TRUE))
scores_list

system.time(
ttest_data <- 
  rbindlist(
  lapply(scores_list, function(score){
    rbindlist(
    lapply(1:dim(pop_comb)[2],
           function(i){
             pop1 <- pop_comb[1,i] 
             pop2 <- pop_comb[2,i]
             ac <- paste0("AC_", c(pop1, pop2))
             an <- paste0("AN_", c(pop1, pop2))
             t_dt <- 
               apply_ttest(rm_fixed_alleles(gad, ac, an),
                           pop1, pop2, score)

             return(t_dt)
           }))
         })
  )
)


ttest_data[pop1 == 'afr' & pop2 == 'nfe' & p.value < 0.05][order(quantiles)]

scols <- c(paste0(score, '_quantile'), pop_freq)

system.time(
gad[, data.table(t(unlist(t.test(AF_afr, AF_nfe), use.names = FALSE))), by = REVEL_quantile]
)


gad[, data.table(t(unlist(wilcox.test(AF_afr, AF_nfe, conf.int = TRUE), use.names = FALSE))), by = REVEL_quantile]



sfs_afr_full <- sfs_1d_from_ac(gad, pop = 'afr', score = score)
system.time(sfs_afr_small <- project_sfs_1d(sfs_afr_full, AN_afr, folded = FALSE))
setkeyv(sfs_afr_small, key(sfs_afr_full))
all.equal(sfs_afr_small, sfs_afr_full)

ac <- paste0("AC_", pop_list)
an <- paste0("AN_", pop_list)

sfs_full <- sfs_from_ac(gad, pop_list, score)
sfs_afr_nfe <- sfs_from_ac(gad, c("afr", "nfe"), score)

sfs_afr_nfe1 <- 
  sfs_afr_nfe[score == score & 
              quantiles == "[0,0.01]"]

sfs_small_afr_nfe <- 
  project_onepop_sfs(sfs_afr_nfe[score == score & 
                     quantiles == "[0,0.01]"], 
                     "nfe", 
                      AN_afr)

sfs_afr_nfe1[, mean(AC_afr/ AN_afr * num_sites)]
sfs_small_afr_nfe[, mean(AC_afr/ AN_afr * num_sites)]
sfs_afr_nfe1[, mean(AC_nfe/ AN_nfe * num_sites)]
sfs_small_afr_nfe[, mean(AC_nfe/ AN_nfe * num_sites)]


sfs_full <- rm_fixed_alleles(sfs_full, ac, an)

sfs_small2 <- copy(sfs_small)
sfs_small2[, num_sites := round(num_sites)]
sfs_small2[num_sites > 0]

AN_afr <- sfs_full[1, AN_afr] 
sfs_small <- sfs_full[quantiles == "[0,0.01]"]
pop = "afr"
system.time(
for(pop in pop_list[2:3]){
  sfs_small <- 
    project_onepop_sfs(sfs_small, pop, AN_afr)
}
)

sfs_small
af_small <- af_from_sfs(sfs_small)
af_small1 <- af_from_sfs(sfs_small1)

sfs_small[, mean(AC_afr/ AN_afr * num_sites)] * 1000
af_small[, .(mean(AF_afr))] * 1000
af_small1[, .(mean(AF_afr))] * 1000

round(0.509098, digits = 1)
