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


summarise_scores_norm_quantile <- 
  function(gad, wd = 0.01)
  {

    norm_scores <- 
      grep("_norm", names(gad), value = TRUE) 

    score_names <- gsub("(^.+)(_norm$)", "\\1", norm_scores)

    nq_sumry <- 
      gad[, lapply(.SD, quantile, seq(0, 1, wd), na.rm = TRUE), 
          .SDcols = norm_scores]

    setnames(nq_sumry, norm_scores, score_names)

    nq_sumry[, quant := seq(0, 1, wd)]
    long_nq <- 
      melt(nq_sumry, id.vars = c("quant"), 
           variable.name = 'annotation', value.name = 'norm_score')

    qual <- 
      long_nq[ quant > 0.75, .(quality = round(100 * sum(quant - norm_score))), 
              by = annotation]

    long_nq <- 
      long_nq[qual, on = .(annotation)]

    setkey(long_nq, quality, annotation)
    long_nq[, quality_rank := .GRP, by = .(quality, annotation)] 
    long_nq[, .GRP, by = .(quality, annotation) ]

    ggplot(long_nq) + 
      geom_line(aes(x = norm_score, y = quant, color = quality_rank))
    ggsave('normalized_scores_quantililes.png')
  }



scols <- 
    c("Chr", "Start", "Ref","Alt", "rsid", 
			"GERP_RS_score", "GERP_RS_rank", "GERP_RS_quantile",
      "AC_afr", "AN_afr", "AC_amr", "AN_amr", "AC_nfe", "AN_nfe")

pop_list <- c("afr", "amr", "nfe")

score <- "Eigen"


sfs_afr_full <- sfs_1d_from_ac(gad, pop = 'afr', score = score)
system.time(sfs_afr_small <- project_sfs_1d(sfs_afr_full, AN_afr, folded = FALSE))
setkeyv(sfs_afr_small, key(sfs_afr_full))
all.equal(sfs_afr_small, sfs_afr_full)




ac <- paste0("AC_", pop_list)
an <- paste0("AN_", pop_list)



sfs_full <- sfs_from_ac(gad, pop_list, score)

sfs_afr_nfe <- sfs_from_ac(gad, c("afr", "nfe"), score)



sfs_afr_nfe1 <- sfs_afr_nfe[score == score & quantiles == "[0,0.01]"]
sfs_small_afr_nfe <- 
  project_onepop_sfs(sfs_afr_nfe[score == score & quantiles == "[0,0.01]"] , "nfe", AN_afr)

sfs_afr_nfe1[, mean(AC_afr/ AN_afr * num_sites)]
sfs_small_afr_nfe[, mean(AC_afr/ AN_afr * num_sites)]
sfs_afr_nfe1[, mean(AC_nfe/ AN_nfe * num_sites)]
sfs_small_afr_nfe[, mean(AC_nfe/ AN_nfe * num_sites)]
sfs_small_afr_nfe


sfs_full <- rm_fixed_alleles(sfs_full, ac, an)
sfs_full

sfs_small2 <- copy(sfs_small)
sfs_small2[, num_sites := round(num_sites)]
sfs_small2[num_sites > 0]
sfs_full

AN_afr <- sfs_full[1, AN_afr] 
sfs_small <- sfs_full[quantiles == "[0,0.01]"]
pop = "afr"
system.time(
for(pop in pop_list[2:3]){
  sfs_small <- project_onepop_sfs(sfs_small, pop, AN_afr)
}
)

sfs_small
af_small <- af_from_sfs(sfs_small)
af_small1 <- af_from_sfs(sfs_small1)

sfs_small[, mean(AC_afr/ AN_afr * num_sites)] * 1000
af_small[, .(mean(AF_afr))] * 1000
af_small1[, .(mean(AF_afr))] * 1000

round(0.509098, digits = 1)

sfs_small[, .(max(AC_afr), max(AC_amr), max(AC_nfe))]
sfs_small1[, .(max(AC_afr), max(AC_amr), max(AC_nfe))]
sfs_small2[, .(max(AC_afr), max(AC_amr), max(AC_nfe))]

sfs_small[,sum(num_sites), AC_afr]
sfs_small[, .N]
sfs_small[num_sites == 1, .N]

an


