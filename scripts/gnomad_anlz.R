library(data.table)
library(parallel)
library(ggplot2)
source("gnomad_ultis.R")

## Main Code
gad_files <-
  grep("chr[0-9]{1,2}.annotated.csv", 
       list.files("../gnomad_data", full.names = TRUE), 
       value = TRUE)

gad1 <- fread("../gnomad_data/gnomad.exomes.r2.1.sites.chr1.annotated.csv")

gad11 <- filter_gad(gad1, 
                    filter_pass = TRUE, 
                    only_nonsynonymous = FALSE, 
                    remove_fixed = TRUE,
                    remove_AA = FALSE,
                    remove_mult_al = TRUE,
                    score_to_numeric = TRUE)

system.time(
gad <- 
  rbindlist(mclapply(gad_files, 
                     function(x){filter_gad(fread(x))},
                     mc.cores = 11)
           )
)

setkey(gad, Chr, Start, Ref)

del_info <-
  grep("_quantile$", names(gad), value = TRUE)

pop_freq <- 
  grep("^AF_[a-z]{3}$", names(gad), value = TRUE)

pop_list <- 
  sort(gsub("(^AF_)([a-z]{3}$)", "\\2", pop_freq))

AN_afr <- 
  gad[, max(AN_afr)]

gad_small <- 
  mc_downsampling(gad, pop_list[2:5], AN_afr)

system.time(rescale_ac(gad))
system.time(make_delrank(gad))
system.time(normalize_del_score(gad))
system.time(make_delquantile(gad, breaks = c(seq(0, 0.5, 0.1), 1)))
system.time(summarise_scores_norm_quantile(gad))

system.time(rescale_ac(gad_small))
system.time(make_delrank(gad_small))
system.time(normalize_del_score(gad_small))
system.time(make_delquantile(gad_small, breaks = c(seq(0, 0.25, 0.05), 1)))
system.time(summarise_scores_norm_quantile(gad_small, fname = 'teste_small.png'))

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
               apply_ttest(rm_fixed_alleles(gad_small, ac, an),
                           pop1, pop2, score, inplace = TRUE)
             return(t_dt)
           }, mc.cores = length(scores_list)))
         })
  )
)

gad1_small <- 
  mc_downsampling(gad11, 
                  pop_list[2:5], 
                  AN_afr)

gad1_small[AA != "", .(mean(AF_afr), mean(AF_nfe), .N), 
           by = ExonicFunc.refGene]

gad1_small[AA  == "" & AF_afr < 0.05 & AF_nfe < 0.05  , 
           .(mean(AF_afr), mean(AF_nfe)),
           by = ExonicFunc.refGene]

gad1_small[, .(mean(AF_afr), mean(AF_nfe)), 
           by = ExonicFunc.refGene]

gad1_small[(Ref == AA | Alt == AA) & 
           ExonicFunc.refGene == 'nonsynonymous SNV', 
           t.test(AF_afr, AF_nfe)]

gad1_small[AA == "" &  AF_afr  < 0.05 & AF_nfe  < 0.05
           ExonicFunc.refGene == 'nonsynonymous SNV', 
           t.test(AF_afr, AF_nfe)]

ttest_afr_nfe <- 
  ttest_data[pop1 == 'afr' & pop2 == 'nfe']

ttest_afr_nfe[mean_AF_pop1 < mean_AF_pop2]

long_tt <- 
  melt(ttest_afr_nfe, 
       measure = patterns(c("var_", "mean_")),
       variable.name = 'pop',
       value.name = c("var_AF", "mean_AF")
       )

slist <- scores_list[c(2,22,15,24)]
slist <- c('VEST3', 'M-CAP', 'MutPred', 'DANN')

long_tt[pop == 1, pop := pop1]
long_tt[pop == 2, pop := pop2]

{
ggplot(long_tt[score %in% slist]) + 
  geom_errorbar(aes(x = quantiles, 
                    ymax = mean_AF + 1.96 * sqrt(var_AF / num_sites),
                    ymin = mean_AF - 1.96 * sqrt(var_AF / num_sites),
                    color = pop), 
                position = position_dodge()) + 
  geom_text(data = long_tt[score %in% slist][pop == 'afr'],
            aes(x = quantiles, 
                 y = mean_AF + 10.0 * sqrt(var_AF / num_sites), 
                 label = plabel(p.value),
                 vjust = 'top'),
                size = 2.5 
             ) + 
  geom_point(aes(x = quantiles, y = mean_AF, color = pop), 
             position = position_dodge(0.9)) + 
  labs(y = 'mean(AF)') + 
  facet_wrap( ~score) + 
  theme_classic()  
  ggsave('mean_AF_afr_nfe_small.png')
}

ttest_afr_nfe[, .N, by = score]
sfs_afr_full <- 
  sfs_1d_from_ac(gad, pop = 'afr', score = score)

system.time(
  sfs_afr_small <- 
    project_sfs_1d(sfs_afr_full, AN_afr, folded = FALSE)
)

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

conditions <- 
  c('all', 'AC_afr > 0', 'AC_nfe > 0', '(AC_afr > 0) & (AC_nfe > 0)')

sfs_small <- 
  rbindlist(
  list(sfs_afr_nfe[quantiles == "[0,0.1]"], 
       sfs_small_afr_nfe))

system.time(
mdt <- 
  rbindlist(
  lapply(
    -1L:4L, 
    function(digits){
    rbindlist(
    lapply(
     conditions, 
     function(cond){
       emv_sfs_afr_nfe(sfs_small, 
                       conditional = cond, 
                       digits = digits)
     }
    ))
    }
  ))
)

setkey(mdt, score, quantiles, AN_nfe, digits, conditional)

mdt_long <- 
  melt(mdt, measure = patterns(c("var_", "mean_")),
       variable.name = 'pop',
       value.name = c("var_AF", "mean_AF")
       )

mdt_long[pop == 1, pop := 'afr']
mdt_long[pop == 2, pop := 'nfe']
mdt_long[digits == -1, digits := 5] 
mdt_long[, digits := factor(digits)]
mdt_long[, conditional := factor(conditional, levels = conditions)]

{
ggplot(mdt_long) + 
  geom_line(aes(x = AN_nfe, y = total_num_sites, color = digits, linetype = pop)) + 
  geom_point(aes(x = AN_nfe, y = total_num_sites, color = digits, shape = pop)) + 
  theme_classic() + 
  scale_color_brewer(palette = 'Spectral') + 
  scale_x_continuous(breaks = seq(AN_afr, AN_nfe, delta_AN)) + 
  facet_wrap( ~ conditional, scales = 'free_y')
  ggsave('numsites_project_nfe.png')
}

