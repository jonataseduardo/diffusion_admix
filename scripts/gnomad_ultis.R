library(data.table)
library(parallel)
library(ggplot2)

## Functions 

rm_fixed_alleles <-
  function(ac_an_dt, ac, an){
   return(ac_an_dt[!Reduce(`&`, Map(`==`, mget(ac), 0))
                  ][!Reduce(`&`, Map(`==`, mget(ac), mget(an)))])
  }

filter_gad <-
  function(gad,
           filter_pass = TRUE, 
           only_nonsynonymous = TRUE, 
           remove_fixed = TRUE,
           remove_AA = TRUE,
           remove_mult_al = TRUE,
           score_to_numeric = TRUE, 
           an_th = 0.7){

    if(filter_pass)
      gad <- gad[FILTER == "PASS"]

    if(only_nonsynonymous)
      gad <- gad[ExonicFunc.refGene == "nonsynonymous SNV"]

    #Remove alleles with ancestral information
    if(remove_AA)
      gad <- gad[AA == ""] 


    #Changing some names in the annotated tables
    phy <- grep("(vertebrate$)|(mammalian$)|(logOdds$)", 
                names(gad), value = TRUE)
    phy_scores <- paste0(phy, "_score")
    setnames(gad, c("GERP++_RS", phy), c("GERP_RS_score", phy_scores))

    snv_info <-
      c("Chr", "Start", "Ref","Alt","AA", "rsid")

    pop_info <- 
      grep("^A[NCF]_[a-z]{3}$", names(gad), value = TRUE)

    del_info <-
      grep("(_score$)|([_-]raw$)", names(gad), value = TRUE)

    gene_info <- 
      c("Gene.refGene")

    cols_info  <- c(snv_info, gene_info, del_info, pop_info)

    gad <- gad[, ..cols_info]

    ac <- grep("AC_", pop_info, value = TRUE)
    an <- grep("AN_", pop_info, value = TRUE)

    #Remove fixed variants
    if(remove_fixed) 
      gad <- rm_fixed_alleles(gad, ac, an)

    #remove multiallelic snps
    if(remove_mult_al){
      gad[, N := .N, by = .(Chr, Start, Ref)]
      gad <- gad[N == 1]
      gad[, N := NULL]
    }

    #remove snps secreend in less an_th proportion of samples
    if(an_th < 0.99){
      max_an <- gad[, lapply(.SD, max), .SDcols = an]
      max_an <- an_th * max_an
      gad <- gad[Reduce(`&`, Map(`>=`, mget(an), max_an[, mget(an)]))]
    }

    if(score_to_numeric)
      gad[, (del_info) := lapply(.SD, as.numeric), .SDcols = del_info]

    return(gad)
  }

make_delrank <- 
  function(gad,
           inplace = TRUE){

    if(!inplace)
      gad <- copy(gad)

    del_info <-
      grep("(_score$)|([_-]raw$)", names(gad), value = TRUE)

    del_rank <- 
      gsub("^(.+)(_score$)|([_-]raw$)", "\\1_rank", del_info)

    gad[, (del_rank) := lapply(.SD, frank, na.last = FALSE), .SDcols = del_info]
    gad[, (del_rank) := 1.0 - .SD / .N , .SDcols = del_rank]

    return(gad)
  }

normalize <-
  function(x){
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    if((x_min == 0) & (x_max == 1)){
      return(x)
    }else{
      x_wd <- x_max - x_min
      return((x - x_min)/x_wd)
    }
  }

normalize_del_score <- 
  function(gad,
           inplace = TRUE){

    if(!inplace)
      gad <- copy(gad)

    del_info <-
      grep("(_score$)|([_-]raw$)", names(gad), value = TRUE)

    del_norm <- 
      gsub("^(.+)(_score$)|([_-]raw$)", "\\1_norm", del_info)

    gad[, (del_norm) := lapply(.SD, normalize), .SDcols = del_info]

    return(gad)
  }

make_delquantile <- 
  function(gad, 
					 breaks = c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 1.0),
           inplace = TRUE){
    if(!inplace)
      gad <- copy(gad)

    del_rank <-
      grep("_rank$", names(gad), value = TRUE)

    del_quant <- 
      gsub("^(.+)(_rank$)", "\\1_quantile", del_rank)

    gad[, (del_quant) := lapply(.SD, cut, breaks = breaks, include.lowest = TRUE), 
				.SDcols = del_rank]
    return(gad)
  }

rescale_ac <- 
  function(gad, 
           eval_maf = TRUE,
           inplace = TRUE){
    if(!inplace)
      gad <- copy(gad)

    ac <- 
      grep("^AC_[a-z]{3}$", names(gad), value = TRUE)

    an <- 
      grep("^AN_[a-z]{3}$", names(gad), value = TRUE)

    af <- 
      grep("^AF_[a-z]{3}$", names(gad), value = TRUE)

    max_an <- gad[, lapply(.SD, max), .SDcols = an]

    gad[, (af) := Map(`/`, mget(ac), mget(an))] #
    gad[, (an) := max_an] 

    if(eval_maf){
      for(i in af){
          gad[get(i) > 0.5, (i) := 1.0 - get(i)]
      }
    }
    gad[, (ac) := Map(`*`, mget(af), mget(an))] 
    gad[, (ac) := lapply(.SD, round), .SDcols = ac] 

    return(gad)
  }

sfs_1d_from_ac <-
  function(gad, pop, score){
    score_quantile = paste0(score, "_quantile")
    ac_pop = paste0("AC_", pop)
    an_pop = paste0("AN_", pop)

    sfs <- gad[, .(score = (score), pop = (pop), num_sites = .N), 
               keyby = c(score_quantile, an_pop, ac_pop)]

    names(sfs) <- c("quantiles", "AN", "AC", "score", "pop", "num_sites")
    setcolorder(sfs, c("score", "quantiles", "pop", "AN", "AC", "num_sites"))
    return(sfs)
  }

sfs_from_ac <-
  function(gad, pop_list, score){
    score_quantile = paste0(score, "_quantile")
    ac = paste0("AC_", pop_list)
    an = paste0("AN_", pop_list)

    sfs <- gad[, .(score = (score), num_sites = .N), 
               keyby = c(score_quantile, an, ac)]

    names(sfs) <- c("quantiles", an, ac, "score", "num_sites")
    setcolorder(sfs, c("score", "quantiles", an, ac, "num_sites"))
    return(sfs)
  }

af_from_sfs_1d <- 
  function(sfs){
    return(sfs[, .(AF = rep(AC/AN, num_sites)), 
               by = .(score, quantiles, pop)])
  }

af_from_sfs <- 
  function(sfs){
    ac <- 
      grep("^AC_[a-z]{3}$", names(sfs), value = TRUE)

    an <- 
      grep("^AN_[a-z]{3}$", names(sfs), value = TRUE)

    af <- 
      gsub("(^AN)(_[a-z]{3}$)", "AF\\2", an)

    sfs[, (af) := Map(`/`, mget(ac), mget(an))] 
    return(sfs[rep(1:.N, num_sites)])
  }

project_sfs_1d <-
  function(sfs_dt, AN_new, groups = "all_groups", folded = TRUE){
    
    if(groups == "all_groups")
      groups <- setdiff(names(sfs_dt), c("AN", "AC", "num_sites"))

    if(folded){
      AN_new_fold <- AN_new / 2
    }else{
      AN_new_fold <- AN_new
    }

    sfs_project <- 
      sfs_dt[, .(AC_new = 0:min((AN_new), AC), 
                 AN_new = (AN_new),
                 num_sites = num_sites), 
              by = c(groups, c("AN", "AC"))
             ][, proj_i := dhyper(AC_new, AN_new, AN - AN_new, AC)
             ][, .(num_sites_new = round(sum(num_sites * proj_i), digits = 0)), 
               keyby = c(groups, c("AN_new", "AC_new"))]

    setnames(sfs_project, 
             paste0(c("AN", "AC", "num_sites"), "_new"),
             c("AN", "AC", "num_sites"))

    setcolorder(sfs_project, names(sfs_dt))

    return(sfs_project[num_sites > 0])
  }

project_onepop_sfs  <- 
  function(sfs_dt, pop, AN_new, folded = TRUE){

    old_names <- names(sfs_dt)
    setnames(sfs_dt, paste0(c("AN_", "AC_"), pop), c("AN", "AC"))

    sfs_project <- 
      project_sfs_1d(sfs_dt, AN_new, folded = folded)

    setnames(sfs_dt, c("AN", "AC"), paste0(c("AN_", "AC_"), pop))
    setnames(sfs_project, c("AN", "AC"), paste0(c("AN_", "AC_"), pop))

    setcolorder(sfs_project, old_names)

    ac <- grep("AC_", names(sfs_project), value = TRUE)
    an <- grep("AN_", names(sfs_project), value = TRUE)
    sfs_project <-
      rm_fixed_alleles(sfs_project, ac, an)  

    return(sfs_project)
  }

summarise_scores_norm_quantile <- 
  function(gad, 
           wd = 0.01, 
           quant_th = 0.1){

    norm_scores <- 
      grep("_norm", names(gad), value = TRUE) 

    score_names <- 
      gsub("(^.+)(_norm$)", "\\1", norm_scores)

    nq_sumry <- 
      gad[, lapply(.SD, quantile, seq(0, 1, wd), na.rm = TRUE), 
          .SDcols = norm_scores]

    setnames(nq_sumry, norm_scores, score_names)

    nq_sumry[, quant := seq(0, 1, wd)]

    long_nq <- 
      melt(nq_sumry, id.vars = c("quant"), 
           variable.name = 'annotation', value.name = 'norm_score')

    qual <- 
      long_nq[quant > 1 - quant_th, 
              .(quality = round(100 * sum(quant - norm_score))), 
              by = annotation]

    long_nq <- 
      long_nq[qual, on = .(annotation)]

    setkey(long_nq, quality, annotation)

    annot_level <- 
      long_nq[, .GRP, keyby = .(quality, annotation)][, annotation]

    long_nq[, annotation := factor(annotation, levels = annot_level)]
    long_nq[, quality_rank := .GRP, keyby = .(quality, annotation)] 

    make_plot <- 
      function(){
        ggplot(long_nq, aes(y = norm_score, x = quant, colour = annotation)) + 
          theme_classic() +
          scale_colour_viridis_d() + 
          guides(colour = guide_legend(ncol=4)) + 
          labs(x = 'score quantile', y = 'normalized score') +
          theme(legend.position = 'bottom', 
                legend.title = element_blank()) + 
          geom_line()
        ggsave('normalized_scores_quantililes.png')
      }

    make_plot()
  }

apply_ttest <-
  function(gad, pop1, pop2, score, inplace = FALSE){

    if(inplace)
      gad <- copy(gad)

    AF_pop1 = paste0('AF_', pop1)
    AF_pop2 = paste0('AF_', pop2)
    score_quantile = paste0(score, '_quantile')

    setnames(gad, c(AF_pop1, AF_pop2), c("AF_x", "AF_y"))
    ttest_dt <- 
      gad[, data.table(t(unlist(t.test(AF_x, AF_y)))), 
            keyby = c(score_quantile)][, c(1:2, 4:8)]

    setnames(ttest_dt, 
             c(score_quantile, "estimate.mean of x", "estimate.mean of y"),
             c('quantiles', 'mean_AF_pop1', 'mean_AF_pop2'))


    old_t_names <- names(ttest_dt)
    new_t_names <- c('score', 'pop1', 'pop2', old_t_names) 

    ttest_dt[, `:=`(score = score, pop1 = pop1, pop2 = pop2)]

    setcolorder(ttest_dt, new_t_names)

    setnames(gad, c("AF_x", "AF_y"), c(AF_pop1, AF_pop2))

    return(ttest_dt[])
  }
