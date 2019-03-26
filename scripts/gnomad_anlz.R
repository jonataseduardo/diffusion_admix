library(data.table)
library(parallel)
install.packages("normalize")

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

norm_scores <- 
  grep("_norm", names(gad), value = TRUE) 

score_names <- gsub("(^.+)(_norm$)", "\\1", norm_scores)

gad[, quantile(GERP_RS_norm, seq(0,1, 0.2), na.rm = TRUE)] 



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


