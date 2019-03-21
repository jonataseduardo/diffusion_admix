library(data.table)
library(parallel)

## Functions 

filter_gad <-
  function(gad,
           filter_pass = TRUE, 
					 only_nonsynonymous = TRUE, 
					 remove_fixed = TRUE,
					 remove_AA = TRUE,
					 remove_mult_al = TRUE,
					 an_th = 0.7){

    setnames(gad, "GERP++_RS", "GERP_RS_score")

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

    if(filter_pass)
			gad <- gad[FILTER == "PASS"]

		if(only_nonsynonymous)
			gad <- gad[ExonicFunc.refGene == "nonsynonymous SNV"]

		#Remove alleles with ancestral information
		if(remove_AA)
			gad <- gad[AA == ""] 

    ac <- grep("AC_", pop_info, value = TRUE)
    an <- grep("AN_", pop_info, value = TRUE)
    
    #Remove fixed variants
	  if(remove_fixed){	
			gad <- 
				gad[!Reduce(`&`, Map(`==`, mget(ac), 0))
						][!Reduce(`&`, Map(`==`, mget(ac), mget(an)))]
		}

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
      grep("^AC_[a-z]{3}$", names(gad), value = TRUE)

    an <- 
      grep("^AN_[a-z]{3}$", names(gad), value = TRUE)

    af <- 
      grep("^AF_[a-z]{3}$", names(gad), value = TRUE)

    sfs[, (af) := Map(`/`, mget(ac), mget(an))] 
    return(sfs[rep(1:.N, num_sites)])
  }

project_sfs_1d <-
  function(sfs_dt, AN_new, groups = "all_groups"){
    
    if(groups == "all_groups")
      groups <- setdiff(names(sfs_dt), c("AN", "AC", "num_sites"))

    sfs_project <- 
      sfs_dt[, .(AC_new = 1:min((AN_new) -1, AC), 
                 AN_new = (AN_new)), 
              by = c(groups, c("AN", "AC", "num_sites"))
            ][, proj_i := dhyper(AC_new, AN_new, AN - AN_new, AC)
            ][, .(num_sites_new = sum(num_sites * proj_i)), 
              keyby = c(groups, c("AN_new", "AC_new"))]

    setnames(sfs_project, 
             paste0(c("AN", "AC", "num_sites"), "_new"),
             c("AN", "AC", "num_sites"))

    setcolorder(sfs_project, names(sfs_dt))

    return(sfs_project[num_sites > 0])
  }

project_onepop_sfs  <- 
  function(sfs_dt, pop, AN_new){

    old_names <- names(sfs_dt)
    setnames(sfs_dt, paste0(c("AN_", "AC_"), pop), c("AN", "AC"))

    sfs_project <- 
      project_sfs_1d(sfs_dt, AN_new)

    setnames(sfs_dt, c("AN", "AC"), paste0(c("AN_", "AC_"), pop))
    setnames(sfs_project, c("AN", "AC"), paste0(c("AN_", "AC_"), pop))

    setcolorder(sfs_project, old_names)
    return(sfs_project)
  }

## Main Code

gad_files <-
  grep("chr[0-9]{1,2}.annotated.csv", 
       list.files("../gnomad_data", full.names = TRUE), 
       value = TRUE)

# gad1 <- fread("../gnomad_data/gnomad.exomes.r2.1.sites.chr1.annotated.csv")
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

gad[,unique(REVEL_quantile)]


scols <- 
    c("Chr", "Start", "Ref","Alt", "rsid", 
			"GERP_RS_score", "GERP_RS_rank", "GERP_RS_quantile",
      "AC_afr", "AN_afr", "AC_amr", "AN_amr", "AC_nfe", "AN_nfe")

pop_list <- 
      c("afr", "amr", "nfe")

score <- "GERP_RS"

sfs <- sfs_from_ac(gad, pop_list, score)
setkeyv(sfs_small, key(sfs))
system.time(sfs_small <- project_onepop_sfs(sfs, "afr", 16256))
all.equal(sfs_small, sfs)


sfs_1d <- sfs_1d_from_ac(gad, pop_list[1], score)
system.time(sfs_small_1d <- project_sfs_1d(sfs_1d, 16256))
setkeyv(sfs_small_1d, key(sfs_1d))
all.equal(sfs_small_1d, sfs_1d)


