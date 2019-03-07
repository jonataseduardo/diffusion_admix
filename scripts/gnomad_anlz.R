library(data.table)
library(parallel)

### Functions 

filter_gad <-
  function(gad){
    snv_info <-
      c("Chr","Start", "Ref","Alt","AA", "rsid")

    pop_info <- 
      grep("^A[NCF]_[a-z]{3}$", names(gad), value = TRUE)

    del_info <-
      grep("(_score$)|([_-]raw$)|(_RS$)", names(gad), value = TRUE)

    gene_info <- 
      c("Gene.refGene")

    cols_info  <- c(snv_info, gene_info, del_info, pop_info)

    gad <- 
      gad[(ExonicFunc.refGene == "nonsynonymous SNV") & 
                     (FILTER == "PASS"), 
                     cols_info, with = FALSE]
    gad[, (del_info) := lapply(.SD, as.numeric), .SDcols = del_info]
    return(gad)
  }



make_delrank <- 
  function(gad){
    del_info <-
      grep("(_score$)|([_-]raw$)|(_RS$)", names(gad), value = TRUE)

    del_rank <- 
      gsub("^(.+)(_score$)|([_-]raw$)", "\\1_rank", del_info)

    drank <- 
      function(column)
        frank(column, na.last = FALSE)
  
    gad[, (del_rank) := lapply(.SD, drank), .SDcols = del_info]
    gad[, (del_rank) := 1 - .SD / .N , .SDcols = del_rank]

    return(gad)
  }


x <- gad[idx >= 2520 & idx <= 2525]

make_delrank(x)
x


##### Main Code

gad_files <-
  grep("chr[0-9]{1,2}.annotated.csv", 
       list.files("../gnomad_data", full.names = TRUE), 
       value = TRUE)

#gad <- filter_gad(fread("../gnomad_data/gnomad.exomes.r2.1.sites.chr1.annotated.csv")) 
system.time(
gad <- 
  rbindlist(mclapply(gad_files, 
                     function(x){filter_gad(fread(x))},
                     mc.cores = 7)
            )
)



setkey(gad, Chr, Start, Ref)
gad

