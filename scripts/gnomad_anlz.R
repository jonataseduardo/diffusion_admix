library(data.table)


gad_files <-
  grep("head.annotated.csv", 
       list.files("../gnomad_data", full.names = TRUE), 
       value = TRUE)
  
gad <- 
  rbindlist(lapply(gad_files, fread))

gad[(ExonicFunc.refGene == "nonsynonymous SNV"), .N, .(FILTER)]
gad[(ExonicFunc.refGene == "nonsynonymous SNV") & (FILTER == "PASS")]
