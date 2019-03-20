library(data.table)


theta <- 1000 

n <- 10
sfs_old <- 
  theta / 1:n

m <- 5 
sfs_new <- 
  sapply(1:(m -1), 
  function(i, sfs_old){
    sum(dhyper(i, m, n - m, i:(n - m + i)) * sfs_old[i:(n - m + i)])
  }, sfs_old)

print(sfs_old)
print(sfs_new)

sfs_dt <- 
  data.table(AN = (AN_old), AC = 1:AN_old, num_sites = (theta)/1:(AN_old))

sfs_dt1 <- 
  data.table(AN = (AN_old), 
             AC = sort(sample(1:AN_old, 6)))[, num_sites := (theta)/AC]


AN_new <- 5 
AN_old <- sfs_dt[1, AN]

sfs_dt2 <- sfs_dt[c(1:.N, 1:.N)][, sel := rep(c("a", "b"), each = 20)]
sfs_dt2


sfs_dt2[,.(AN = (AN_new),
          AC = 1:(AN_new - 1),
          num_sites = {
            sfs_old <- data.table(AC = AC, num_sites = num_sites);
            n <- (AN_old);
            m <- (AN_new);
            sfs_new <- 
              sapply(1:(m - 1), 
                     function(i){
                       j_idx <- sfs_old[, .I[which(AC >= i & AC <= (n - m + i))]]
                       j <- sfs_old[j_idx, AC]
                       sfs_old_in  <- sfs_old[j_idx, num_sites]
                       sum(dhyper(i, m, n - m, j) * sfs_old_in)
                     });
            sfs_new
          }), by = sel]


project_sfs_1d <-
  function(sfs_dt, AN_new, groups = "all_groups"){
    
    if(groups == "all_groups")
      groups <- setdiff(names(sfs_dt), c("AN", "AC", "num_sites"))

    sfs_project <- 
      sfs_dt[,{
              sfs_old <- data.table(AC = AC, num_sites = num_sites);
              n <- max(AN);
              m <- (AN_new);
              max_ac <- min((AN_new) - 1, max(AC));
              sfs_new <- 
                unlist(lapply(1:max_ac, 
                       function(i){
                         j_idx <- sfs_old[, .I[which(AC >= i & AC <= (n - m + i))]]
                         j <- sfs_old[j_idx, AC]
                         sfs_old_in  <- sfs_old[j_idx, num_sites]
                         sum(dhyper(i, m, n - m, j) * sfs_old_in)
                       }));
              .(AN = rep((AN_new), max_ac),
                AC = 1:max_ac,
                num_sites = sfs_new)
              }
              , by = groups]

    return(sfs_project[num_sites > 0])
  }
