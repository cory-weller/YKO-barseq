library(data.table)
allchr <- readRDS("/Users/giovanettism/Documents/R directory/1011 genomes/allchr.RDS")
candidates <- fread("/Users/giovanettism/Documents/R directory/1011 genomes/candidates.txt")
 
#objective: count differences between all strains (locations where A="0/0", B="1/1" or A="1/1", B-"0/0")
setkey(allchr, CHROM, POS)
# extract strain names as every column from 3 and on
strain_names <- colnames(allchr)[7:length(colnames(allchr))]

window_step <- 5e3    # 500 e3 is 500 kb
window_size <- 100e3      # 1e6 is 1 mb
last_window <- 11.5e6 + 1 

strain_pairs <- CJ("A" = strain_names, "B" = strain_names)[A != B]

# initialize trigger
markIdentical <- FALSE

# outer loop: for each chromosome
for(chromosome in unique(allchr$CHROM)) {
  # middle loop: for every window
  for(window_start in seq(1, last_window, window_step)) {
    # calculate end of window based on defined paramaters
    window_stop <- window_start + window_size
    
    # create a temporary subset object based on these positions
    all.sub <- allchr[CHROM == chromosome & POS >= window_start & POS < window_stop]
    # innermost loop: for every pair of strains
    # iterate over every row number of strain_pairs
    for (i in 1:nrow(strain_pairs)) {
      # extract strain names from the strain_pairs table
      strainA = strain_pairs[i, A]
      strainB = strain_pairs[i, B]
      
      all.sub2 <- all.sub[get(strainA) %in% c("0/0", "1/1")  & get(strainB) %in% c("0/0", "1/1")]
    # test whether the two columns are identical
      strains_identical <- identical(all.sub2[[strainA]], all.sub2[[strainB]])
      
      # if the two columns ARE identical
      if(strains_identical == TRUE) {
        # set trigger 
        markIdentical <- TRUE
        # break out of current `for` loop
        # this means stop iterating over strain pairs
        break
      }
    }
    
    if(markIdentical == TRUE) {
      # reset trigger
      markIdentical <- FALSE
      # break out of current `for` loop
      # e.g. stop iterating over this window and go to the next
      next
    }
    # You only get here if you go through all pairs and none are identical
    # because markIdentical has stayed FALSE the entire time
    # e.g. this window uniquely identifies every strain
    print(paste(chromosome, window_start, window_stop))
  }
}


library(foreach)
#generate strain pairs without duplicates and set column names
strain_pairs <- as.data.table(t(combn(strain_names, 2)))
setnames(strain_pairs, c("A", "B"))
#sim modifies to count 
divergence <- foreach(i =1:nrow(candidates), .combine="rbind") %do% {
    chromosome <- candidates[[i, 1]]
    start <- candidates[[i, 2]]
    stop <- candidates[[i, 3]]
  # create a temporary subset object based on these positions
  all.sub <- allchr[CHROM == chromosome & POS >= start & POS < stop]
  # innermost loop: for every pair of strains
  # iterate over every row number of strain_pairs
  foreach (j = 1:nrow(strain_pairs), .combine="rbind") %do% {
    # extract strain names from the strain_pairs table
    strainA = strain_pairs[j, A]
    strainB = strain_pairs[j, B]
    
    ndifferences <- nrow(all.sub[(get(strainA) == "0/0" & get(strainB) == "1/1") | (get(strainA) == "1/1" & get(strainB) == "0/0")])
    data.table("CHROM" = chromosome,
               "start" = start,
               "stop" = stop, 
               "strainA" = strainA,
               "strainB" = strainB,
               "N" = ndifferences)
      }
}
divergence.ag <- divergence[, list("meanN" = mean(N)), by=list(start, CHROM)]
ggplot(data = divergence.ag, aes(x=start, y=meanN)) + geom_point(shape=21, alpha=0.4) + facet_grid(CHROM~.)

#identifying the minimum number of differences between any strains comparison 
divergence[, list("min_divergence"=min(N, na.rm=TRUE)), by=list(CHROM,start)]

#to count number of locations in chromosome that have heterozygous genotype
het1 <- sub1[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het2 <- sub2[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het3 <- sub3[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het4 <- sub4[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het5 <- sub5[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het6 <- sub6[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het7 <- sub7[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het8 <- sub8[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het9 <- sub9[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))),  .SDcols=strain_names]
het10 <- sub10[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het11 <- sub11[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het12 <- sub12[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het13 <- sub13[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het14 <- sub14[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het15 <- sub15[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het16 <- sub16[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]

#make all het files data tables
het1  <- as.data.table(t(het1))
het2  <- as.data.table(t(het2))
het3  <- as.data.table(t(het3))
het4  <- as.data.table(t(het4))
het5  <- as.data.table(t(het5))
het6  <- as.data.table(t(het6))
het7  <- as.data.table(t(het7))
het8  <- as.data.table(t(het8))
het9  <- as.data.table(t(het9))
het10 <- as.data.table(t(het10))
het11 <- as.data.table(t(het11))
het12 <- as.data.table(t(het12))
het13 <- as.data.table(t(het13))
het14 <- as.data.table(t(het14))
het15 <- sub15[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]
het16 <- sub16[, apply(.SD, 2, function(x) sum(x %in% c("0/1","1/0"))), .SDcols=strain_names]

#make all het files data tables
het1  <- as.data.table(t(het1))
het2  <- as.data.table(t(het2))
het3  <- as.data.table(t(het3))
het4  <- as.data.table(t(het4))
het5  <- as.data.table(t(het5))
het6  <- as.data.table(t(het6))
het7  <- as.data.table(t(het7))
het8  <- as.data.table(t(het8))
het9  <- as.data.table(t(het9))
het10 <- as.data.table(t(het10))
het11 <- as.data.table(t(het11))
het12 <- as.data.table(t(het12))
het13 <- as.data.table(t(het13))
het14 <- as.data.table(t(het14))
het15 <- as.data.table(t(het15))
het16 <- as.data.table(t(het16))

#collect counts and combine
het1[, chromosome := 1 ]
het2[, chromosome := 2 ]
het3[, chromosome := 3 ]
het4[, chromosome := 4 ]
het5[, chromosome := 5 ]
het6[, chromosome := 6 ]
het7[, chromosome := 7 ]
het8[, chromosome := 8 ]
het9[, chromosome := 9 ]
het10[, chromosome := 10]
het11[, chromosome := 11]
het12[, chromosome := 12]
het13[, chromosome := 13]
het14[, chromosome := 14]
het15[, chromosome := 15]
het16[, chromosome := 16]
all_het <- rbindlist(list(
  het1,
  het2,
  het3,
  het4,
  het5,
  het6,
  het7,
  het8,
  het9,
  het10,
  het11,
  het12,
  het13,
  het14,
  het15,
  het16)
)

fwrite(all_het, file = "all_het.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")