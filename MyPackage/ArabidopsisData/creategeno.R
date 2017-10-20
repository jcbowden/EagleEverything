## purpose: to put Aribadopsis genogtype data into a raw form WMAM can read
DIR <- "/home/geo047/gitHUB_WMAM/MyPackage/ArabidopsisData/"

setwd(DIR)

## read data
## contains genotypes for 1000+ plants
cat(" Reading in geno data ... \n")
dat <- read.csv("all_chromosomes_binary.csv")

cat(" Processing ... \n")
## read pheno file so that id's can be matched 
pheno <- read.table("./pheno.txt", header=TRUE)

## match pheno id's with id's in genotype file (remembering columns are id's in genotype file)
keep.cindx <- match(pheno[["id"]], as.numeric(substring(names(dat)[-(1:2)], first=2)) )  + 2
subdat <- dat[, keep.cindx]



## take transpose of subdat so that plants become the rows
subdat <- t(subdat)


write.table(subdat , file="./geno.txt", row.names=FALSE, col.names=FALSE)


