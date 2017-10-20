## read SNP data
cat("reading in genoytpe data ... \n")
geno <- read.table("geno.txt")  ## Arabidopsis genotype data 

## take 100 x 200 block of genotypes
cat(" Forming data sets ... \n")
set.seed(101)
geno_small <- geno[1:100, sample(1:ncol(geno), 200, FALSE)]


## create generic map
library(mpMap)
map <- sim.map(len=rep(100,10), n.mar=20, eq.spacing=FALSE, include.x=FALSE)
df <- data.frame(Mrk=NULL, Chr=NULL, Pos=NULL)
for(ii in 1:10)
{
  chrm.df <- data.frame(Mrk=names(map[[ii]]), Chr=ii, Pos=as.vector(map[[ii]]))
  df <- rbind.data.frame(df, chrm.df)

}



## simulate phenotypic data 
## id, cov1, cov2, fac, trait1, trait2, trait3

cov1 <- rnorm(100, mean=10, sd=2)
cov2 <- runif(100, 0, 5)
fac <- as.factor(sample(c("M","F"), 100, TRUE))

trait1 <- 2*cov1 +  2*cov2 + 3*as.numeric(fac) + 5.0*geno_small[,100] + rnorm(100, 0, 1)

logp1 <- rep(NA, 200)
logp2 <- rep(NA, 200)
logp3 <- rep(NA, 200)
for(ii in 1:ncol(geno_small)){
    mrk <- geno_small[, ii]
    mdl <- lm(formula(trait1 ~ cov1 + cov2 + fac + mrk) )
    mat <- summary(mdl)$coefficients
    rindx <- which(rownames(mat)=="mrk")
     logp1[ii] <- -log( mat[rindx, 4])
}
plot(logp1)
max(logp1)


trait2 <- 2*cov1 +  5*geno_small[,100] + rnorm(100, 0, 1)
logp2 <- rep(NA, 200)
for(ii in 1:ncol(geno_small)){
    mrk <- geno_small[, ii]
    mdl <- lm(formula(trait2 ~  cov2 +  mrk) )
    mat <- summary(mdl)$coefficients
    rindx <- which(rownames(mat)=="mrk")
     logp2[ii] <- -log( mat[rindx, 4])
}   
plot(logp2)
max(logp2)

logp3 <- rep(NA, 200)
trait3 <- 0.8 * geno_small[,100] + rnorm(100, 0, 1)
for(ii in 1:ncol(geno_small)){
    mrk <- geno_small[, ii]
    mdl <- lm(formula(trait3 ~ mrk) )
    mat <- summary(mdl)$coefficients
    rindx <- which(rownames(mat)=="mrk")
     logp3[ii] <- -log( mat[rindx, 4])
}   
plot(logp3)
max(logp3)

### create data frame and save geno and pheno files

write.table(geno_small, file = "../genoexampleCwise.txt", row.names=FALSE, col.names=FALSE)
write.table(t(geno_small), file = "../genoexampleRwise.txt", row.names=FALSE, col.names=FALSE)

pheno <-  data.frame(cov1=round(cov1,3), cov2=round(cov2,3), fac=fac,  trait1=round(trait1,3), 
                                                     trait2=round(trait2,3), 
                                                     trait3=round(trait3,3))

write.table(pheno, file="../phenoexample.csv", row.names=FALSE, col.names=TRUE, sep=",")


write.table(df, file="../mapexample.txt", row.names=FALSE, col.names=TRUE)



