

dat <- read.csv("199_phenotypes.csv")

uniq_ids <- unique(dat[["ecotype_id"]])
uniq_phens <- as.character(unique(dat[["phenotype_name"]]))

pheno_mat <- matrix(data=NA, nrow=length(uniq_ids), 
                             ncol=(1+ length(uniq_phens)))
pheno_df <- data.frame(pheno_mat)


# add names
colnames(pheno_df) <- c("id", uniq_phens)
pheno_df[["id"]] <-  uniq_ids

for(ii in uniq_phens){
  print(ii)
  rindx <- which(dat[["phenotype_name"]] == ii)
  sub_data <- dat[rindx, ]

  cindx_p <- which(colnames(pheno_df) == ii)
  rindx_p <- match(sub_data[["ecotype_id"]], pheno_df[, "id"])

  pheno_df[rindx_p, cindx_p] <- sub_data[["value"]]

}

print(head(pheno_df))

write.table(pheno_df, file="../pheno.txt")

