install.packages("rCNV")
library(rCNV)


vcf <- readVCF("~/Desktop/genotype_diplo_chr10_hardfilter_bcffiltered.vcf.gz")

pp<-substr(colnames(vcf$vcf)[-c(1:9)],1,8)




?maf()
?hetTgen()

h.table <- hetTgen(
  vcf,
  info.type = c("AD", "AD-tot", "GT", "GT-012", "GT-AB", "DP"),
  verbose = TRUE
)

h.table <- maf(h.table, AD = TRUE, verbose = TRUE)
pp<-substr(colnames(h.tabl)[-c(1:9)],1,8)
hzygots<-h.zygosity(h.table,plot=TRUE,pops=pp)
