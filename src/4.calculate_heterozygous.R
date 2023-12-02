install.packages("rCNV")
library(rCNV)
library(vcfR)

vcf <- readVCF("~/Downloads/genotype_diplo_chr1_hardfilter_bcffiltered_bialleles.vcf.gz")

pp<-substr(colnames(vcf$vcf)[-c(1:9)],1,8)
hzygots <- h.zygosity(vcf,plot=TRUE,pops=pp)





str(hzygots)

hzygots$heterozygotes <- hzygots$total - hzygots$`O(Hom)`
hzygots$heterozygosity <- hzygots$heterozygotes/hzygots$total


vcf <- readVCF("~/Downloads/genotype_diplo_chr1_hardfilter_bcffiltered_bialleles.vcf.gz")

vcf2 <- read.vcfR("~/Downloads/genotype_diplo_chr1_hardfilter_bcffiltered_bialleles.vcf.gz")
geno <- extract.gt(vcf2, element = "GT")

str(geno)

# Calculate the percentage of NA in each row
na_percentage <- apply(geno, 1, function(x) {
  na_count <- sum(is.na(x))
  total_count <- length(x)
  return((na_count / total_count) * 100)
})

# Plot a histogram of the NA percentages
hist(na_percentage, main = "Histogram of Missing Data Percentage", 
     xlab = "Percentage of Missing Data", 
     ylab = "Frequency", 
     col = "blue", 
     breaks = 10) # Adjust 'breaks' as needed for bin size


# Function to count genotypes
count_genotypes <- function(column) {
  homo_ref_count <- sum(column == "0/0", na.rm = TRUE)
  homo_alt_count <- sum(column == "1/1", na.rm = TRUE)
  hetero_count <- sum(column == "0/1" | column == "1/0", na.rm = TRUE)
  return(c(homo_ref = homo_ref_count, homo_alt = homo_alt_count, hetero = hetero_count))
}

# Apply the function to each sample (column)
genotype_counts <- apply(geno, 2, count_genotypes)

# Convert to a data frame
genotype_counts_df <- as.data.frame(t(genotype_counts))

# Naming the rows after the samples
rownames(genotype_counts_df) <- colnames(geno)

# Display the data frame
genotype_counts_df


# Calculate the percentages and add them to the data frame
genotype_counts_df$homo_ref_perc <- (genotype_counts_df$homo_ref / rowSums(genotype_counts_df)) * 100
genotype_counts_df$homo_alt_perc <- (genotype_counts_df$homo_alt / rowSums(genotype_counts_df)) * 100
genotype_counts_df$hetero_perc <- (genotype_counts_df$hetero / rowSums(genotype_counts_df)) * 100

# Display the updated data frame
genotype_counts_df



# Extract REF and ALT alleles
ref_alleles <- vcf2@fix[, "REF"]
alt_alleles <- vcf2@fix[, "ALT"]

# Extract Allele Depth (AD)
# This part can vary depending on how the AD is stored in your VCF
# Assuming AD is the first field in the genotype (GT:AD:DP:GQ:PL format)
ad_info <- extract.gt(vcf2, element = "AD")

# Combine the extracted information into a data frame
allele_info_df <- data.frame(REF = ref_alleles, ALT = alt_alleles, AD = ad_info)
allele_info_df[1:5,1:10]
