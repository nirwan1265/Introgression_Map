# Load Packages
library(rCNV)
library(vcfR)

# Directory containing VCF files
dir <- "/Users/nirwantandukar/Documents/BZea_files/diplo"

# Function to process each VCF file
process_vcf <- function(file_path) {
  # Read VCF file
  vcf <- readVCF(file_path)
  
  # Sample Names
  pp <- substr(colnames(vcf$vcf)[-c(1:9)], 1, 8)
  
  # Calculate Homozygosity and Heterozygosity
  homo_hzygots <- h.zygosity(vcf, plot=TRUE, pops=pp)
  homo_hzygots$heterozygotes <- homo_hzygots$total - homo_hzygots$`O(Hom)`
  homo_hzygots$heterozygosity <- homo_hzygots$heterozygotes/homo_hzygots$total
  homo_hzygots <- homo_hzygots[,-6]
  homo_hzygots$chromosome <- sub(".*_chr(\\d+).*", "\\1", file_path)
  
  # Read VCF for allele frequency
  vcf2 <- read.vcfR(file_path)
  
  # Extracting Genotype Information
  geno <- extract.gt(vcf2, element = "GT")
  
  # Calculate the percentage of NA in each row
  na_percentage <- apply(geno, 1, function(x) {
    na_count <- sum(is.na(x))
    total_count <- length(x)
    return((na_count / total_count) * 100)
  })
  
  # Function to count genotypes
  count_genotypes <- function(column) {
    homo_ref_count <- sum(column == "0/0" | column == "0|0", na.rm = TRUE)
    homo_alt_count <- sum(column == "1/1" | column == "1|1", na.rm = TRUE)
    hetero_count <- sum(column == "0/1" | column == "1/0", na.rm = TRUE)
    return(c(homo_ref = homo_ref_count, homo_alt = homo_alt_count, hetero = hetero_count))
  }
  
  # Apply the function to each sample (column)
  genotype_counts <- apply(geno, 2, count_genotypes)
  
  # Convert to a data frame
  genotype_counts_df <- as.data.frame(t(genotype_counts))
  
  # Naming the rows after the samples
  rownames(genotype_counts_df) <- colnames(geno)
  genotype_counts_df$total <- rowSums(genotype_counts_df)
  
  # Calculate the percentages and add them to the data frame
  genotype_counts_df$homo_ref_perc <- (genotype_counts_df$homo_ref / genotype_counts_df$total) * 100
  genotype_counts_df$homo_alt_perc <- (genotype_counts_df$homo_alt / genotype_counts_df$total) * 100
  genotype_counts_df$hetero_perc <- (genotype_counts_df$hetero / genotype_counts_df$total) * 100
  
  # Extract REF and ALT alleles
  ref_alleles <- vcf2@fix[, "REF"]
  alt_alleles <- vcf2@fix[, "ALT"]
  
  # Extract Allele Depth (AD)
  ad_info <- extract.gt(vcf2, element = "AD")
  
  # Combine the extracted information into a data frame
  allele_info_df <- data.frame(REF = ref_alleles, ALT = alt_alleles, AD = ad_info)
  
  # Return the result as a list
  return(list(na_percentage=na_percentage, homo_hzygots=homo_hzygots, genotype_counts_df=genotype_counts_df, allele_info_df=allele_info_df))
}

# List all VCF files
vcf_files <- list.files(dir, pattern="genotype_diplo_chr\\d+_hardfilter_bcffiltered_bialleles\\.vcf\\.gz", full.names=TRUE)

# Process each file and store results in a list
results <- lapply(setNames(vcf_files, gsub(".*_(chr\\d+)_.*", "diplo_\\1", vcf_files)), process_vcf)

# Set up plot layout and file
png(file="missing_data_histograms.png", width=2400, height=1200)
par(mfrow=c(2, 5))

# Plot the histogram for each file
names(results) <- gsub(".*_(chr\\d+)_.*", "\\1", names(results)) # Adjusting names to be chr1, chr2, etc.
for(i in seq_along(results)) {
  hist(results[[i]]$na_percentage, main=paste("Histogram of Missing Data for", names(results)[i]), 
       xlab="Percentage of Missing Data", 
       ylab="Frequency", 
       col="blue", 
       breaks=10)
}
# Close the PNG device
dev.off()


