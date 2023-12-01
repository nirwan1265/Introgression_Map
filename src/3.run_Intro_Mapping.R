# Get paths to allele count files
file_paths = list.files(path = paste0(getwd(),"/data/filtered"), pattern = "\\.txt$", full.names=T)

# Get sample names
sampleIDs <- basename(file_paths)
# Remove the .txt extension and "allelicCounts" part from sampleIDs
sampleIDs <- gsub("\\.txt$", "", sampleIDs)
sampleIDs <- gsub("\\_AD_processed$", "", sampleIDs)

# Create the expDesign object
expDesign = data.frame(files=file_paths, name=sampleIDs)


# Get chromosome lengths for the example data included in the package
# chromosome length from the alleliccount file
# Create the "SN" and "chr" names
#chr_len <- c(308452471,243675191,238017767,250330460,226353449,181357234,185808916,182411202,163004744,152435371)
chr_len <- c(308452471,243675191,238017767)
#names(chr_len) <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
names(chr_len) <- c('chr1','chr2','chr3')

# Running the RTIGER function that identifies Introgression regions
output <- paste0(getwd(),"/result")
myres2 = RTIGER(expDesign = expDesign,
               outputdir = output,
               seqlengths = chr_len,
               rigidity = 700,
               autotune = FALSE,
               average_coverage = 0.8,
               nstates = 3,
               crossovers_per_megabase = 1, 
               post.processing = FALSE,
               save.results = TRUE)
RTIGER::calcCOnumber(myres)
RTIGER::plotCOs(myres)
RTIGER::optimize_R(myres)
post_post.processing=NULL
?RTIGER()
  