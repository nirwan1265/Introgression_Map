# Get paths to allele count files
file_paths = list.files(path = paste0(getwd(),"/data/filtered"), pattern = "\\.txt$", full.names=T)

# Get sample names
sampleIDs <- basename(file_paths)
# Remove the .txt extension and "allelicCounts" part from sampleIDs
sampleIDs <- gsub("\\.txt$", "", sampleIDs)
sampleIDs <- gsub("\\.allelicCounts$", "", sampleIDs)

# Create the expDesign object
expDesign = data.frame(files=file_paths, name=sampleIDs)


# Get chromosome lengths for the example data included in the package
# chromosome length from the alleliccount file
# Create the "SN" and "chr" names
#chr_len <- c(308452471,243675191,238017767,250330460,226353449,181357234,185808916,182411202,163004744,152435371)
chr_len <- c(182411202,163004744,152435371)
#names(chr_len) <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
names(chr_len) <- c('chr8','chr9','chr10')
152435371*0.00000662
# Running the RTIGER function that identifies Introgression regions
output <- paste0(getwd(),"/result")
myres = RTIGER(expDesign = expDesign,
               outputdir = output,
               seqlengths = chr_len,
               rigidity = 20,
               autotune = TRUE,
               average_coverage = 0.8,
               nstates = 3,
               crossovers_per_megabase = 3, 
               post.processing = TRUE,
               save.results = TRUE)


