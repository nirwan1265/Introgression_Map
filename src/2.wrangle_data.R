# Set the directory path to your data folder
data_folder <- paste0(getwd(), "/data")

# Create the "data/filtered" directory if it doesn't exist
filtered_folder <- file.path(data_folder, "filtered")
if (!dir.exists(filtered_folder)) {
  dir.create(filtered_folder)
}

# List all .tsv files in the data folder with full paths
tsv_files <- list.files(path = data_folder, pattern = "\\.tsv$", full.names = TRUE)

# Define the valid chromosome names you want to keep
valid_chromosomes <- paste0("chr", 1:3)

# Loop through each .tsv file, remove header lines, reorder columns, remove the first row, and save as .txt separately
for (file_path in tsv_files) {
  # Read the .tsv file using read.delim() without row names and skip lines starting with "@"
  data <- read.delim(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL, comment.char = "@")
  
  # Remove header lines (lines starting with "@") from the data frame
  #data <- data[!grepl("^@", data[, 1]), ]
  
  # Reorder columns
  #data <- data[, c(1, 2, 5, 3, 6, 4)]
  
  # Remove the first row (which might contain the column names)
  #data <- data[-1, ]
  
  # Filter rows to keep only those with valid chromosome names
  data <- data[data[, 1] %in% valid_chromosomes, ]
  
  # Get the original file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Define the path to save the modified file as .txt in the "data/filtered" directory
  save_path <- file.path(filtered_folder, paste0(file_name, ".txt"))
  
  # Save the modified data as a tab-delimited text file with no colnames or rownames without quotes
  write.table(data, file = save_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Print a message indicating the file has been processed and saved
  cat(paste("Processed and saved:", save_path), "\n")
}
str(data)
rm(data)
