# A function to read in data from a folder where all the files are csv
# Author: Sam Johnson

read_in_data <- function(file_path) {
  temp_files <- Sys.glob(file_path)
  temp_metrics <- sapply(temp_files, read.csv, simplify = FALSE)
  temp_metrics <- do.call(rbind, temp_metrics)
  temp_metrics
}
