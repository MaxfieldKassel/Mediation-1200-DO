
library(data.table)
library(readxl)
source("logging.R")

# --- Load CSV File ---
load_csv <- function(filename) {
  log_msg("Loading CSV file:", filename)
  #If the file is a CSV, use read.csv
  if (grepl(".csv$", filename)) {
    return(as.data.frame(fread(filename, showProgress = FALSE, blank.lines.skip = TRUE)))
  } else if (grepl(".xlsx$", filename)) {
    return(as.data.frame(read_xlsx(filename)))
  } else {
    stop_msg("File type not supported.")
  }
}

#--- Create folder --- 
create_folder <- function(folder, add_timestamp = TRUE) {
  if (add_timestamp) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    folder <- paste0(folder, "_", timestamp)
  }
  if (!dir.exists(folder)) {
    dir.create(folder)
  }
  return(folder)
}