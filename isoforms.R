source("logging.R")

required_columns <- c("mouse", "GenLit", "Sex", "DietName")
option_column <- c("Batch")

# --- Check isoform file ---
check_isoforms <- function(isoforms, trait) {
  # Check to make sure "mouse, GenLit, Sex, DietName" are in the CSV file as column names
  missing_columns <- setdiff(required_columns, colnames(isoforms))

  if (length(missing_columns) > 0) {
    stop_msg("The following required columns are missing:", paste(missing_columns, collapse = ", "))
  } else {
    log_msg("All required columns are present.")
  }


  # Check to make sure the trait is in the CSV file as a column name
  if (!(trait %in% colnames(isoforms))) {
    stop_msg("The trait", trait, "is not present in the CSV file.")
  } else {
    log_msg("Trait", trait, "is present in the CSV file.")
  }
}

# --- Get possible traits ---
get_possible_traits <- function(csv) {
  # Get possible traits
  possible_traits <- setdiff(colnames(csv), c(required_columns, option_column))
  log_msg("Possible traits:", paste(possible_traits, collapse = ", "))
  return(possible_traits)
}

# --- Remove unused columns ---
remove_unused_columns <- function(csv, trait) {
  #log columns to be removed
  unused_columns <- setdiff(colnames(csv), c(required_columns, option_column, trait))
  log_msg("Removing unused columns:", paste(unused_columns, collapse = ", "))
  # Remove unused columns
  csv <- csv[, c(required_columns, option_column, trait)]
  log_msg("Columns after removal:", paste(colnames(csv), collapse = ", "))
  return(csv)
}

# --- Filter isoforms ---
filter_isoforms <- function(isoforms, chromosome, position, flank, minexp) {
  log_msg("Isoforms before filtering for position:", nrow(isoforms))
  # Filter isoforms based on chromosome, position, and flank
  #
  isoforms <- subset(isoforms,
                            chromosome_name == chromosome &
                            transcript_start >= position - flank &
                            transcript_end <= position + flank)
  log_msg("Isoforms after filtering for position:", nrow(isoforms))

  # Extract columns starting with 'DO'
  DO_cols_logical <- grepl("^DO", names(isoforms))
  DO_columns <- isoforms[, DO_cols_logical]

  # Identify the smallest number in columns starting with 'DO'
  if (ncol(DO_columns) == 0) {
    stop_msg("No columns starting with 'DO' found.")
  }

  smallest <- min(apply(DO_columns, 2, min, na.rm = TRUE), na.rm = TRUE)

  log_msg("Smallest value in DO columns:", smallest)
  log_msg("Isoforms before filtering for minimum expression:", nrow(isoforms))
  # Filter isoforms based on minimum non-zero expression value
  isoforms <- isoforms[which(rowSums(DO_columns > smallest, na.rm = TRUE) >= minexp),]
  log_msg("Isoforms after filtering for minimum expression:", nrow(isoforms))
  return(isoforms)
}

# --- Filter isoforms for trait ---
filter_isoforms_for_trait <- function(isoforms, trait) {
  log_msg("Isoforms before filtering for trait:", nrow(isoforms))
  # Filter isoforms for the 'external_gene_name' that contains the trait
  isoforms <- subset(isoforms, grepl(trait, external_gene_name))
  log_msg("Isoforms after filtering for trait:", nrow(isoforms))
  return(isoforms)
}

# --- Format isoforms ---
format_isoforms <- function(isoforms) {
  # Create a new column 'unique_name' combining 'X' and 'external_gene_name'
  unique_name <- paste(isoforms$X, isoforms$external_gene_name, sep = ".")

  # Remove all non-DO columns
  DO_cols_logical <- grepl("^DO", names(isoforms))
  DO_column_names <- names(isoforms)[DO_cols_logical]
  isoforms <- isoforms[, DO_column_names]

  #Copy the header to a new object
  header <- names(isoforms)

  # Transpose data and convert to a data frame
  isoforms <- as.data.frame(t(isoforms))

  # Set 'unique_name' as the header
  colnames(isoforms) <- unique_name

  # Add 'Mouse' column and reorganize data
  isoforms$mouse <- header
  isoforms <- isoforms[, c(ncol(isoforms), 1:(ncol(isoforms) - 1))]

  # Sort data by the 'Mouse' column
  isoforms <- isoforms[order(isoforms$mouse),]
}

# --- Combine files ---
combine_files <- function(csv, isoforms) {
  log_msg("Combining CSV and isoforms files.")
  mouse <- unique(isoforms$mouse)
  mouse <- mouse[!is.na(mouse)]
  mouse <- mouse[!mouse %in% csv$mouse]
  if (length(mouse) > 0) {
    log_msg("The following mouse are not in the CSV file:", paste(mouse, collapse = ", "))
    isoforms <- isoforms[!isoforms$mouse %in% mouse,]
  }
  mouse <- unique(csv$mouse)
  mouse <- mouse[!is.na(mouse)]
  mouse <- mouse[!mouse %in% isoforms$mouse]
  if (length(mouse) > 0) {
    log_msg("The following mouse are not in the isoforms file:", paste(mouse, collapse = ", "))
    csv <- csv[!csv$mouse %in% mouse,]
  }
  log_msg("Columns after combining:", ncol(isoforms))
  #Combine csv and isoforms by mouse
  return(merge(csv, isoforms, by = "mouse"))
}

