library(data.table)

# Read in the CSV file
csv <- read.csv(file = "20230627_DOplasmalipids.csv")

# Check to make sure "mouse, GenLit, Sex, DietName" are in the CSV file as column names
required_columns <- c("mouse", "GenLit", "Sex", "DietName")
missing_columns <- setdiff(required_columns, colnames(csv))

if (length(missing_columns) > 0) {
  print(paste("The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  #Exit
  q()
} else {
  print("All required columns are present.")
}

#Get other columns from csv
other_columns <- setdiff(colnames(csv), required_columns)
#Prompt user for column to use as trait
trait_name <- readline(prompt = paste("Enter the column name to use as the trait (", paste(other_columns, collapse = ", "), "): "))
#Check to make sure the trait column is in the CSV file
if (!trait_name %in% colnames(csv)) {
  print(paste("The trait column", trait_name, "is not in the CSV file."))
  #Exit
  q()
} else {
  print("The trait column is present.")
}

#Remove all other columns from csv (except for required, selected and "Batch" if present)
batch_present <- "Batch" %in% colnames(csv)
if (batch_present) {
  csv <- csv[, c(required_columns, "Batch", trait_name)]
} else {
  csv <- csv[, c(required_columns, trait_name)]
}

csv <- as.data.table(csv)
print("The following columns will be used:")
print(colnames(csv))
#print header




# Define a function to rank data
rankz = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

# Set parameters
chromosome <- "5"
position <- 125000000
flanking <- 2000000
minexp <- 100
percent <- 0.5

# Load isoforms data
Isoforms <- fread("vsd_DO_isoforms_new_app.csv")

# Extract the isoforms for each region based on user input (chromosome, position, flanking region)
# Filter: If the chromosome matches and the position is within the flanking region, keep the row
isoforms1 <- Isoforms[chromosome_name == chromosome & transcript_start >= position - flanking & transcript_end <= position + flanking]

# Clean up by removing the original isoforms object
rm(Isoforms)

# Extract columns starting with 'DO'
DO_cols_logical <- grepl("^DO", names(isoforms1))
DO_columns <- isoforms1[, ..DO_cols_logical]

# Identify the smallest number in columns starting with 'DO'
if (ncol(DO_columns) > 0) {
  smallest <- min(apply(DO_columns, 2, min, na.rm = TRUE), na.rm = TRUE)
} else {
  smallest <- NA # Placeholder for no 'DO' columns
}

# Filter isoforms based on minimum non-zero expression value
isoforms1 <- isoforms1[which(rowSums(DO_columns > smallest, na.rm = TRUE) >= minexp),]

# Create a new column 'unique_name' combining 'X' and 'external_gene_name'
unique_name <- paste(isoforms1$X, isoforms1$external_gene_name, sep = ".")

# Remove all non-DO columns
DO_column_names <- names(isoforms1)[DO_cols_logical]
isoforms1 <- isoforms1[, ..DO_column_names]

#Copy the header to a new object
header <- names(isoforms1)

# Transpose data and convert to a data frame
isoforms1 <- as.data.frame(t(isoforms1))

# Set 'unique_name' as the header
colnames(isoforms1) <- unique_name

# Add 'Mouse' column and reorganize data
isoforms1$mouse <- header
isoforms1 <- isoforms1[, c(ncol(isoforms1), 1:(ncol(isoforms1) - 1))]

# Sort data by the 'Mouse' column
isoforms1 <- isoforms1[order(isoforms1$mouse),]

# Write the transformed data to a CSV file
fwrite(isoforms1, file = "isoforms1.csv")

mouse <- unique(isoforms1$mouse)
mouse <- mouse[!is.na(mouse)]
mouse <- mouse[!mouse %in% csv$mouse]
if (length(mouse) > 0) {
  print(paste("The following mouse are not in the CSV file:", paste(mouse, collapse = ", ")))
  isoforms1 <- isoforms1[!isoforms1$mouse %in% mouse,]
}
mouse <- unique(csv$mouse)
mouse <- mouse[!is.na(mouse)]
mouse <- mouse[!mouse %in% isoforms1$mouse]
if (length(mouse) > 0) {
  print(paste("The following mouse are not in the isoforms1 file:", paste(mouse, collapse = ", ")))
  csv <- csv[!csv$mouse %in% mouse,]
}
#Combine csv and isoforms1 by mouse
combined <- merge(csv, isoforms1, by = "mouse")

#Save combined as csv
fwrite(combined, file = "combined.csv")





# --- Libraries ---
library(gplots)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(qtl2)
library(qtl2convert)
library(parallel)

# --- Configuration Variables ---
filename <- "combined.csv"

generate_plots <- TRUE
x_limit <- c(0, 120)
y_limit <- c(0, 13)

folder <- "ResultsFolder"
add_timestamp <- TRUE # Set to TRUE to add a timestamp, or FALSE to exclude it

numCores <- detectCores() - 1 # Keeping one core free for other processes

# --- Load Data ---
genoprobs <- readRDS("final_genoprobs_1176.rds")
map <- readRDS("gigamuga_map_test_v2.RDS")
K <- readRDS("K_1176_DietDO.rds")
pheno <- as.data.frame(combined)
colnames(pheno)[1] <- "mouse"
rownames(pheno) <- pheno$mouse
cat("Data loaded.\n")

# --- Create Folder ---
if (add_timestamp) {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  folder <- paste0(folder, "_", timestamp)
}
if (!dir.exists(folder)) {
  dir.create(folder)
}
cat("Folder for results created:", folder, "\n")

# --- Get Gene Names ---
gene_names <- unique_name

# --- Make a kinship matrix and subset for the specified chromosome ---
K_chr <- K[[chromosome]]
gp <- genoprobs[, chromosome]
cat("Kinship matrix and genoprobs subsetted.\n")

# --- Reclassify covars ---
pheno$Diet <- as.factor(pheno$Diet)
pheno$GenLit <- as.factor(pheno$GenLit)
pheno$Sex <- as.factor(pheno$Sex)
base_covar <- model.matrix(~Sex * Diet + GenLit, data = pheno)[, -1]
interactive_covar <- t(t(model.matrix(~Diet, data = pheno)[, -1]))
cat("Covariates reclassified.\n")


#------ check LOD profile for trait without a transcript as an additive covar -----------
pheno[[paste0("rankz.", trait_name)]] <- rankz(pheno[[trait_name]])
qtl_no_additive <- scan1(genoprobs = gp, pheno = pheno[, paste0("rankz.", trait_name), drop = FALSE], kinship = K_chr, addcovar = base_covar, intcovar = interactive_covar)
cat("QTL scan for trait alone completed.\n")

# Plot and save to PNG
if (generate_plots) {
  png_filename <- paste0(folder, "/", trait_name, "_alone.png")
  png(png_filename)
  plot_scan1(x = qtl_no_additive, xlim = x_limit, ylim = y_limit, map = map, main = "QTL for trait alone")
  dev.off()
}
qtl_no_additive[qtl_no_additive == 0] <- 1e-6 # Avoid division by zero


# --- Function to calculate LOD difference and perform original QTL analysis ---
analyze_gene <- function(gene, genoprobs, pheno, K_chr, base_covar, interactive_covar, map, qtl_no_additive) {
  # Transform data
  pheno[[paste0("rankz.", gene)]] <- rankz(pheno[[gene]])

  # Perform scan with additive covariate
  addcovar_with_gene <- model.matrix(~Sex * Diet + GenLit + pheno[[paste0("rankz.", gene)]], data = pheno)[, -1]
  qtl_with_additive <- scan1(genoprobs = genoprobs, pheno = pheno[, trait_name, drop = FALSE], kinship = K_chr, addcovar = addcovar_with_gene, intcovar = interactive_covar, cores = 0)

  # Calculate percent change in LOD scores and handle division by zero
  lod_diff <- qtl_with_additive - qtl_no_additive
  percent_change_lod <- lod_diff / qtl_no_additive * 100

  # Find the index and value of the maximum percent change
  max_index <- which(abs(percent_change_lod) == max(abs(percent_change_lod), na.rm = TRUE))
  max_percent_change <- percent_change_lod[max_index]

  # Extract chromosome-specific positions from map
  chromosome_positions <- map[[chromosome]]

  # Validate max_index and extract max_change_position
  if (length(max_index) > 0 && max_index <= length(chromosome_positions)) {
    max_change_position <- chromosome_positions[max_index]
  } else {
    max_change_position <- NA # Assign NA if max_index is not valid
  }

  # Find peaks for the additive model
  peaks <- find_peaks(scan1_output = qtl_with_additive, map = map, threshold = 9)

  # Plot and save to PNG
  if (generate_plots) {
    png_filename <- paste0(folder, "/", gene, ".png")
    png(png_filename)
    plot_scan1(x = qtl_with_additive, xlim = x_limit, ylim = y_limit, map = map, main = paste("QTL Scan for ", trait_name, " + ", gene))
    dev.off()
  }

  # Add gene name, max percent change, and its location to peak results
  peaks$gene <- gene
  peaks$max_percent_change <- max_percent_change
  peaks$max_change_position <- max_change_position

  return(peaks)
}

# --- Set up Cluster ---
cl <- makeCluster(numCores)
clusterExport(cl, c("pheno", "gp", "K_chr", "base_covar", "interactive_covar", "rankz", "map", "find_peaks", "scan1",
 "analyze_gene", "trait_name", "qtl_no_additive", "generate_plots", "plot_scan1", "folder", "x_limit", "y_limit", "chromosome"))
cat("Cluster set up and functions exported.\n")

# --- Parallel Loop for Gene Analysis ---
analysis_results <- parLapply(cl, gene_names, function(gene) {
  analyze_gene(gene, gp, pheno, K_chr, base_covar, interactive_covar, map, qtl_no_additive)
})

# --- Close Cluster ---
stopCluster(cl)

# --- Aggregate Results ---
analysis_df <- do.call(rbind, analysis_results)

# --- Remove max_diff and max_diff_loc from analysis_results before saving ---
analysis_df_no_max_change <- analysis_df[, !(names(analysis_df) %in% c("max_percent_change", "max_change_position"))]
write.csv(analysis_df_no_max_change, file.path(folder, paste0(trait_name, "_analysis_results.csv")), row.names = FALSE)

# --- Extract and Save LOD Difference Results with Percent Change ---
lod_change_df <- analysis_df[, c("gene", "max_percent_change", "max_change_position")]
lod_change_df <- lod_change_df[order(-lod_change_df$max_percent_change),] # Sorting by largest percent change
write.csv(lod_change_df, file.path(folder, paste0(trait_name, "_LOD_percent_changes.csv")), row.names = FALSE)