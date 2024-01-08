# Load data table library
source("logging.R")
source("utils.R")
source("isoforms.R")
# --- Libraries ---
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(qtl2)
library(qtl2convert)
library(parallel)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the number of arguments is correct
expected_num_args <- 8
if (length(args) != expected_num_args) {
  stop("Incorrect number of arguments. Expected ", expected_num_args, ", got ", length(args), ".")
}

# Assign arguments to variables
csv_filename <- args[1]
trait_name <- args[2]
chromosome <- args[3]
position <- as.numeric(args[4])
flanking <- as.numeric(args[5])
minexp <- as.numeric(args[6])
percent <- as.numeric(args[7])
threshold <- as.numeric(args[8])

# Define a function to rank data
rankz = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

numCores <- min(100, detectCores() - 1) # Keeping one core free for other processes


# Non-interactive parameters
output_folder <- "Output/"
input_folder <- "Input/"
data_folder <- "Data/"
folder <- paste0(output_folder, trait_name, "_", chromosome, "_", position) # Folder to save results to
add_timestamp <- TRUE # Set  to TRUE to add a timestamp, or FALSE to exclude it
isoforms_filename <- paste0(data_folder, "vsd_DO_isoforms_new_app.csv")
gp_filename <- paste0(data_folder, "Genoprobs/chromosome_")
map_filename <- paste0(data_folder, "gigamuga_map_test_v2.RDS")
k_filename <- paste0(data_folder, "K_1176_DietDO.rds")
csv_filename <- paste0(input_folder, csv_filename)

# Create folder and set log file
folder <- create_folder(folder, add_timestamp)
set_log_file(folder)
log_msg("Folder for results created:", folder)
log_msg("Log file set.")
log_msg("File name:", csv_filename)
log_msg("Trait name:", trait_name)
log_msg("Chromosome:", chromosome)
log_msg("Position:", position, "mbp")
position <- position * 1000000 # Convert to base pairs
log_msg("Flanking region:", flanking, "mbp")
flanking <- flanking * 1000000 # Convert to base pairs
log_msg("Minimum expression value:", minexp)
log_msg("Percent LOD change:", percent)
log_msg("Threshold:", threshold)
log_msg("Number of cores:", numCores)

log_msg("Loading CSV file...")



# Read in the CSV file
csv <- load_csv(csv_filename)
# Check the CSV file for required columns, stop if any are missing
check_isoforms(csv, trait_name)
# Remove unused columns
csv <- remove_unused_columns(csv, trait_name)

# Load isoform data
isoforms <- load_csv(isoforms_filename)
# Filter isoforms based on chromosome, position, and minexp
isoforms <- filter_isoforms(isoforms, chromosome, position, flanking, minexp)
# Format isoforms for merging
isoforms <- format_isoforms(isoforms)
# Merge isoforms with CSV file
isoforms <- combine_files(csv, isoforms)

# output the isoforms to a csv file
write.csv(isoforms, file.path(folder, paste0(trait_name, "_isoforms.csv")), row.names = FALSE)

# --- Load Data ---
log_msg("Loading genoprobs...", paste0(gp_filename, chromosome, ".rds"))
gp <- readRDS(paste0(gp_filename, chromosome, ".rds"))
#gp <- readRDS("Data/final_genoprobs_1176.rds")[,chromosome]
#test2 <- readRDS("Data/final_genoprobs_1176.rds")
log_msg("Loading map...")
map <- readRDS(map_filename)
log_msg("Loading kinship matrix...")
K <- readRDS(k_filename)
log_msg("Loading isoforms...")
pheno <- as.data.frame(isoforms)
colnames(pheno)[1] <- "mouse"
rownames(pheno) <- pheno$mouse
log_msg("All data loaded.")


# --- Get Gene Names ---
#get it from the isoforms header - the first x columns are not genes after the trait
trait_index <- which(colnames(pheno) == trait_name)
gene_names <- colnames(pheno)[(trait_index + 1):ncol(pheno)]
log_msg("Gene names extracted.")

# --- Make a kinship matrix and subset for the specified chromosome ---
K_chr <- K[[chromosome]]
log_msg("Kinship matrix and genoprobs subsetted.")

# --- Reclassify covars ---
pheno$Diet <- as.factor(pheno$Diet)
pheno$GenLit <- as.factor(pheno$GenLit)
pheno$Sex <- as.factor(pheno$Sex)
if ("Batch" %in% colnames(pheno)) {
  pheno$Batch <- as.factor(pheno$Batch)
}

# --- Create covariates ---
base_covar <- model.matrix(~Sex * Diet + GenLit, data = pheno)[, -1]
interactive_covar <- model.matrix(~Diet, data = pheno)[, -1, drop = FALSE]
log_msg("Covariates reclassified.")


#------ check LOD profile for trait without a transcript as an additive covar -----------
pheno[[paste0("rankz.", trait_name)]] <- rankz(pheno[[trait_name]])
qtl_no_additive <- scan1(genoprobs = gp, pheno = pheno[, paste0("rankz.", trait_name), drop = FALSE], kinship = K_chr, addcovar = base_covar, intcovar = interactive_covar)
log_msg("QTL scan for trait alone completed.")

qtl_no_additive[qtl_no_additive == 0] <- 1e-6 # Avoid division by zero
log_msg("LOD profile for trait without a transcript as an additive covar completed.")

peaks <- find_peaks(scan1_output = qtl_no_additive, map = map, threshold = threshold)
log_msg("Peaks found.")

#Check to make sure there is only one peak
if (nrow(peaks) > 1) {
  log_msg("There are multiple peaks, choosing the highest one within 2 Mb of the original peak.")
  # Find the highest peak within 2 Mb of the original peak
  peak_position <- peaks$pos[which.min(abs(peaks$pos - position))]
  peaks <- peaks[peaks$pos == peak_position,]
  log_msg("Peak position:", peak_position)
} else if (nrow(peaks) == 0) {
  stop_msg("There are no peaks.")
} else {
  log_msg("There is only one peak.")
}

#Get the peak position
peak_lod_position <- peaks$pos
# Get the highest LOD score
peak_lod <- peaks$lod

log_msg("Peak position:", peak_lod_position)
log_msg("Peak LOD score:", peak_lod)

# Set x_limit to the peak position +/- 10 Mb
x_limit <- c(max(peak_lod_position - 10, 0), min(peak_lod_position + 10, max(map[[chromosome]])))

# Set y_limit to 0 to the highest LOD score + 30%
y_limit <- c(0, peak_lod * 1.3)

log_msg("X_limit:", x_limit)
log_msg("Y_limit:", y_limit)

# Plot and save to PNG
png_filename <- paste0(folder, "/", trait_name, "_alone.png")
png(png_filename)
plot_scan1(x = qtl_no_additive, xlim = x_limit, ylim = y_limit, map = map, main = "QTL for trait alone")
dev.off()


# --- Function to calculate LOD difference and perform original QTL analysis ---
analyze_gene <- function(gene, genoprobs, pheno, K_chr, base_covar, interactive_covar, map, qtl_no_additive, initial_lod = NA) {
  log_msg("Analyzing gene:", gene)
  # Transform data
  pheno[[paste0("rankz.", gene)]] <- rankz(pheno[[gene]])

  # Perform scan with additive covariate
  addcovar_with_gene <- model.matrix(~Sex * Diet + GenLit + pheno[[paste0("rankz.", gene)]], data = pheno)[, -1]
  qtl_with_additive <- scan1(genoprobs = genoprobs, pheno = pheno[, trait_name, drop = FALSE], kinship = K_chr, addcovar = addcovar_with_gene, intcovar = interactive_covar, cores = 0)

  peaks <- find_peaks(scan1_output = qtl_with_additive, map = map, sort = "lod")

  # Check if there are any peaks
  if (nrow(peaks) > 1) {
    warning(paste("There are multiple peaks for gene", gene, "choosing the highest one within 2 Mb of the original peak."))
    # Find the highest peak within 2 Mb of the original peak
    peak_position <- peaks$pos[which(abs(peaks$pos - peak_lod_position) <= 2000000)]
    peaks <- peaks[peaks$pos == peak_position,]
    if (nrow(peaks) > 1) {
      warning(paste("There are still multiple peaks for gene", gene, "choosing the highest one."))
      peaks <- peaks[peaks$lod == max(peaks$lod),]
    }
  }

  if (nrow(peaks) == 0) {
    warning(paste("There are no peaks for gene", gene))
    peaks <- data.frame(
      lodindex = 1,
      lodcolumn = trait_name,
      chr = chromosome,
      pos = -1,
      lod = -1,
      gene = gene,
      percent_change = -1
    )
  } else {
    if (!is.na(initial_lod)) {
      # Add percent change
      peaks$percent_change <- (initial_lod - max(peaks$lod)) / initial_lod * 100
    }
    # Add gene name
    peaks$gene <- gene
  }


  # Plot and save to PNG if the LOD score has changed by more than the specified percent
  if (!is.null(initial_lod) && (nrow(peaks) == 0 || peaks$percent_change > percent)) {
    png_filename <- paste0(folder, "/", gene, ".png")
    png(png_filename)
    plot_scan1(x = qtl_with_additive, xlim = x_limit, ylim = y_limit, map = map, main = paste("QTL Scan for ", trait_name, " + ", gene))
    dev.off()
  }
  log_msg("Gene analysis completed:", gene)
  # Return peaks
  return(peaks)
}

# --- Set up Cluster ---
log_msg("Setting up cluster and exporting functions.")
cl <- makeCluster(numCores)
clusterExport(cl, c(
  "pheno", "gp", "K_chr", "base_covar", "interactive_covar", "rankz", "map", "find_peaks", "scan1",
  "analyze_gene", "trait_name", "qtl_no_additive", "plot_scan1", "folder", "x_limit", "y_limit",
  "chromosome", "threshold", "percent", "peak_lod", "peak_lod_position", "log_msg"
), envir = environment())

log_msg("Cluster set up and functions exported. ")

# --- Parallel Loop for Gene Analysis ---
analysis_results <- parLapply(cl, gene_names, function(gene) {
  source("logging.R")
  analyze_gene(gene, gp, pheno, K_chr, base_covar, interactive_covar, map, qtl_no_additive, initial_lod = peak_lod)
})

# --- Close Cluster ---
stopCluster(cl)
log_msg(" Cluster closed. ")

# --- Aggregate Results ---
analysis_df <- do.call(rbind, analysis_results)
log_msg(" Results aggregated. ")

# --- Remove first column if it is not lodindex ---
if (names(analysis_df)[1] != "lodindex") {
  analysis_df <- analysis_df[, -1]
}

# --- Add initial LOD score to analysis_results at the top ---
new_df <- data.frame(
  lodindex = 1,
  lodcolumn = trait_name,
  chr = chromosome,
  pos = peak_lod_position,
  lod = peak_lod,
  gene = "",
  percent_change = 0
)

# Sort analysis_df by percent_change in descending order
analysis_df <- analysis_df[order(-analysis_df$percent_change),]

# Add new_df to the top of sorted analysis_df
analysis_df <- rbind(new_df, analysis_df)

#Split gene into 2 columns, transcript and gene symbol by the first period
split_gene <- strsplit(as.character(analysis_df$gene), "\\.")

# Extracting the first part (transcript)
analysis_df$transcript <- sapply(split_gene, function(x) ifelse(length(x) > 0, x[1], NA))

# Extracting the second part (gene symbol) and handling rows without ' . '
analysis_df$gene_symbol <- sapply(split_gene, function(x) ifelse(length(x) > 1, x[2], NA))

# Remove gene column
analysis_df$gene <- NULL

#move the transcript and gene symbol to the front
analysis_df <- analysis_df[, c("transcript", "gene_symbol", "lodindex", "lodcolumn", "chr", "pos", "lod", "percent_change")]

# --- Remove max_diff and max_diff_loc from analysis_results  before saving ---
write.csv(analysis_df, file.path(folder, paste0(trait_name, "_analysis_results.csv")), row.names = FALSE)
log_msg("Results saved.")