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

mode <- "expression" # Set to "physiological" or "expression" to change logging mode
# if expression use the csv file to get the trait data
csv_filename <- "20230627_DOplasmalipids.csv"
trait_name <- "TC" #Apobec3 
chromosome <- "5"
position <- 125000000 # Position in base pairs
flanking <- 2000000 # Flanking region in base pairs
minexp <- 100 # Minimum non-zero expression value for isoforms
percent <- 0.3 #Percent LOD change to plot graph
threshold <- 20 # Threshold for LOD score

# Define a function to rank data
rankz = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

numCores <- detectCores() - 1 # Keeping one core free for other processes


# Non-interactive parameters
folder <- "Results"
add_timestamp <- TRUE # Set to TRUE to add a timestamp, or FALSE to exclude it
isoforms_filename <- "vsd_DO_isoforms_new_app.csv"
generate_plots <- TRUE #TODO: Remove this

# Create folder and set log file
folder <- create_folder(folder, add_timestamp)
set_log_file(folder)
log_msg("Folder for results created:", folder)
if (mode == "expression") {
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
} else {
  # Load isoform data
  isoforms <- load_csv(isoforms_filename)
  # Filter isoforms for those containing the trait
  isoforms <- filter_isoforms_for_trait(isoforms, trait_name)
}

# --- Load Data ---
log_msg("Loading genoprobs...")
genoprobs <- readRDS("final_genoprobs_1176.rds")
log_msg("Loading map...")
map <- readRDS("gigamuga_map_test_v2.RDS")
log_msg("Loading kinship matrix...")
K <- readRDS("K_1176_DietDO.rds")
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
gp <- genoprobs[, chromosome]
log_msg("Kinship matrix and genoprobs subsetted.")

# --- Reclassify covars ---
pheno$Diet <- as.factor(pheno$Diet)
pheno$GenLit <- as.factor(pheno$GenLit)
pheno$Sex <- as.factor(pheno$Sex)

base_covar <- model.matrix(~Sex * Diet + GenLit, data = pheno)[, -1]
interactive_covar <- t(t(model.matrix(~Diet, data = pheno)[, -1]))
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
  log_msg("There is more than one peak.")
  #Exit
  q()
} else {
  log_msg("There is only one peak.")
}

#Get the peak position
peak_position <- peaks$pos
# Get the highest LOD score
peak_lod <- peaks$lod

log_msg("Peak position:", peak_position)
log_msg("Peak LOD score:", peak_lod)

#Set x_limit to the peak position +/- 10 Mb
x_limit <- c(peak_position - 10, peak_position + 10)
# Set y_limit to 0 to the highest LOD score + 30%
y_limit <- c(0, peak_lod * 1.3)

log_msg("X_limit:", x_limit)
log_msg("Y_limit:", y_limit)

# Plot and save to PNG
if (generate_plots) {
  png_filename <- paste0(folder, "/", trait_name, "_alone.png")
  png(png_filename)
  plot_scan1(x = qtl_no_additive, xlim = x_limit, ylim = y_limit, map = map, main = "QTL for trait alone")
  dev.off()
}

# --- Function to calculate LOD difference and perform original QTL analysis ---
analyze_gene <- function(gene, genoprobs, pheno, K_chr, base_covar, interactive_covar, map, qtl_no_additive, initial_lod = NA) {
  # Transform data
  pheno[[paste0("rankz.", gene)]] <- rankz(pheno[[gene]])

  # Perform scan with additive covariate
  addcovar_with_gene <- model.matrix(~Sex * Diet + GenLit + pheno[[paste0("rankz.", gene)]], data = pheno)[, -1]
  qtl_with_additive <- scan1(genoprobs = genoprobs, pheno = pheno[, trait_name, drop = FALSE], kinship = K_chr, addcovar = addcovar_with_gene, intcovar = interactive_covar, cores = 0)

  # Find peaks for the additive model
  peaks <- find_peaks(scan1_output = qtl_with_additive, map = map, threshold = threshold)

  # Plot and save to PNG if the LOD score has changed by more than the specified percent
  if (!is.null(initial_lod) && abs(max(peaks$lod) - initial_lod) / initial_lod > percent) {
    png_filename <- paste0(folder, "/", gene, ".png")
    png(png_filename)
    plot_scan1(x = qtl_with_additive, xlim = x_limit, ylim = y_limit, map = map, main = paste("QTL Scan for ", trait_name, " + ", gene))
    dev.off()
  }

  # Add gene name
  peaks$gene <- gene

  #if there is a percent change in the LOD score, return the percent change
  if (!is.null(initial_lod)) {
    peaks$percent_change <- abs(max(peaks$lod) - initial_lod) / initial_lod
  }
  return(peaks)
}
log_msg("Setting up cluster and exporting functions.")
# --- Set up Cluster ---
cl <- makeCluster(numCores)
clusterExport(cl, c(
  "pheno", "gp", "K_chr", "base_covar", "interactive_covar", "rankz", "map", "find_peaks", "scan1",
  "analyze_gene", "trait_name", "qtl_no_additive", "generate_plots", "plot_scan1", "folder", "x_limit", "y_limit", "chromosome",
  "threshold", "percent"))
log_msg("Cluster set up and functions exported.")

# --- Parallel Loop for Gene Analysis ---
analysis_results <- parLapply(cl, gene_names, function(gene) {
  analyze_gene(gene, gp, pheno, K_chr, base_covar, interactive_covar, map, qtl_no_additive, initial_lod = peak_lod)
})

# --- Close Cluster ---
stopCluster(cl)
log_msg("Cluster closed.")

# --- Aggregate Results ---
analysis_df <- do.call(rbind, analysis_results)
log_msg("Results aggregated.")

# --- Add initial LOD score to analysis_results at the top ---
analysis_df <- rbind(data.frame(gene = trait_name, pos = peak_position, lod = peak_lod, percent_change = 0), analysis_df)
log_msg("Initial LOD score added.")

# --- Remove max_diff and max_diff_loc from analysis_results before saving ---
write.csv(analysis_df, file.path(folder, paste0(trait_name, "_analysis_results.csv")), row.names = FALSE)
log_msg("Results saved.")

log_msg("LOD difference results saved.")
