logging_env <- new.env()
logging_env$log_file <- NULL
logging_env$mode <- "debug" # Example modes: "log", "debug", or "off"
# debug = print all messages and save to log file
# log = save all messages to log file
# silent = do not print or save any messages


# --- Set log file --- 
set_log_file <- function(folder = NULL) {
  if (!is.null(logging_env$log_file)) {
    close(logging_env$log_file)
    logging_env$log_file <- NULL
  }

  log_file_path <- if (!is.null(folder) && dir.exists(folder)) {
    paste0(folder, "/log.txt")
  } else {
    "log.txt"
  }

  logging_env$log_file_path <- log_file_path
  log("Log file path set to: ", log_file_path)

  logging_env$log_file <- file(log_file_path, open = "a")
}




# --- Log function ---
log <- function(..., collapse = " ", add_time = TRUE) {
  if (logging_env$mode %in% c("log", "debug")) {
    args <- list(...)
    message <- paste(sapply(args, as.character), collapse = collapse)



    if (logging_env$mode == "debug") {
      cat(message, "\n")
    }

    if (add_time) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      message <- paste(timestamp, message)
    }

    if (!is.null(logging_env$log_file)) {
      writeLines(message, con = logging_env$log_file)
      flush(logging_env$log_file)
    } else {
      cat("Log file not set.\n")
    }
  }
}
