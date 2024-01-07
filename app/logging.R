# Environment to hold logging settings and state
logging_env <- new.env()
logging_env$log_file_path <- NULL
logging_env$mode <- "debug" # Modes: "log", "debug", "off"
# debug = print & save messages
# log = save messages only
# off = neither print nor save

# --- Set log file path and initialize log file ---
set_log_file <- function(folder = NULL) {
  # Close existing log file if open
  if (!is.null(logging_env$log_file)) {
    close(logging_env$log_file)
    logging_env$log_file <- NULL
  }

  # Determine log file path
  log_file_path <- if (!is.null(folder) && dir.exists(folder)) {
    paste0(folder, "/log.txt")
  } else {
    "log.txt"
  }

  logging_env$log_file_path <- log_file_path

  # Open new log file for appending
  
}

# --- Function to record log messages ---
log_msg <- function(..., collapse = " ", add_time = TRUE) {
  # Check if logging mode requires message recording
  if (logging_env$mode %in% c("log", "debug")) {
    args <- list(...)
    message <- paste(sapply(args, as.character), collapse = collapse)

    # Print message in debug mode
    if (logging_env$mode == "debug") {
      cat(message, "\n")
    }

    # Add timestamp if required to file message
    if (add_time) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      message <- paste(timestamp, message)
    }

    # Write message to log file if set
    if (!is.null(logging_env$log_file_path)) {
      log_file <- file(logging_env$log_file_path, open = "a")
      writeLines(message, con = log_file)
      flush(log_file)
      close(log_file)
    } else {
      warning("Log file not set.\n")
    }
  }
}

# --- function to stop and print error message ---
stop_msg <- function(..., collapse = " ", add_time = TRUE) {
  args <- list(...)
  message <- paste(sapply(args, as.character), collapse = collapse)

  # Add timestamp if required to file message
  if (add_time) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message <- paste(timestamp, message)
  }

  # Write message to log file if set
  if (!is.null(logging_env$log_file_path)) {
    log_file <- file(logging_env$log_file_path, open = "a")
    writeLines(message, con = log_file)
    flush(log_file)
    close(log_file)
  } else {
    warning("Log file not set.\n")
  }

  stop(message)
}
