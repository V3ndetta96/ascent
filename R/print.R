# ==============================================================================
# PRINT METHODS FOR ALL S3 CLASSES
# ==============================================================================

#' @export
print.ascfcd <- function(x, ...) {
  cat("ASC-CFD Paired Analysis\n")
  cat(sprintf("  Contrasts: %d | Dimensions: %d (%.1f%% var | Quality: %.1f%%)\n",
              length(x$rDelta_C), x$k_retained, x$var_retained * 100, x$quality * 100))
  cat(sprintf("  Null models: %s\n", if(!is.null(x$null_models)) "evaluated" else "not evaluated"))
  cat("Use summary() for full results.\n")
  invisible(x)
}

#' @export
print.ascfcd_pw <- function(x, ...) {
  cat("ASC-CFD Pairwise Spatial Analysis\n")
  cat(sprintf("  Communities: %d | Contrasts: %d | Dimensions: %d (%.1f%% var | Quality: %.1f%%)\n",
              nrow(x$cwm_global), nrow(x$pairwise_results), x$k_retained,
              x$var_retained * 100, x$quality * 100))
  cat(sprintf("  Null models: %s\n", if(!is.null(x$null_models)) "evaluated" else "not evaluated"))
  cat("Use summary() for full results.\n")
  invisible(x)
}

#' @export
print.ascfcd_base <- function(x, ...) {
  cat("ASC-CFD Baseline Functional Topology\n")
  cat(sprintf("  Entities: %d | Dimensions: %d (%.1f%% var | Quality: %.1f%%)\n",
              nrow(x$entities_results), x$k_retained, x$var_retained * 100, x$quality * 100))
  cat("Use summary() for full results.\n")
  invisible(x)
}

#' @export
print.ascfcd_entities <- function(x, ...) {
  n_ent <- length(unique(x$entity_classification$Entity_ID))
  n_spp <- nrow(x$entity_classification)
  cat("ASC-CFD Functional Entities Classification\n")
  cat(sprintf("  Species: %d | Entities: %d | Metric: %s | Method: %s\n",
              n_spp, n_ent, x$metric, x$method))
  cat("Use summary() for full results.\n")
  invisible(x)
}

#' @export
print.ascfcd_transitions <- function(x, ...) {
  cat("ASC-CFD Functional Leverage Analysis\n")
  cat(sprintf("  Contrasts: %d\n", length(x)))
  cat("Use summary() for detailed results.\n")
  invisible(x)
}
