#' Example data, simulated from a piecewise exponential model.
#'
#' Data is simulated for a concurrent trial with three columns named "tte" 
#' (time-to-event), "event" (event indicator), and "X_trt" (treatment indicator).
#' It was simulated using the following parameters: 
#'
#' @docType data
#'
#' @usage data(piecewise_exp_cc)
#'
#' @keywords datasets
#'
#' @examples
#' data(piecewise_exp_cc)
#' survival_model <- survival::survfit(survival::Surv(tte, event) ~ X_trt, data = piecewise_exp_cc)
#' line_colors <- c("blue", "red")  # Adjust colors as needed
#' line_types <- 1:length(unique(piecewise_exp_cc$X_trt))
#' plot(survival_model, col = line_colors, lty = line_types, 
#'      xlab = "Time (tte)", ylab = "Survival Probability", 
#'      main = "Kaplan-Meier Survival Curves by Treatment")
"piecewise_exp_cc"