#' Example data, simulated from a Weibull distribution.
#'
#' Data is simulated for a concurrent trial with three columns named "tte" 
#' (time-to-event), "event" (event indicator), and "X_trt" (treatment indicator).
#' It was simulated by drawing samples from a Weibull with kappa = 1.5 (shape) 
#' and nu = 0.4 (scale)
#'
#' @docType data
#'
#' @usage data(weibull_cc)
#'
#' @keywords datasets
#'
#' @examples
#' data(weibull_cc)
#' survival_model <- survival::survfit(survival::Surv(tte, event) ~ X_trt, data = weibull_cc)
#' line_colors <- c("blue", "red")  # Adjust colors as needed
#' line_types <- 1:length(unique(weibull_cc$X_trt))
#' plot(survival_model, col = line_colors, lty = line_types, 
#'      xlab = "Time (tte)", ylab = "Survival Probability", 
#'      main = "Kaplan-Meier Survival Curves by Treatment")
"weibull_cc"