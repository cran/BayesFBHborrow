#' Example data, simulated from a Weibull distribution
#'
#' Data is simulated for a historical trial with two columns named "tte" 
#' (time-to-event) and "event" (event indicator).
#' It was simulated using the following parameters: 
#'
#' @docType data
#'
#' @usage data(weibull_hist)
#'
#' @keywords datasets
#'
#' @examples
#' data(weibull_cc)
#' data(weibull_hist)
#' weibull_hist$X_trt <- 0
#' survival_model <- survival::survfit(survival::Surv(tte, event) ~ X_trt, 
#'                                     data = rbind(weibull_cc, 
#'                                     weibull_hist))
#' line_colors <- c("blue", "red", "green")  # Adjust colors as needed
#' line_types <- 1:length(unique(weibull_cc$X_trt))
#' plot(survival_model, col = line_colors, lty = line_types, 
#'      xlab = "Time (tte)", ylab = "Survival Probability", 
#'      main = "Kaplan-Meier Survival Curves by Treatment")
"weibull_hist"