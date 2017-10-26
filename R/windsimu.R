#' Simulate Datasets with User-Defined Patterns of Within-Individual Variation
#' 
#' Function to create simulated dataset given user-defined parameters and data
#' structure for growth curve model with optional random effect for L2 variance
#' 
#' @param sample_size List of named vectors specifying sample size at level 2 and level 1; order not important but need to be named (as \code{n2} and \code{n1}, respectively).
#' @param var_names Character vector containing names for intercept and time variable (order not important, but elements must be
#' named \code{x0} (corresponding to intercept) and \code{x1} (time); e.g. \code{c("x0" = "cons", "x1" = "age")}).
#' @param x1_range Vector with min and max values for time variable (in that order; e.g. \code{c(12, 23)} to range from 12 to 23).
#' @param design_matrices List of logical vectors indicating whether each variable to be included in design matrices or not; order of elements in each logical vector not important, but must be
#' named \code{x0} (intercept) and \code{x1} (time). \code{TRUE} indicates variable to be included, \code{FALSE} otherwise. List elements themselves must be named \code{x_A} (corresponding to X matrix (fixed part) of model for mean of y),
#' \code{z_A} (Z matrix (random part) of model for mean of y), \code{x_B} (X matrix (fixed part) of model for level 1 variance function),
#' \code{z_B} (Z matrix (random part) of model for level 1 variance function). E.g. \code{list(x_A = c("x0" = TRUE, "x1" = TRUE), z_A = c("x0" = TRUE, "x1" = TRUE), x_B = c("x0" = TRUE, "x1" = TRUE), z_B = c("x0" = TRUE, "x1" = FALSE))}.
#' @param beta Vector containing mean for each beta (i.e. coefficients of fixed effects for model for mean of y), in order of x0, x1.
#' @param alpha Vector containing mean for each alpha (i.e. coefficients of fixed effects for level 1 variance function), in order of x0, x1.
#' @param sigma2u Level 2 covariance matrix. E.g. the vector in \code{matrix(c(93, 12, 2.5, 12, 3.8, 0.8, 2.5, 0.8, 0.4), nrow = 3, ncol = 3, byrow = TRUE)}
#' would correspond to the following: sigma2u00 (\code{93}), sigma2u01 (\code{12}), sigma2u02 (\code{2.5}), sigma2u01 (\code{12}), sigma2u11 (\code{3.8}), sigma2u12 (\code{0.8}), sigma2u02 (\code{2.5}),
#' sigma2u12 (\code{0.8}), sigma2u22 (\code{0.4}).
#' @param sigma2e Numerical vector (one element); only specify if 1-level model; defaults to \code{NULL}.
#' @param log_sigma2e Logical vector indicating whether log link used for level 1 variance function (\code{TRUE}) or not (\code{FALSE});
#' defaults to \code{TRUE}.
#' @param linear_sigma2e_attempts Numerical vector (one element) corresponding to maximum number of attempts made to generate non-negative sigma2e (only applicable when \code{log_sigma2e = FALSE});
#' defaults to \code{10}.
#' @param random_seed Numerical vector indicating value of random seed to be set at outset.
#' @return List containing simulated dataset (as \code{simu_dataframe}), plus other relevant generated objects and user-defined parameters.
#' @examples
#' \dontrun{
#' library(windsimu)
#'
#' ## Level 2 random effect in level 1 variance function
#' ## with this term also allowed to covary with random
#' ## intercept & random slope in level 2 covariance matrix:
#' simu_u2j_covary <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = c("x0" = "cons", "x1" = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = c("x0" = TRUE, "x1" = TRUE),
#'     z_A = c("x0" = TRUE, "x1" = TRUE),
#'     x_B = c("x0" = TRUE, "x1" = TRUE),
#'     z_B = c("x0" = TRUE, "x1" = FALSE)
#'   ),
#'   beta = c(150, 6.5),
#'   alpha = c(-.95, 0.48),
#'   sigma2u =
#'     matrix(
#'       c(93, 12, 2.5,
#'         12, 3.8, 0.8,
#'         2.5, 0.8, 0.4),
#'       nrow = 3,
#'       ncol = 3,
#'       byrow = TRUE
#'     ),
#'   log_sigma2e = TRUE,
#'   linear_sigma2e_attempts = 10,
#'   random_seed = 1
#' )
#' 
#' ## inspect structure of returned object:
#' str(simu_u2j_covary)
#' 
#' ## export simulated dataset (to current working directory) as .dta:
#' library(haven)
#' write_dta(simu_u2j_covary$simu_dataframe, "simu_u2j_covary.dta", version = 14)
#' 
#' 
#' ## Level 2 random effect in level 1 variance function
#' ## but covariance of this term set to 0 in level 2 covariance matrix:
#' simu_u2j_no_covary <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = c("x0" = "cons", "x1" = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = c("x0" = TRUE, "x1" = TRUE),
#'     z_A = c("x0" = TRUE, "x1" = TRUE),
#'     x_B = c("x0" = TRUE, "x1" = TRUE),
#'     z_B = c("x0" = TRUE, "x1" = FALSE)
#'   ),
#'   beta = c(150, 6.5),
#'   alpha = c(-.95, 0.48),
#'   sigma2u =
#'     matrix(
#'       c(93, 12, 0,
#'         12, 3.8, 0,
#'         0, 0, 0.4),
#'       nrow = 3,
#'       ncol = 3,
#'       byrow = TRUE
#'     ),
#'   sigma2e = NULL,
#'   log_sigma2e = TRUE,
#'   linear_sigma2e_attempts = 10,
#'   random_seed = 1
#' )
#' 
#' 
#' ## No additional level 2 random effect:
#' simu_no_u2j <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = c("x0" = "cons", "x1" = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = c("x0" = TRUE, "x1" = TRUE),
#'     z_A = c("x0" = TRUE, "x1" = TRUE),
#'     x_B = c("x0" = TRUE, "x1" = TRUE),
#'     z_B = c("x0" = FALSE, "x1" = FALSE)
#'   ),
#'     beta = c(150, 6.5),
#'     alpha = c(-.95, 0.48),
#'     sigma2u =
#'       matrix(
#'         c(93, 12,
#'           12, 3.8),
#'         nrow = 2,
#'         ncol = 2,
#'         byrow = TRUE
#'       ),
#'     log_sigma2e = TRUE,
#'     linear_sigma2e_attempts = 10,
#'     random_seed = 1
#'   )
#' 
#' 
#' ## No log link for level 1 variance function.
#' ## Note age uncentred (level 1 variance function more likely
#' ## to stay positive); also note linear_sigma2e_attempts:
#' simu_u2j_covary_linearL1varfun_Revised <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = c("x0" = "cons", "x1" = "age"),
#'   x1_range = c(11.25, 13.25),
#'   design_matrices = list(
#'     x_A = c("x0" = TRUE, "x1" = TRUE),
#'     z_A = c("x0" = TRUE, "x1" = TRUE),
#'     x_B = c("x0" = TRUE, "x1" = TRUE),
#'     z_B = c("x0" = TRUE, "x1" = FALSE)
#'   ),
#'   beta = c(150, 6.5),
#'   alpha = c(-.95, 0.48),
#'   sigma2u =
#'     matrix(
#'       c(93, 12, 2.5,
#'         12, 3.8, 0.8,
#'         2.5, 0.8, 0.4),
#'       nrow = 3,
#'       ncol = 3,
#'       byrow = TRUE
#'     ),
#'   sigma2e = NULL,
#'   log_sigma2e = FALSE,
#'   linear_sigma2e_attempts = 100,
#'   random_seed = 1
#' )
#' }
#'
#' @export
windsimu <- function(sample_size,
                     var_names,
                     x1_range,
                     design_matrices,
                     beta,
                     alpha,
                     sigma2u,
                     sigma2e = NULL,
                     log_sigma2e = TRUE,
                     linear_sigma2e_attempts = 10,
                     random_seed) {

  set.seed(random_seed)
  
  # ensure order of elements consistent ---------------------------
  # (not necessary to order sample_size as only refer to elements of it by name)
  var_names <- var_names[c("x0", "x1")]
  design_matrices <- lapply(design_matrices, function(x) x[c("x0", "x1")])
  
  # create vectors ---------------------------
  total_L1 <- sample_size$n2 * sample_size$n1
  L2_ID <- rep(1:sample_size$n2, each = sample_size$n1)
  time_var <- rep(seq(from = x1_range[1], to = x1_range[2], length.out = sample_size$n1), times = sample_size$n2)
  how_many_variables_x_matrix_A <- sum(design_matrices$x_A)
  how_many_variables_z_matrix_A <- sum(design_matrices$z_A)
  how_many_variables_x_matrix_B <- sum(design_matrices$x_B)
  how_many_variables_z_matrix_B <- sum(design_matrices$z_B)
  predictors <- matrix(nrow = total_L1, ncol = 2)
  predictors[, 1] <- 1 #constant
  predictors[, 2] <- time_var
  
  # create design matrices, and multiply out (latter for fixed part only) ---------------------------
  if (any(design_matrices$x_A)){
    x_matrix_A <- as.matrix(predictors[, design_matrices$x_A])
    colnames(x_matrix_A) <- var_names[design_matrices$x_A]
    fixpart_A <- x_matrix_A %*% beta
  } else {
    x_matrix_A <- NULL
  }
  
  if (any(design_matrices$z_A)){
    z_matrix_A <- as.matrix(predictors[, design_matrices$z_A])
    colnames(z_matrix_A) <- var_names[design_matrices$z_A]
  } else {
    z_matrix_A <- NULL
  }
  
  if (any(design_matrices$x_B)){
    x_matrix_B <- as.matrix(predictors[, design_matrices$x_B])
    colnames(x_matrix_B) <- var_names[design_matrices$x_B]
    fixpart_B <- x_matrix_B %*% alpha
  } else {
    x_matrix_B <- NULL
  }
  
  if (any(design_matrices$z_B)){
    z_matrix_B <- as.matrix(predictors[, design_matrices$z_B])
    colnames(z_matrix_B) <- var_names[design_matrices$z_B]
  } else {
    z_matrix_B <- NULL
  }

  # Create sigma_e ---------------------------
  # This is all wrapped up in a function (called later),
  # as may need to run whole thing repeatedly to get non-negative
  # sigma2_e if log_sigma2e = FALSE

  create_sigma2_e <- function(n2,
                                how_many_variables_z_matrix_A,
                                how_many_variables_z_matrix_B,
                                sigma2u,
                                sigma2e,
                                z_matrix_B = NULL,
                                fixpart_B) {
    u <- MASS::mvrnorm(n = sample_size$n2,
                 mu = rep(0, how_many_variables_z_matrix_A + how_many_variables_z_matrix_B),
                 Sigma = sigma2u
    )
    
    u_names <- NULL
    for (i in seq(ncol(u))) {
      u_names[i] <- paste0("u", i-1, "j")
    }
    colnames(u) <- u_names # for readability
    
    if(is.null(z_matrix_B)){
      randpart_B <- NULL
      sigma2_e <- fixpart_B
    } else {
      randpart_B <- rowSums(z_matrix_B * u[L2_ID, "u2j"])
      sigma2_e <- fixpart_B + randpart_B
    }
    
    list(u = u,
         randpart_B = randpart_B,
         sigma2_e = sigma2_e
    )
  }

  # if 1-level model:
  if (how_many_variables_z_matrix_A + how_many_variables_z_matrix_B == 0) {
    sigma_e <- sqrt(sigma2e)
  } else {
  # if 2-level model:
    if (log_sigma2e == FALSE) {
      for (i in linear_sigma2e_attempts) {
        function_call <- create_sigma2_e(
          n2 = sample_size$n2,
          how_many_variables_z_matrix_A = how_many_variables_z_matrix_A,
          how_many_variables_z_matrix_B = how_many_variables_z_matrix_B,
          sigma2u = sigma2u,
          z_matrix_B = z_matrix_B,
          fixpart_B = fixpart_B
        )
        if (all(function_call$sigma2_e >= 0)) {
          break
        }
        if (i == linear_sigma2e_attempts)
          stop(
            "No viable simulations of sigma2e achieved within max number of iterations specified."
          )
      }
      sigma_e <- sqrt(function_call$sigma2_e)
    } else {
      function_call <- create_sigma2_e(
        n2 = sample_size$n2,
        how_many_variables_z_matrix_A = how_many_variables_z_matrix_A,
        how_many_variables_z_matrix_B = how_many_variables_z_matrix_B,
        sigma2u = sigma2u,
        z_matrix_B = z_matrix_B,
        fixpart_B = fixpart_B
      )
      sigma_e <- sqrt(exp(function_call$sigma2_e))
    }
  }
  
  # Create level 1 residuals ---------------------------
  
  e <- rnorm(n = total_L1, mean = 0, sd = sigma_e)

  # Create randpart_A; generate y ---------------------------
  
  # if 1-level model:
  if (how_many_variables_z_matrix_A + how_many_variables_z_matrix_B == 0) {
    
    y <- fixpart_A + e
  
    } else {
  
  # if 2-level model:
  randpart_A <- rowSums(z_matrix_A * function_call$u[L2_ID, c("u0j", "u1j")])
  
  y <- fixpart_A + randpart_A + e
  
  }
  
  # Create data frame ---------------------------
  L1_ID <- 1:total_L1
  simu_dataframe <- data.frame(L2_ID, L1_ID, y, predictors)
  colnames(simu_dataframe)[4:5] <- var_names
  
  # Return list containing relevant outputs ---------------------------
  list(
    simu_dataframe = simu_dataframe,
    x_matrix_A = x_matrix_A,
    z_matrix_A = z_matrix_A,
    x_matrix_B = x_matrix_B,
    z_matrix_B = z_matrix_B,
    fixpart_A = fixpart_A,
    fixpart_B = fixpart_B,
    randpart_A = randpart_A,
    randpart_B = function_call$randpart_B,
    u = function_call$u,
    sigma2_e = function_call$sigma2_e,
    sigma_e = sigma_e,
    log_sigma2e = log_sigma2e,
    alpha = alpha,
    beta = beta,
    sample_size = sample_size,
    total_L1 = total_L1
  )
}
