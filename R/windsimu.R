#' Simulate Datasets with User-Defined Patterns of Within-Individual Variation
#' 
#' Function to create simulated dataset given user-defined parameters and data
#' structure for growth curve model with optional random effect for L2 variance
#' 
#' @param how_many_L2_per_L3 Sample size at level 2 (there isn't an L3, but future compatible in case there ever is!)
#' @param how_many_L1_per_L2 Sample size at level 1, per level 2 unit
#' @param time_var_range Vector with min and max for time variable (e.g. \code{c(12, 23)} to range from 12 to 23).
#' @param var_names Character vector containing names for intercept and time variable (e.g. \code{c("cons", "age")}).
#' @param x_matrix_A_cons_time Logical vector indicating whether intercept (cons) and time variables to be included in X design matrix A
#' (i.e. fixed part of model for mean of y); e.g. \code{c(TRUE, FALSE)} indicates intercept, but not time, to be included.
#' @param z_matrix_A_cons_time As for \code{x_matrix_A_cons_time}, but Z matrix.
#' @param x_matrix_B_cons_time As for \code{x_matrix_A_cons_time}, but X matrix for level 1 variance function.
#' @param z_matrix_B_cons_time As for \code{x_matrix_A_cons_time}, but Z matrix for level 1 variance function.
#' @param beta Vector containing mean for each beta (i.e. i.e. coefficients of fixed effects for model for mean of y).
#' @param alpha Vector containing mean for each alpha (i.e. i.e. coefficients of fixed effects for level 1 variance function).
#' @param sigma2u Level 2 covariance matrix; e.g. \code{matrix(c(93, 12, 2.5, 12, 3.8, 0.8, 2.5, 0.8, 0.4), nrow = 3, ncol = 3, byrow = TRUE)}.
#' @param sigma2e Vector (one element); only specify if 1-level model; defaults to \code{NULL}.
#' @param log_sigma2e Logical vector indicating whether log link used for level 1 variance function (\code{TRUE}) or not (\code{FALSE}); defaults to \code{TRUE}.
#' @param linear_sigma2e_attempts Maximum number of attempts made to generate non-negative sigma2e (only applicable when \code{log_sigma2e = FALSE}); defaults to \code{10}.
#' @param random_seed Random seed adopted.
#' 
#' @examples
#' \dontrun{
#' library(windsimu)
#' 
#' simu_u2j_covary <- within_ind_var_simu_A(
#' how_many_L2_per_L3 = 1000,
#' how_many_L1_per_L2 = 9,
#' x1_range = c(-1, 1),
#' var_names = c("cons", "age"),
#' x_matrix_A_x0_x1 = c(TRUE, TRUE),
#' z_matrix_A_x0_x1 = c(TRUE, TRUE),
#' x_matrix_B_x0_x1 = c(TRUE, TRUE),
#' z_matrix_B_x0_x1 = c(TRUE, FALSE),
#' beta = c(150, 6.5),
#' alpha = c(-.95, 0.48),
#' sigma2u =
#'   matrix(
#'     c(93, 12, 2.5,
#'       12, 3.8, 0.8,
#'       2.5, 0.8, 0.4),
#'     nrow = 3,
#'     ncol = 3,
#'     byrow = TRUE
#'   ),
#' log_sigma2eij = TRUE,
#' linear_sigma2eij_attempts = 10,
#' random_seed = 1
#' )
#' 
#' simu_u2j_no_covary <- windsimu(
#' how_many_L2_per_L3 = 1000,
#' how_many_L1_per_L2 = 9,
#' time_var_range = c(-1, 1),
#' var_names = c("cons", "age"),
#' x_matrix_A_cons_time = c(TRUE, TRUE),
#' z_matrix_A_cons_time = c(TRUE, TRUE),
#' x_matrix_B_cons_time = c(TRUE, TRUE),
#' z_matrix_B_cons_time = c(TRUE, FALSE),
#' beta = c(150, 6.5),
#' alpha = c(-.95, 0.48),
#' sigma2u =
#'   matrix(
#'     c(93, 12, 0,
#'       12, 3.8, 0,
#'       0, 0, 0.4),
#'     nrow = 3,
#'     ncol = 3,
#'     byrow = TRUE
#'   ),
#' sigma2e = NULL,
#' log_sigma2e = TRUE,
#' linear_sigma2e_attempts = 10,
#' random_seed = 1
#' )
#' 
#' simu_no_u2j <- within_ind_var_simu_A(
#'   how_many_L2_per_L3 = 1000,
#'   how_many_L1_per_L2 = 9,
#'   x1_range = c(-1, 1),
#'   var_names = c("cons", "age"),
#'   x_matrix_A_x0_x1 = c(TRUE, TRUE),
#'   z_matrix_A_x0_x1 = c(TRUE, TRUE),
#'   x_matrix_B_x0_x1 = c(TRUE, TRUE),
#'   z_matrix_B_x0_x1 = c(FALSE, FALSE),
#'   beta = c(150, 6.5),
#'   alpha = c(-.95, 0.48),
#'   sigma2u =
#'     matrix(
#'       c(93, 12,
#'         12, 3.8),
#'       nrow = 2,
#'       ncol = 2,
#'       byrow = TRUE
#'     ),
#'   log_sigma2eij = TRUE,
#'   linear_sigma2eij_attempts = 10,
#'   random_seed = 1
#' )
#' 
#' ## Note age uncentred
#' simu_u2j_covary_linearL1varfun_Revised <- windsimu(
#'   how_many_L2_per_L3 = 1000,
#'   how_many_L1_per_L2 = 9,
#'   time_var_range = c(11.25, 13.25),
#'   var_names = c("cons", "age"),
#'   x_matrix_A_cons_time = c(TRUE, TRUE),
#'   z_matrix_A_cons_time = c(TRUE, TRUE),
#'   x_matrix_B_cons_time = c(TRUE, TRUE),
#'   z_matrix_B_cons_time = c(TRUE, FALSE),
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

windsimu <- function(how_many_L2_per_L3,
                     how_many_L1_per_L2,
                     time_var_range,
                     var_names = c("cons", "age"),
                     x_matrix_A_cons_time = c(TRUE, TRUE),
                     z_matrix_A_cons_time = c(TRUE, TRUE),
                     x_matrix_B_cons_time = c(TRUE, TRUE),
                     z_matrix_B_cons_time = c(TRUE, FALSE),
                     beta,
                     alpha,
                     sigma2u,
                     sigma2e = NULL,
                     log_sigma2e = TRUE,
                     linear_sigma2e_attempts = 10,
                     random_seed) {
  
  library(MASS)
  set.seed(random_seed)
  
  ## ----create scalars, vectors
  total_L1 <- how_many_L2_per_L3 * how_many_L1_per_L2
  L2_ID <- rep(1:how_many_L2_per_L3, each = how_many_L1_per_L2)
  time_var <- rep(seq(from = time_var_range[1], to = time_var_range[2], length.out = how_many_L1_per_L2), times = how_many_L2_per_L3)
  how_many_variables_x_matrix_A <- sum(x_matrix_A_cons_time)
  how_many_variables_z_matrix_A <- sum(z_matrix_A_cons_time)
  how_many_variables_x_matrix_B <- sum(x_matrix_B_cons_time)
  how_many_variables_z_matrix_B <- sum(z_matrix_B_cons_time)
  predictors <- matrix(nrow = total_L1, ncol = 2)
  predictors[, 1] <- 1 #constant
  predictors[, 2] <- time_var
  
  ## ----create design matrices, and multiply out (latter for fixed part only):
  if (any(x_matrix_A_cons_time)){
    x_matrix_A <- as.matrix(predictors[, x_matrix_A_cons_time])
    colnames(x_matrix_A) <- var_names[x_matrix_A_cons_time]
    fixpart_A <- x_matrix_A %*% beta
  } else {
    x_matrix_A <- NULL
  }
  
  if (any(z_matrix_A_cons_time)){
    z_matrix_A <- as.matrix(predictors[, z_matrix_A_cons_time])
    colnames(z_matrix_A) <- var_names[z_matrix_A_cons_time]
  } else {
    z_matrix_A <- NULL
  }
  
  if (any(x_matrix_B_cons_time)){
    x_matrix_B <- as.matrix(predictors[, x_matrix_B_cons_time])
    colnames(x_matrix_B) <- var_names[x_matrix_B_cons_time]
    fixpart_B <- x_matrix_B %*% alpha
  } else {
    x_matrix_B <- NULL
  }
  
  if (any(z_matrix_B_cons_time)){
    z_matrix_B <- as.matrix(predictors[, z_matrix_B_cons_time])
    colnames(z_matrix_B) <- var_names[z_matrix_B_cons_time]
  } else {
    z_matrix_B <- NULL
  }

  ## ----Create sigma_e; this is all wrapped up in a function (called later),
  ## as may need to run whole thing repeatedly to get non-negative
  ## sigma2_e if log_sigma2e = FALSE

  create_sigma2_e <- function(how_many_L2_per_L3,
                                how_many_variables_z_matrix_A,
                                how_many_variables_z_matrix_B,
                                sigma2u,
                                sigma2e,
                                z_matrix_B = NULL,
                                fixpart_B) {
    u <- mvrnorm(n = how_many_L2_per_L3,
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
          how_many_L2_per_L3 = how_many_L2_per_L3,
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
        how_many_L2_per_L3 = how_many_L2_per_L3,
        how_many_variables_z_matrix_A = how_many_variables_z_matrix_A,
        how_many_variables_z_matrix_B = how_many_variables_z_matrix_B,
        sigma2u = sigma2u,
        z_matrix_B = z_matrix_B,
        fixpart_B = fixpart_B
      )
      sigma_e <- sqrt(exp(function_call$sigma2_e))
    }
  }
  
  ## ----Create level 1 residuals
  
  e <- rnorm(n = total_L1, mean = 0, sd = sigma_e)

  ## ----Create randpart_A; generate y
  
  # if 1-level model:
  if (how_many_variables_z_matrix_A + how_many_variables_z_matrix_B == 0) {
    
    y <- fixpart_A + e
  
    } else {
  
  # if 2-level model:
  randpart_A <- rowSums(z_matrix_A * function_call$u[L2_ID, c("u0j", "u1j")])
  
  y <- fixpart_A + randpart_A + e
  
  }
  
  ## ----Create data frame
  L1_ID <- 1:total_L1
  simu_dataframe <- data.frame(L2_ID, L1_ID, y, predictors)
  colnames(simu_dataframe)[4:5] <- var_names
  
  ## ----Return list containing relevant outputs
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
    alpha = alpha,
    beta = beta,
    total_L1 = total_L1
  )
}
