#' Simulate Datasets with User-Defined Patterns of Within-Individual Variation
#' 
#' Function to create simulated dataset given user-defined parameters and data
#' structure for growth curve model with optional random effect for L2 variance
#' 
#' @param sample_size List of named vectors specifying sample size at level 2 (named \code{n2}) and level 1 (named \code{n1}).
#' @param var_names List of named character vectors indicating names for intercept (\code{x0}) and time variable (\code{x1}),
#' e.g. \code{list(x0 = "cons", x1 = "age")}.
#' @param x1_range Vector with min and max values for time variable (in that order: e.g. \code{c(12, 23)} to range from 12 to 23).
#' @param design_matrices List of five named lists, with each of the latter relating to a design matrix. These five lists are \code{x_A}
#' (corresponding to X matrix (fixed part) of model for mean of y), \code{z_A} (Z matrix (random part) of model for mean of y),
#' \code{x_B} (X matrix (fixed part) of model for level 1 variance function), \code{z_B} (Z matrix (random part) of model for
#' level 1 variance function) and \code{x_C} (X matrix (fixed part) for model of distal outcome). Each of the first four lists (i.e. \code{x_A},
#' \code{z_A}, \code{x_B} and \code{z_B}) to contain two named logical vectors (\code{x0} for intercept, \code{x1} for time), indicating whether each variable to be included
#' in design matrices or not; \code{TRUE} indicates variable to be included, \code{FALSE} otherwise. The last list (\code{x_C}) to contain up to four named logical vectors indicating
#' whether constant (\code{x0}), and u2j terms (\code{x1}, \code{x2}, \code{x3}) to be included in model for distal outcome or not.
#' E.g. \code{list(x_A = list(x0 = TRUE, x1 = TRUE), z_A = list(x0 = TRUE, x1 = TRUE), x_B = list(x0 = TRUE, x1 = TRUE), z_B = list(x0 = TRUE, x1 = FALSE))} for each variable
#' in each design matrix, bar time variable excluded from \code{z_B}. Note if \code{gamma = NULL}, \code{x_C} is set to \code{NULL}
#' (i.e. no need to include \code{x_C} in \code{design_matrices} if do not wish to simulate distal outcome).
#' @param beta Vector containing coefficients of fixed effects for model for mean of y, in order of x0, x1.
#' @param alpha Vector containing coefficients of fixed effects for level 1 variance function, in order of x0, x1.
#' @param gamma Vector containing coefficients of fixed effects for model of distal outcome, in order of x0, x1, x2, x3. Defaults to \code{NULL}.
#' @param distal_type Character vector specifying type for distal outcome: either \code{"Continuous"} or \code{"Binary"};
#' defaults to \code{NULL}, i.e. assumes no distal outcome required.
#' @param sigma2e_distal Residual variance of distal outcome (applies when \code{distal_type = "Continuous"} only).
#' @param sigma2u Level 2 covariance matrix. E.g. \code{matrix(c(93, 12, 2.5, 12, 3.8, 0.8, 2.5, 0.8, 0.4), nrow = 3)} 
#' would correspond to the following: sigma2u00 (\code{93}), sigma2u01 (\code{12}), sigma2u02 (\code{2.5}), sigma2u01 (\code{12}), sigma2u11 (\code{3.8}),
#' sigma2u12 (\code{0.8}), sigma2u02 (\code{2.5}), sigma2u12 (\code{0.8}), sigma2u22 (\code{0.4}).
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
#'   var_names = list(x0 = "cons", x1 = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = list(x0 = TRUE, x1 = TRUE),
#'     z_A = list(x0 = TRUE, x1 = TRUE),
#'     x_B = list(x0 = TRUE, x1 = TRUE),
#'     z_B = list(x0 = TRUE, x1 = FALSE)
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
#' ## Level 2 random effect in level 1 variance function
#' ## but covariance of this term set to 0 in level 2 covariance matrix:
#' simu_u2j_no_covary <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = list(x0 = "cons", x1 = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = list(x0 = TRUE, x1 = TRUE),
#'     z_A = list(x0 = TRUE, x1 = TRUE),
#'     x_B = list(x0 = TRUE, x1 = TRUE),
#'     z_B = list(x0 = TRUE, x1 = FALSE)
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
#'   log_sigma2e = TRUE,
#'   linear_sigma2e_attempts = 10,
#'   random_seed = 1
#' )
#' 
#' 
#' ## No additional level 2 random effect:
#' simu_no_u2j <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = list(x0 = "cons", x1 = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = list(x0 = TRUE, x1 = TRUE),
#'     z_A = list(x0 = TRUE, x1 = TRUE),
#'     x_B = list(x0 = TRUE, x1 = TRUE),
#'     z_B = list(x0 = FALSE, x1 = FALSE)
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
#' simu_u2j_covary_linearL1varfun <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = list(x0 = "cons", x1 = "age"),
#'   x1_range = c(11.25, 13.25),
#'   design_matrices = list(
#'     x_A = list(x0 = TRUE, x1 = TRUE),
#'     z_A = list(x0 = TRUE, x1 = TRUE),
#'     x_B = list(x0 = TRUE, x1 = TRUE),
#'     z_B = list(x0 = TRUE, x1 = FALSE)
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
#'   log_sigma2e = FALSE,
#'   linear_sigma2e_attempts = 100,
#'   random_seed = 1
#' )
#' 
#' ## As simu_u2j_covary (above), but with distal outcome;
#' ## note x_C in design_matrices, and gamma:
#' simu_u2j_covary_distal <- windsimu(
#'   sample_size = list(n2 = 1000, n1 = 9),
#'   var_names = list(x0 = "cons", x1 = "age"),
#'   x1_range = c(-1, 1),
#'   design_matrices = list(
#'     x_A = list(x0 = TRUE, x1 = TRUE),
#'     z_A = list(x0 = TRUE, x1 = TRUE),
#'     x_B = list(x0 = TRUE, x1 = TRUE),
#'     z_B = list(x0 = TRUE, x1 = FALSE),
#'     x_C = list(x0 = TRUE, x1 = TRUE, x2 = TRUE, x3 = TRUE)
#'   ),
#'   beta = c(150, 6.5),
#'   alpha = c(-.95, 0.48),
#'   gamma = c(0.5, 0.2, 0.7, 1.4),
#'   distal_type = "Continuous",
#'   sigma2e_distal = 1,
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
#' }
#'
#' @export
windsimu <- function(sample_size,
                     var_names,
                     x1_range,
                     design_matrices,
                     beta,
                     alpha,
                     gamma = NULL,
                     distal_type = NULL,
                     sigma2e_distal = NULL,
                     sigma2u,
                     log_sigma2e = TRUE,
                     linear_sigma2e_attempts = 10,
                     random_seed) {

  set.seed(random_seed)
  
  # return error if incorrect number of sample sizes ---------------------------  
  if (length(sample_size) != 2){
    stop("Incorrect number of sample sizes specified: see sample_size via ?windsimu", call. = FALSE)
  }

  # ensure elements of design matrices in expected order ---------------------------  

  # ensure x0, x1 are in that order for submodel A and B
  # (NB not necessary to order sample_size as refer to elements of that by name only,
  # likewise x_A, z_A. etc., in design_matrices)
  var_names <- unlist(var_names[c("x0","x1")])
  design_matrices[c("x_A", "z_A", "x_B", "z_B")] <- lapply(design_matrices[c("x_A", "z_A", "x_B", "z_B")], function(x) x[c("x0", "x1")])

  # for submodel C, ensure x_C is NA if no gamma specified; else ensure x0, x1, x2, x3 are in that order
  if(is.null(gamma)){
    design_matrices["x_C"] <- NA
    } else {
      if (is.null(distal_type)){
        stop("Need to specify response type (Binary or Continuous) for distal outcome, via distal_type argument", call. = FALSE)
      }
      design_matrices["x_C"] <- lapply(design_matrices["x_C"], function(x) x[c("x0", "x1", "x2", "x3")])
    }

  # unlist ---------------------------
  design_matrices <- lapply(design_matrices, function(x) unlist(x))
  
  # create vectors ---------------------------
  total_L1 <- sample_size$n2 * sample_size$n1
  L2_ID <- rep(1:sample_size$n2, each = sample_size$n1)
  time_var <- rep(seq(from = x1_range[1], to = x1_range[2], length.out = sample_size$n1), times = sample_size$n2)
  how_many_variables_x_matrix_A <- sum(design_matrices$x_A)
  how_many_variables_z_matrix_A <- sum(design_matrices$z_A)
  how_many_variables_x_matrix_B <- sum(design_matrices$x_B)
  how_many_variables_z_matrix_B <- sum(design_matrices$z_B)

  # return error if no variables in z matrix ---------------------------
  if (how_many_variables_z_matrix_A + how_many_variables_z_matrix_B == 0) {
    stop("Z matrices empty: expects at least one variable (specify via TRUE)", call. = FALSE)
  }
  
  # create predictors ---------------------------
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

  # Create sigma2_e, u, randpart_B ---------------------------
  # This is all wrapped up in a function (called later),
  # as may need to run whole thing repeatedly to get non-negative
  # sigma2_e if user has specified identity link for level 1
  # variance function (i.e. if log_sigma2e = FALSE)

  create_sigma2_e <- function(n2,
                              how_many_variables_z_matrix_A,
                              how_many_variables_z_matrix_B,
                              sigma2u,
                              z_matrix_B = NULL,
                              fixpart_B) {
      u <- MASS::mvrnorm(n = n2,
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
    
  } # end create_sigma2_e()

  # Generate sigma_e when log link for level 1 variance function ---------------------------
    if (log_sigma2e == TRUE) {
      sigma2_e_function_call <- create_sigma2_e(
        n2 = sample_size$n2,
        how_many_variables_z_matrix_A = how_many_variables_z_matrix_A,
        how_many_variables_z_matrix_B = how_many_variables_z_matrix_B,
        sigma2u = sigma2u,
        z_matrix_B = z_matrix_B,
        fixpart_B = fixpart_B
        )
      sigma_e <- sqrt(exp(sigma2_e_function_call$sigma2_e))
    
      } else {
      
      # Generate sigma_e when identity link for level 1 variance function ---------------------------
      for (i in linear_sigma2e_attempts) {
        sigma2_e_function_call <- create_sigma2_e(
          n2 = sample_size$n2,
          how_many_variables_z_matrix_A = how_many_variables_z_matrix_A,
          how_many_variables_z_matrix_B = how_many_variables_z_matrix_B,
          sigma2u = sigma2u,
          z_matrix_B = z_matrix_B,
          fixpart_B = fixpart_B
        )

        if (all(sigma2_e_function_call$sigma2_e >= 0)) {
          break
        }
        if (i == linear_sigma2e_attempts)
          stop(
            "No viable simulations of sigma2e achieved within max number of iterations specified."
          )
      }
      sigma_e <- sqrt(sigma2_e_function_call$sigma2_e)
    }

  # Create level 1 residuals ---------------------------
  
  e <- rnorm(n = total_L1, mean = 0, sd = sigma_e)

  # Create randpart_A; generate y ---------------------------
  
  randpart_A <- rowSums(z_matrix_A * sigma2_e_function_call$u[L2_ID, c("u0j", "u1j")])
  
  y <- fixpart_A + randpart_A + e

  # Create L1 data frame ---------------------------
  L1_ID <- 1:total_L1
  simu_dataframe <- data.frame(L2_ID, L1_ID, y, predictors)
  colnames(simu_dataframe)[4:5] <- var_names

  # Create distal outcome; add to L1 data frame and create L2 data frame ---------------------------
  if(!is.null(gamma)){
      L2_cons <- rep(1, times = sample_size$n2)
      predictors_C <- cbind(L2_cons, sigma2_e_function_call$u) # append u to cons
      x_matrix_C <- predictors_C[, design_matrices$x_C] # only keep those predictors user requested
      # multiply matrix with vector (gamma) (row-wise, via sweep()),
      # then add up (row-wise, via rowSums()):
      fixpart_C <- rowSums(
        sweep(
          x_matrix_C, MARGIN = 2, gamma, `*`
          )
        )
      if (distal_type == "Binary"){
        pr = exp(fixpart_C) / (1 + exp(fixpart_C))
        distal_y = rbinom(sample_size$n2, 1, pr)
      } else { # "Continuous" distal_type
        randpart_C <- rnorm(n = sample_size$n2, mean = 0, sd = sqrt(sigma2e_distal))
        distal_y <- fixpart_C + randpart_C
      }
      simu_dataframe <- cbind(simu_dataframe, distal_y = distal_y[L2_ID])
      L2_ID_B <- unique(L2_ID)
      simu_dataframe_distal <- data.frame(cbind(L2_ID = L2_ID_B, distal_y, cons = L2_cons, sigma2_e_function_call$u))
      
  } else {
    # if user does not require distal outcome:
    x_matrix_C <- NULL
    fixpart_C <- NULL
    randpart_C <- NULL # what about if distal_type == "Binary"?
    simu_dataframe_distal <- NULL
  }
  
  # Return list containing relevant outputs ---------------------------
  list(
    simu_dataframe = simu_dataframe,
    simu_dataframe_distal = simu_dataframe_distal,
    sample_size = sample_size,
    total_L1 = total_L1,
    x_matrix_A = x_matrix_A,
    z_matrix_A = z_matrix_A,
    x_matrix_B = x_matrix_B,
    z_matrix_B = z_matrix_B,
    x_matrix_C = x_matrix_C,
    fixpart_A = fixpart_A,
    randpart_A = randpart_A,
    fixpart_B = fixpart_B,
    randpart_B = sigma2_e_function_call$randpart_B,
    fixpart_C = fixpart_C,
    randpart_C = randpart_C,
    u = sigma2_e_function_call$u,
    sigma2_e = sigma2_e_function_call$sigma2_e,
    sigma_e = sigma_e,
    beta = beta,
    alpha = alpha,
    gamma = gamma,
    sigma2u = sigma2u,
    log_sigma2e = log_sigma2e,
    linear_sigma2e_attempts = linear_sigma2e_attempts
  )
  
}
