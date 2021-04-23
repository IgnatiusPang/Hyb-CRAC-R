## This code implements Gaussian process analysis, in order to identify
## transcripts that significantly change their binding profiles in response to
## a nutrient shift.
##
## @authors Alina Selega, Sander Granneman, Guido Sanguinetti

## Functions ------------------------------------------------------------------

## For internal use
## This function computes likelihood for multiple replicates
likelihood_replicates <- function(X, y, M, precision, sigma2) {
    N <- length(X)

    ## Compute log marginal likelihood p(y) for treatment GP
    exponent1 <- 0
    exponent2 <- 0
    C <- solve(precision)
    for (j in 1:M) {
        exponent1 <- exponent1 + t(y)[(1 + (j-1)*N):(N + (j-1)*N)] %*%
          y[(1 + (j-1)*N):(N + (j-1)*N)]
        for (k in 1:M) {
            exponent2 <- exponent2 + t(y)[(1 + (j-1)*N):(N + (j-1)*N)] %*% C %*%
              y[(1 + (k-1)*N):(N + (k-1)*N)]
        }
    }
    exponent1 <- (-1/(2*sigma2)) * exponent1
    exponent2 <- 1/(2*(sigma2^2)) * exponent2
    exponent <- exponent1 + exponent2

    return(exponent)
}

## For internal use
## This function computes the covariance matrix for GP (with RBF kernel)
sigma_f_noise_covariance <- function(X, amplitude, lengthscale) {
    sigma_f <- matrix(, nrow = length(X), ncol = length(X))
    for (i in 1:dim(sigma_f)[1]) {
        for (j in 1:dim(sigma_f)[2]) {
            sigma_f[i,j] <- amplitude^2 * exp(-((X[i] - X[j])^2) /
                                                (2 * lengthscale^2))
        }
    }
    return(sigma_f)
}

## For internal use
## This function computes the precision matrix of GP by completing the square
precision_matrix <- function(sigma_f, sigma2, N, M) {
    return(solve(sigma_f) + (M/sigma2) * diag(N))
}

## For internal use
## This function computes constants
constants <- function(sigma2, N, M, precision) {
    ## Constant of p(y|f)
    const2 <- 1 / (sqrt(2*pi * sigma2))^(N*M)
    ## Constant of the completed square
    const3 <- 1 / sqrt((2*pi)^N * det(solve(precision)))
    return(list("const2" = const2, "const3" = const3))
}

## This function computes marginal log-likelihoods of:
## the control observations;
## the treatment observations;
## all observations,
## marginalising out a hidden Gaussian process with a squared-exponential
## covariance function.
##
## For the first two log-likelihoods, two distinct Gaussian processes are fitted
## to control and treatment observations correspondingly, representing two
## different underlying processes generating binding profiles in control and in
## treatment (alternative hypothesis).
##
## For the last log-likelihood, one Gaussian process is fitted to all
## observations, representing a single underlying process generating
## observations in control and in treatment (null hypothesis: no change in
## response to nutrient shift).
##
## Parameters include:
## @param matrix_final
##    A matrix of binding profiles. Rows correspond to transcripts, columns
##    correspond to time points. Replicates for a transcript should be
##    concatenated column-wise. Control replicates should come first, followed
##    by treatment replicates.
## @param lengthscale
##    The lengthscale hyperparameter of the squared-exponential covariance
##    function used to model covariance between unobserved function values.
## @param X
##    Vector of time points.
## @param num_controls
##    Number of control replicates.
## @param num_treatments
##    Number of treatment replicates.
##
## @value
##    The functon returns a list of three vectors, each containing the
##    marginal log-likelihoods (or data evidence) of:
##    control observations,
##    treatment observations,
##    and all observations for each transcript.
##
marginal_loglik <- function(matrix_final, lengthscale, X, num_controls,
                            num_treatments) {

    ## Matrix for control GP log-likelihoods
    loglik_control <- array(dim(matrix_final)[1])

    ## Matrix for treatment GP log-likelihoods
    loglik_stress <- array(dim(matrix_final)[1])

    ## Matrix for M0 model GP log-likelihoods
    loglik_m0 <- array(dim(matrix_final)[1])

    ## Number of time points
    N <- length(X)

    for (transcript in 1:dim(matrix_final)[1]) {

        ## Noise variance set as var(control)
        sigma2 <- var(matrix_final[transcript, 1:N])

        ## RBF kernel parameter set as
        ## amplitude^2 = 0.5 * var(all time-series)
        amplitude <- 0.71 * sd(matrix_final[transcript, ])

        ## Compute the covariance matrix for GP (with RBF kernel)
        sigma_f <- sigma_f_noise_covariance(X, amplitude, lengthscale)

        ## Constant of unobserved p(f)
        const1 <- 1 / sqrt((2*pi)^N * det(sigma_f))

        ## Complete the square for control GP: precision matrix
        precision_control <- precision_matrix(sigma_f, sigma2, N, num_controls)

        ## Compute constants
        const2 <- constants(sigma2, N, num_controls, precision_control)$const2
        const3 <- constants(sigma2, N, num_controls, precision_control)$const3

        ## Extract control observations y
        y <- matrix_final[transcript, 1:N]
        dim(y) <- c(length(y), 1)

        ## Compute marginal likelihood p(y) for control GP
        loglik_control[transcript] <- log(const1) + log(const2) +
                                      log(const3) +
                                      likelihood_replicates(X, y, num_controls,
                                                            precision_control,
                                                            sigma2)


        ## Complete the square for treatment GP: precision matrix
        precision_stress <- precision_matrix(sigma_f, sigma2, N, num_treatments)

        ## Compute constants
        const2 <- constants(sigma2, N, num_treatments, precision_stress)$const2
        const3 <- constants(sigma2, N, num_treatments, precision_stress)$const3

        ## Extract treatment observations y
        y <- matrix_final[transcript, (N+1):(dim(matrix_final)[2])]
        dim(y) <- c(length(y), 1)

        ## Compute marginal likelihood p(y) for treatment GP
        loglik_stress[transcript] <- log(const1) + log(const2) +
                                     log(const3) +
                                     likelihood_replicates(X, y, num_treatments,
                                                           precision_stress,
                                                           sigma2)


        ## Complete the square for M0 model GP: precision matrix
        precision_m0 <- precision_matrix(sigma_f, sigma2, N,
                                         num_controls + num_treatments)

        ## Compute constants
        const2 <- constants(sigma2, N, num_controls + num_treatments,
                            precision_m0)$const2
        const3 <- constants(sigma2, N, num_controls + num_treatments,
                            precision_m0)$const3

        ## Extract all observations y
        y <- matrix_final[transcript, 1:dim(matrix_final)[2]]
        dim(y) <- c(length(y), 1)

        ## Compute marginal likelihood p(y) for M0 model GP
        loglik_m0[transcript] <- log(const1) + log(const2) +
                                 log(const3) +
                                 likelihood_replicates(X, y,
                                                      num_controls + num_treatments,
                                                      precision_m0, sigma2)

        print(transcript)
    }

    return(list("loglik_control" = loglik_control,
                "loglik_stress" = loglik_stress, "loglik_m0" = loglik_m0))
}

## ----------------------------------------------------------------------------

## Set the lengthscale hyperparameter of the RBF kernel
## (default value: 1 minute)
lengthscale <- 1

## Specify time points
X <- c(0, 1, 2, 4, 8, 14, 20)

## Read in data matrix of binding profiles
bindingProfiles <- read.table('/path/')

## Specify numbers of control and treatment replicates
numControls <- 1
numTreatments <- 2

## Compute log-likelihoods
loglikelihoods <- marginal_loglik(bindingProfiles, lengthscale, X, numControls,
                                  numTreatments)

## Compute Bayes Factors as the sum of log-likelihoods of independent processes
## minus the log-likelihood of the null hypothesis model
bayesFactors <- (loglikelihoods$loglik_control + loglikelihoods$loglik_stress) -
                 loglikelihoods$loglik_m0

## Select those transcripts that have Bayes Factor >= 10
highBF <- which(bayesFactors >= log(10))
