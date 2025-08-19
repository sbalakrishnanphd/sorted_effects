# NEW SPE function - edited

#' @importFrom Hmisc wtd.quantile
#' @importFrom boot boot
#' @importFrom stats quantile rexp qnorm
#' @importFrom parallel detectCores
#' @importFrom pbapply setpb startpb closepb
#' @export


spe1 <- function(fm, data, method = c("ols", "logit", "probit", "QR"),
                var_type = c("binary", "continuous", "categorical"), var,
                compare, subgroup = NULL, samp_weight = NULL, us = c(1:9)/10,
                alpha = 0.1, taus = c(5:95)/100, b = 500, parallel = FALSE,
                ncores = detectCores(), seed = 1, bc = TRUE,
                boot_type = c("nonpar", "weighted")) {
  # ----- Stopping Condition
  if (alpha >= 1 || alpha <= 0) stop("Please specify a correct size for
                                     hypothesis testing between 0 and 1.")
  # ----- Replace Null samp_weight Specification
  if (is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  samp_weight <- samp_weight/mean(samp_weight) # renormalize
  # ----- Matching Arguments
  method <- match.arg(method)
  var_type <- match.arg(var_type)
  boot_type <- match.arg(boot_type)
  # ----- 1. Call to estimate PE
  output <- suppressWarnings(peestimate(fm, data, samp_weight, var_type, var, compare, method,
                                        subgroup, taus))
  pe_est <- output$pe_est

    # ----- 2. Now get estimated SPE (full and subgroup samples)
  if (method != "QR") {
    spe_est <- wtd.quantile(x=pe_est, weight=NULL, probs=us)
  } else {
    spe_est <- wtd.quantile(pe_est, matrix(samp_weight, ncol = 1,
                                           nrow = nrow(pe_est),
                                           byrow = FALSE), us)
  }
  if (!is.null(subgroup)) {
    pesub_est <- output$pesub_est
    pesub_w <- output$samp_weight_sub
    spesub_est <- wtd.quantile(x=pesub_est, weight=pesub_w, probs=us)
  }
  # -----3. Bootstrap Samples
  #  statistic in one bootstrap for SPE
  #  Accommodate two types of bootstraps:
  #  (1) Nonparametric: draw multinomial weights
  #  (2) Weighted:      draw exponential weights
  
  # set a bootstrap counting variable for the purpose of showing a progress bar
  rep_count <- 1
  
  # The resampling function for weighted bootstrap
  data_rg <- function(data, mle) {
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight
    data$.w <- weight/mean(weight)
    return(data)
  }
  # Resampling function for nonparametric bootstrap
  # Note: nonpar is a special type of weighted bootstrap with multinomial weight
  data_non <- function(data, mle) {
    D        <- data %>% nest(-CLUSTER)
    bs       <- bootstraps(D, times =1)  
    as_tibble(bs$splits[[1]]) %>% arrange(CLUSTER)
    data <- unnest(as_tibble(bs$splits[[1]]))
    data$samp_weight = rep(1, nrow(data))
    data$.w = rep(1, nrow(data))
    return(data)
  }
  # Function that computes bootstrap statistics in each draw
  stat_boot_weight <- function(data) {
    setpb(pb, rep_count)
    rep_count <<- rep_count + 1
    samp_weight = rep(1, nrow(data))
    output_bs <- suppressWarnings(peestimate(fm, data, samp_weight = rep(1, nrow(data)), var_type, var,
                                             compare, method, subgroup, taus))
    est_pe_bs <- output_bs$pe_est
    est_ape_bs <- mean(est_pe_bs) # scalar, length = 1
    if (method != "QR") {
      est_spe_bs <- wtd.quantile(x=est_pe_bs, weight=rep(1, nrow(data)), probs=us)
    } else {
      est_spe_bs <- wtd.quantile(est_pe_bs, matrix(rep(1, nrow(data)), ncol = 1,
                                                   nrow = nrow(pe_est),
                                                   byrow = FALSE), us)
    } # length = length(us)
    if (!is.null(subgroup)) {
      est_pesub_bs <- output_bs$pesub_est
      pesub_w_bs <- output_bs$samp_weight_sub
      est_apesub_bs <- mean(est_pesub_bs) # length = 1
      est_spesub_bs <- wtd.quantile(est_pesub_bs, weight=pesub_w_bs, probs=us) # length(us)
      return(c(est_apesub_bs, est_spesub_bs)) # Each has length: 1+length(us)
    } else {
      return(c(est_ape_bs, est_spe_bs))
    }
  }
  
  # Use boot command
  set.seed(seed)
  if (parallel == FALSE) ncores <- 1
  if (boot_type == "nonpar") {
    data$.w <- samp_weight
    # print a message showing how many cores are used
    cat(paste("Using", ncores, "CPUs now.\n"))
    # set up a progress bar
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, sim = "parametric", ran.gen = data_non, mle=0,
                        statistic = stat_boot_weight, parallel = "multicore", 
                        ncpus = ncores, R = b)
    data$.w <- NULL
    closepb(pb)
  } else if (boot_type == "weighted") {
    data$.w <- samp_weight
    cat(paste("Using", ncores, "CPUs now.\n"))
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = stat_boot_weight,
                        sim = "parametric", ran.gen = data_rg, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    data$.w <- NULL
    closepb(pb)
  }
  # ------ 4. Inference
  ##############################
  ### Full Sample
  ##############################
  # SPE inference
  draws_spe_bs <- result_boot$t[, 2:(length(us) + 1)]
    #print(draws_spe_bs)
  # bundle function implements algorithm 2.1, remark 2.2 and 2.3
  inf_spe <- bundle(draws_spe_bs, spe_est, alpha)
    #print(inf_spe)
  # APE inference
  est_ape <- mean(pe_est)
    #print(est_ape)
  draws_ape_bs <- result_boot$t[, 1]
    #print(draws_ape_bs)
  inf_ape <- ape(draws_ape_bs, est_ape, alpha)
    #print(inf_ape)
  ##############################
  ### Subgroup Sample (if subgroup is not NULL)
  ##############################
  if (!is.null(subgroup)) {
    # SPE
    draws_spesub_bs <- result_boot$t[, 2:(length(us) + 1)]
    inf_spesub <- bundle(draws_spesub_bs, spesub_est, alpha)
    # APE
    est_apesub <- mean(pesub_est)
    draws_apesub_bs <- result_boot$t[, 1]
    inf_apesub <- ape(draws_apesub_bs, est_apesub, alpha)
  }
  # ----- 5. Return Results
  # Depends on whether subgroup is NULL and whether bias correction is wanted
  # The resulting outputs are APE(sub) & SPE(sub) with correpsonding confidence
  # intervals. All are stored in an "inf" bundle list
  if (!is.null(subgroup)) {
    if (bc == TRUE) {
      output <- list(spe = inf_spesub[c(1:3, 7)], ape = inf_apesub[c(1:3, 7)],
                     us = us, alpha = alpha)
    } else {
      output <- list(spe = inf_spesub[4:7], ape = inf_apesub[4:7], us = us,
                     alpha = alpha)
    }
  } else {
    if (bc == TRUE) {
      output <- list(spe = inf_spe[c(1:3, 7)], ape = inf_ape[c(1:3, 7)],
                     us = us, alpha = alpha, individual = pe_est)
    } else {
      output <- list(spe = inf_spe[4:7], ape = inf_ape[4:7], us = us,
                     alpha = alpha, individual = pe_est)
    }
  }
  # claim output as an object of class "spe" so as to use S3
  output <- structure(output, class = "spe")
  return(output)
}

# ------- Auxiliary Functions ----------------
# Implementing algorithm 2.1 and get (bias corrected) estimate and confidence
# bands for sorted effects
bundle <- function(bs, est, alpha) {
  z_bs <- bs - matrix(est, nrow = nrow(bs), ncol = ncol(bs), byrow = TRUE)
  sigma <- (apply(z_bs, 2, quantile, .75, na.rm = TRUE) -
              apply(z_bs, 2, quantile, .25, na.rm = TRUE)) / (qnorm(0.75) -
                                                                qnorm(.25))
  t_hat <- apply(abs(z_bs / matrix(sigma, nrow = nrow(bs), ncol = ncol(bs),
                                   byrow = TRUE)), 1, max)
  crt <- quantile(t_hat, 1 - alpha, na.rm=TRUE)
  # bias-correction
  est_bc <- sort(2*est - apply(bs, 2, mean))
  ubound_est_bc <- sort(est_bc + crt * sigma)
  lbound_est_bc <- sort(est_bc - crt * sigma)
  # uncorrected
  ubound_est <- sort(est + crt * sigma)
  lbound_est <- sort(est - crt * sigma)
  out <- list(est_bc = est_bc, ubound_est_bc = ubound_est_bc,
              lbound_est_bc = lbound_est_bc, est = est, ubound_est = ubound_est,
              lbound_est = lbound_est, sigma_spe = sigma)
  return(out)
}

# Get (biased corrected) average partial effects and confidence intervals
ape <- function(bs, est, alpha) {
  sigma <- (quantile(bs, .75, na.rm = TRUE) -
              quantile(bs, .25, na.rm = TRUE)) / (qnorm(0.75) - qnorm(.25))
  mzs <- abs(bs - est) / sigma
  crt2 <- quantile(mzs, 1 - alpha, na.rm=T)
  # bias-correction
  est_bc <- 2 * est - mean(bs)
  ubound_ape_bc <- est_bc + crt2 * sigma
  lbound_ape_bc <- est_bc - crt2 * sigma
  # uncorrected
  ubound_ape <- est + crt2 * sigma
  lbound_ape <- est - crt2 * sigma
  out <- list(est_bc = est_bc, ubound_est_bc = ubound_ape_bc,
              lbound_est_bc = lbound_ape_bc, est = est, ubound_ape = ubound_ape,
              lbound_ape = lbound_ape, sigma_ape = sigma)
  return(out)
}

# ------- Plotting Function ----------------
# Since plot() is a generic, we write a new plotting method for the SPE class.
# --Plotting (Output: average effect, sorted effect, correspondent confidence bands) --
#' Plot output of \code{\link{spe}} command. The x-axis limits are set to the
#' specified range of percentile index.
#' @param x           Output of \code{\link{spe}} command.
#' @param ylim        y-axis limits. Default is NULL.
#' @param main        Main title of the plot. Defualt is NULL.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is "Percentile Index".
#' @param ylab        y-axis label. Default is "Sorted Effects".
#' @param ...         graphics parameters to be passed to the plotting
#'                    routines.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#' test <- spe(fm = fm, data = mortgage, var = "black", method = "logit",
#' us = c(2:98)/100, b = 50)
#'
#' plot(x = test, main="APE and SPE of Being Black on the prob of
#' Mortgage Denial", sub="Logit Model", ylab="Change in Probability")
#'
#' @importFrom graphics plot polygon lines abline points legend
#' @export
plot.spe <- function(x, ylim = NULL, main = NULL, sub = NULL,
                     xlab = "Percentile Index", ylab = "Sorted Effects",
                     ...) {
  # SPE and confidence bands
  SE <- x$spe
  AE <- x$ape
  us <- x$us
  alpha <- x$alpha
  xlim <- range(us)
  plot(us, as.numeric(SE[[1]]), type = "l", xlim, ylim, log = "", main, sub, xlab, ylab,
       col = 4, lwd = 2)
  polygon(c(us, rev(us)),c(SE[[2]], rev(SE[[3]])), density = 60, border = F,
          col = 'light blue', lty = 1, lwd = 1)
  lines(us, SE[[1]], lwd = 2, col = 4 )
  # APE and CI
  abline(h = AE[[1]], col = 1)
  abline(h = AE[[2]], col = 1, lty = 2)
  abline(h = AE[[3]], col = 1, lty = 2)
  points(c(min(us),.2,.4,.6,.8,max(us)), rep(AE[[1]], 6), col = 1, pch = 15)
  legend(x = "topleft", col = c(4, 1, "light blue", 1), lwd = c(1, 1, 5, 1),
         lty = c(1, 1, 1, 2), pch = c(NA, 15, NA, NA), pt.cex = c(2, 1),
         bty = 'n', legend = c("ITE","ATE", paste0((1 - alpha)*100,"% CI(ITE)"),
                               paste0((1 - alpha)*100,"% CI(ATE)")))
}

# ------- Summary Function ----------------
#' Tabulate the output of \code{\link{spe}} function.
#'
#' The option \code{result} allows user to tabulate either sorted estimates or
#' average estimates. For sorted estimates, the table shows user-specified
#' quantile indices, sorted estimates, standard errors, point-wise confidence
#' intervals, and uniform confidence intervals. For average estimates, the
#' table shows average estiamtes, standard errors, and confidence intervals.
#'
#' @param object   The output of \code{\link{spe}} function.
#' @param result   Whether the user wants to see the sorted or the average
#'                 estimates. Default is \code{sorted}, which shows the
#'                 sorted estimates.
#' @param ...      additional arguments affecting the summary produced.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#' test <- spe(fm = fm, data = mortgage, var = "black", method = "logit",
#' us = c(2:98)/100, b = 50)
#' summary(test)
#' @export
summary.spe <- function(object, result = c("sorted", "average"), ...) {
  spe <- object$spe
  ape <- object$ape
  us <- object$us
  alpha <- object$alpha
  result <- match.arg(result)
  if (result == "sorted") {
    table <- matrix(0, nrow = length(us), ncol = 6)
    table[, 1] <- unlist(spe[1])
    table[, 2] <- unlist(spe[4])
    table[, 3] <- table[, 1] - qnorm(1 - alpha/2)*table[, 2]
    table[, 4] <- table[, 1] + qnorm(1 - alpha/2)*table[, 2]
    table[, 5] <- unlist(spe[3])
    table[, 6] <- unlist(spe[2])
    rownames(table) <- us
    colnames(table) <- c("Est", "SE", paste0((1 - alpha)*100,"% PLB"),
                         paste0((1 - alpha)*100,"% PUB"),
                         paste0((1 - alpha)*100,"% ULB"),
                         paste0((1 - alpha)*100,"% UUB"))
  } else {
    table <- matrix(0, nrow = 1, ncol = 4)
    table[, 1] <- unlist(ape[1])
    table[, 2] <- unlist(ape[4])
    table[, 3] <- unlist(ape[3])
    table[, 4] <- unlist(ape[2])
    rownames(table) <- "APE"
    colnames(table) <- c("Est", "SE", paste0((1 - alpha)*100,"% LB"),
                         paste0((1 - alpha)*100,"% UB"))
  }
  return(table)
}