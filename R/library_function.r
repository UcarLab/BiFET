#' Function to reformat the output matrix of footprint calls from
#'  PIQ or CENTIPEDE for use in BiFET
#' @param motif an output matrix of footprint calls from PIQ or CENTIPEDE
#'  algorithm where each row represents the location where each PWM
#'   was detected, its ID, purity score from the algorithm and the strand: +/-
#' @return Returns a matrix of footprint calls where the suffix ".RC" in the
#'  id of PWM detected on the reverse strand is removed and redundant counts
#'   of the same PWM on the same location are removed.
#' @author Ahrim Youn
#' @examples
#' # Load in the output matrix of footprint calls from PIQ or CENTIPEDE
#' footprints <- system.file("extdata",
#' "PBMC_PIQ.bed", package = "BiFET")
#' PBMCmotif <- read.delim(footprints, header = FALSE)
#' # Process the output matrix of footprint calls from PIQ
#' # or CENTIPEDE for use in BiFET
#' PBMCmotif <- motif_processing(PBMCmotif)
#' head(PBMCmotif)
#' @export

motif_processing <- function(motif) {
    motif[, 4] <- removing_RC(motif[, 4])
    motif <- motif[, -6]
    motif <- unique(motif)
    return(motif)
}

#' Function to generate a TF binding matrix M where i_th row represents
#'  i_th PWM and j_th column represetns j_th peak with M_i, j=1 if the
#'   footprint of i_th PWM overlaps j_th peak and 0 otherwise
#' @param peaks A matrix of ATAC-seq or DNase-seq peak whose first column
#'  is the chromosome number, the second and third column is the starting
#'   and ending position of the peak.
#' @param motif A matrix of footprint calls outputted from
#' motif_processing function
#' @return Returns a TF binding matrix M where i_th row represents i_th PWM
#'  and j_th column represetns j_th peak with M_i, j=1 if the footprint of
#'  i_th PWM overlaps j_th peak and 0 otherwise
#' @seealso \code{\link{motif_processing}}
#' @author Ahrim Youn
#' @examples
#' peak_file <- system.file("extdata",
#' "islet_PBMC_consensus_peak.Rdata", package = "BiFET")
#' load(peak_file)
#' # Load in the output matrix of footprint calls from PIQ or CENTIPEDE
#' footprints <- system.file("extdata",
#' "PBMC_PIQ.bed", package = "BiFET")
#' PBMCmotif <- read.delim(footprints, header = FALSE)
#' # Process the output matrix of footprint calls from PIQ
#' # or CENTIPEDE for use in BiFET
#' PBMCmotif <- motif_processing(PBMCmotif)
#' # Generate a TF binding matrix M where i_th row represents
#' # i_th PWM and j_th column represetns j_th peak with M_i, j=1
#' # if the footprint of i_th PWM overlaps j_th peak and 0 otherwise
#' TFbinding.mat <- bindingTF_per_peak(peaks, PBMCmotif)
#' head(TFbinding.mat)
#' @export

bindingTF_per_peak <- function(peaks, motif) {
    alltf <- unique(motif[, 4])
    TFbinding.mat <- matrix(0, nrow = length(alltf), ncol = nrow(peaks))
    rownames(TFbinding.mat) <- alltf

    for (i in 1:nrow(peaks)) {
        print(paste("Peak ", i, " out of ", nrow(peaks), " peaks", sep = ""))
        ipeak.motif <- motif[as.character(peaks[i, 1]) ==
                               as.character(motif[, 1]), ]
        overlap <- peaks[i, 2] < as.numeric(ipeak.motif[, 3]) & peaks[i, 3] >
          as.numeric(ipeak.motif[, 2])
        bindingtf <- unique(as.character(ipeak.motif[overlap, 4]))
        TFbinding.mat[bindingtf, i] <- 1
    }
    TFbinding.mat <- TFbinding.mat[rowSums(TFbinding.mat) > 0, ]
    return(TFbinding.mat)
}


#' Function to calculate p-value testing if footprints of a TF are
#'  over-represented in the target set of peaks compared to the background
#'   set of peaks correcting for the bias arising from the imbalance of
#'    GC-content and read counts between target and background set
#' @param TFbinding.mat A TF binding matrix M outputted from the function
#'  bindingTF_per_peak where i_th row represents i_th PWM and j_th column
#'   represetns j_th peak with M_i,j=1 if the footprint of i_th PWM
#'    overlaps j_th peak and 0 otherwise
#' @param reads A vector of read counts per peak
#' @param GCcontent A vector of GC contents per peak
#' @param targetpeak A vector of peak IDs belonging to the target set
#' @param backgroundpeak A vector of peak IDs belonging to the background set
#' @return Returns a list of parameter alpha_k, p values from BiFET algorithm
#'  and p values from the hypergeometric test
#' @seealso \code{\link{bindingTF_per_peak}}
#' @author Ahrim Youn
#' @examples
#' # Load in a peak file, a matrix of ATAC-seq or DNase-seq peak whose
#' # first column is the chromosome number, the second and third column
#' # are the starting and ending position of the peak, the fourth and
#' # fifth column are the read counts and GC content of the peak and
#' # the sixth column is the type of the peak: "target","background","no"
#'
#' # GC content of the peak was obtained using peak annotation program
#' # "annotatePeaks.pl" from the HOMER software
#' \dontrun{
#' peak_file <- system.file("extdata",
#' "islet_PBMC_consensus_peak.Rdata", package = "BiFET")
#' load(peak_file)
#'
#' # Load in the output matrix of footprint calls from PIQ or CENTIPEDE
#' footprints <- system.file("extdata",
#' "PBMC_PIQ.bed", package = "BiFET")
#' PBMCmotif <- read.delim(footprints, header = FALSE)
#' # Process the output matrix of footprint calls from PIQ
#' # or CENTIPEDE for use in BiFET
#' PBMCmotif <- motif_processing(PBMCmotif)
#' # Generate a TF binding matrix M where i_th row represents
#' # i_th PWM and j_th column represetns j_th peak with M_i, j=1
#' # if the footprint of i_th PWM overlaps j_th peak and 0 otherwise
#' TFbinding.mat <- bindingTF_per_peak(peaks, PBMCmotif)
#' # Other input parameter values for the function calculate_enrich_p
#' targetpeak <- which(peaks[, "peaktype"] == "target")
#' backgroundpeak <- which(peaks[, "peaktype"] == "background")
#' reads <- peaks[, "reads"]
#' GCcontent <- peaks[, "GCcontent"]
#' result <- calculate_enrich_p(TFbinding.mat, reads,
#'  GCcontent, targetpeak, backgroundpeak)
#'  }
#' @export

calculate_enrich_p <- function(TFbinding.mat, reads,
                               GCcontent, targetpeak, backgroundpeak) {
    backgroundpeak <- c(targetpeak, backgroundpeak)
    upperlimit <- stats::quantile(reads, 0.999)
    reads[reads > upperlimit] <- upperlimit

    q <- rep(NA, nrow(TFbinding.mat))
    alpha <- matrix(NA, nrow = nrow(TFbinding.mat), ncol = 2)
    p.mat <- matrix(NA, nrow(TFbinding.mat), ncol = 2)
    rownames(p.mat) <- rownames(TFbinding.mat)

    for (i in (1:nrow(TFbinding.mat))) {
        print(paste("PWM ", i, " out of ",
                    nrow(TFbinding.mat), " PWMs", sep = ""))
        res <- calculate_enrich_p_per_TF(rbind(TFbinding.mat[i, ]),
                    reads, GCcontent, targetpeak, backgroundpeak)
        q[i] <- res[[1]]
        alpha[i, ] <- res[[2]]
        p.mat[i, ] <- res[[3]]
    }

    res <- list(alpha_k = 1 / alpha[, 1],
                p.bifet = p.mat[, 1], p.hyper = p.mat[, 2])
    return(res)
}


calculate_enrich_p_per_TF <- function(TFbinding.mat, reads,
                        GCcontent, targetpeak, backgroundpeak) {
    GCbias.cutoff <- 0.05
    GCcontent.same <- stats::t.test(GCcontent[targetpeak],
            GCcontent[backgroundpeak[
              !(backgroundpeak %in% targetpeak)]])$p.value > 0.01

    model <- stats::glm(TFbinding.mat[, backgroundpeak]~reads[backgroundpeak] +
                   GCcontent[backgroundpeak],
                   family = stats::binomial(link = "logit"))
    GCcoef <- summary(model)$coef[3, ]

    if ( !GCcontent.same & GCcoef[1] > 0 & GCcoef[4] < GCbias.cutoff) {
        res <- main_optim(rbind(TFbinding.mat[, backgroundpeak]),
                  reads[backgroundpeak], GCcontent[backgroundpeak], c(300, 3))
        q <- res[[1]]
        alpha <- res[[2]]
        p.mat <- enrich_p(reads, GCcontent, TFbinding.mat,
                          q, alpha, targetpeak, backgroundpeak)

    } else {
        res <- main_optim(rbind(TFbinding.mat[, backgroundpeak]),
                reads[backgroundpeak],
                rep(Inf, length(backgroundpeak)), c(300, 3))
        q <- res[[1]]
        alpha <- res[[2]]
        p.mat <- enrich_p(reads, rep(Inf, length(GCcontent)),
                TFbinding.mat, q, alpha, targetpeak, backgroundpeak)
        alpha[2] <- 10 ^ (-10)
    }
    res <- list(q, alpha, p.mat)
    return(res)
}

main_optim <- function(TFbinding.mat, reads, GCcontent, alphainit) {
  total_lik_old <- -Inf
  total_lik_new <- -10 ^ 10
  alpha <- alphainit
  count <- 0
  while (total_lik_old - total_lik_new < -10 ^ (-10) *
         abs(total_lik_new) & count < 100) {
      count <- count + 1
      total_lik_old <- total_lik_new

      init <- min(sum(TFbinding.mat) / sum(f(reads, alpha, GCcontent)),
                  1 / max(f(reads, alpha, GCcontent)) - 10 ^ (-6))

      temp <- stats::optim(init, fn = loglik_per_TF, gr = deriv_per_TF,
              lower = 10 ^ (-16),
              upper = 1 / max(f(reads, alpha, GCcontent)) - 10 ^ (-10),
              method = "L-BFGS-B", reads = reads, alpha = alpha,
              GCcontent = GCcontent, TFbinding.mat = TFbinding.mat,
              control = list(fnscale = -1))

      q <- temp$par

      temp <- stats::optim(par = alpha[1], fn = total_lik1, gr = total_deriv1,
              lower = 10 ^ (-5) + max(0, (reads / (log(q * f0(GCcontent,
                    alpha[2]) + 2) - log(abs(q * f0(GCcontent, alpha[2]) -
                        2))))[q * f0(GCcontent, alpha[2]) > 2]),
              upper = Inf, method = "L-BFGS-B", alpha2 = alpha[2], qi = q,
              reads = reads, GCcontent = GCcontent,
              TFbinding.mat = TFbinding.mat, control = list(fnscale = -1))
      alpha[1] <- temp$par

      if (GCcontent[1] != Inf) {
          temp <- stats::optim(par = alpha[2],
                  fn = total_lik2, gr = total_deriv2,
                  lower = 10 ^ (-5) + max(0, (GCcontent / (log(q *
                      f0(reads, alpha[1]) + 2) - log(abs(q *
                        f0(reads, alpha[1]) - 2))))[q *
                          f0(reads, alpha[1]) > 2]),
                  upper = Inf, method = "L-BFGS-B", alpha1 = alpha[1],
                  qi = q, reads = reads, GCcontent = GCcontent,
                  TFbinding.mat = TFbinding.mat, control = list(fnscale = -1))
          alpha[2] <- temp$par
      }

      total_lik_new <- total_lik(alpha, q, reads, GCcontent, TFbinding.mat)

  }
  return(list(q, alpha, total_lik_new))
}

enrich_p <- function(reads, GCcontent, TFbinding.mat,
                     q, alpha, targetpeak, backgroundpeak) {
    pvalue.bifet <- pvalue.hyper <- rep(NA, nrow(TFbinding.mat))
    for (ii in 1:nrow(TFbinding.mat)) {
        no.footprint <- TFbinding.mat[ii, ]
        pp <- q[ii] * f(reads[targetpeak], alpha, GCcontent[targetpeak])
        pvalue.bifet[ii] <- 1 -
          poibin::ppoibin(kk = sum(no.footprint[targetpeak]) - 1, pp = pp)

        pick <- sum(no.footprint[backgroundpeak])
        whitepick <- sum(no.footprint[targetpeak]) - 1
        white <- length(targetpeak)
        black <- length(backgroundpeak) - length(targetpeak)

        pvalue.hyper[ii] <- 1 - stats::phyper(whitepick, white, black, pick)
    }
    res <- cbind(pvalue.bifet, pvalue.hyper)
    rownames(res) <- rownames(TFbinding.mat)
    return(res)
}


removing_RC <- function(TF) {
  TF <- as.vector(as.matrix(as.data.frame(strsplit(as.character(TF), ".RC"))))
  return(TF)
}

f0 <- function(x, y) {
    return(1 / (1 + exp(-x / y)) - 1 / 2)
}

df0 <- function(x, y) {
    return(-x / y ^ 2 * exp(-x / y) / (1 + exp(-x / y)) ^ 2)
}

f <- function(reads, alpha, GCcontent) {
    return(f0(reads, alpha[1]) * f0(GCcontent, alpha[2]))
}

df <- function(reads, alpha, GCcontent) {
    return(cbind(df0(reads, alpha[1]) * f0(GCcontent, alpha[2]),
        f0(reads, alpha[1]) * df0(GCcontent, alpha[2])))
}

loglik_per_TF <- function(qi, reads, alpha, GCcontent, TFbinding.mat) {
    qf <- qi * f(reads, alpha, GCcontent)
    return(sum(log(qf) * TFbinding.mat + log(1 - qf) * (1 - TFbinding.mat)))
}
deriv_per_TF <- function(qi, reads, alpha, GCcontent, TFbinding.mat) {
  f <- f(reads, alpha, GCcontent)
  qf <- qi * f
  return(sum(1 / qi * TFbinding.mat - f / (1 - qf) * (1 - TFbinding.mat)))
}

total_lik <- function(alpha, qi, reads, GCcontent, TFbinding.mat) {
    return(loglik_per_TF(qi, reads, alpha, GCcontent, TFbinding.mat))
}

total_lik1 <- function(alpha1, alpha2, qi, reads, GCcontent, TFbinding.mat) {
    alpha <- c(alpha1, alpha2)
    return(loglik_per_TF(qi, reads, alpha, GCcontent, TFbinding.mat))
}

total_lik2 <- function(alpha2, alpha1, qi, reads, GCcontent, TFbinding.mat) {
    alpha <- c(alpha1, alpha2)
    return(loglik_per_TF(qi, reads, alpha, GCcontent, TFbinding.mat))
}

total_deriv <- function(alpha, qi, reads, GCcontent, TFbinding.mat) {
    qf <- qi * f(reads, alpha, GCcontent)
    qdf <- qi * df(reads, alpha, GCcontent)
    return(colSums(qdf / qf * as.vector(TFbinding.mat) -
          qdf / (1 - qf) * (1 - as.vector(TFbinding.mat))))
}

total_deriv1 <- function(alpha1, alpha2, qi,
                         reads, GCcontent, TFbinding.mat) {
    alpha <- c(alpha1, alpha2)
    qf <- qi * f(reads, alpha, GCcontent)
    qdf <- qi * df(reads, alpha, GCcontent)

    return( colSums(qdf / qf * as.vector(TFbinding.mat) -
          qdf / (1 - qf) * (1 - as.vector(TFbinding.mat)))[1])
}

total_deriv2 <- function(alpha2, alpha1, qi, reads, GCcontent, TFbinding.mat) {
    alpha <- c(alpha1, alpha2)
    qf <- qi * f(reads, alpha, GCcontent)
    qdf <- qi * df(reads, alpha, GCcontent)

    return(colSums(qdf / qf * as.vector(TFbinding.mat) -
          qdf / (1 - qf) * (1 - as.vector(TFbinding.mat)))[2])
}
