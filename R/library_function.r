#' Function to generate a TF binding matrix M where i_th row represents
#'  i_th PWM and j_th column represetns j_th peak with M_i, j=1 if the
#'   footprint of i_th PWM overlaps j_th peak and 0 otherwise
#' @param GRpeaks ATAC-seq or DNase-seq peaks in GRanges class where each row
#' represents the location, read counts, GC content and type
#' ("target","background","no") of the peak
#' @param GRmotif Footprint calls from PIQ or CENTIPEDE algorithm in GRanges
#' class where each row represents the location of each PWM occurrence and its
#'  purity score from the algorithm. The footprint calls in the forward strand
#'  and those in the backward strand from the same PWM are not differentiated.
#' @return Returns a TF binding matrix M where i_th row represents i_th PWM
#'  and j_th column represetns j_th peak with M_i, j=1 if the footprint of
#'  i_th PWM overlaps j_th peak and 0 otherwise
#' @author Ahrim Youn
#' @examples
#' # Load in the peak file and footprint calls from PIQ or CENTIPEDE algorithm
#' peak_file <- system.file("extdata", "input_peak_motif.Rdata",
#'  package = "BiFET")
#' load(peak_file)
#' # Generate a TF binding matrix M where i_th row represents
#' # i_th PWM and j_th column represetns j_th peak with M_i, j=1
#' # if the footprint of i_th PWM overlaps j_th peak and 0 otherwise
#' TFbinding.mat <- bindingTF_per_peak(GRpeaks, GRmotif)
#' TFbinding.mat[1:4, 1:4]
#' @export

bindingTF_per_peak <- function(GRpeaks, GRmotif) {

    singles <- split(GRmotif, names(GRmotif))
    singles <- GenomicRanges::GRangesList(singles)
    motifoverlap <- GenomicRanges::findOverlaps(singles, GRpeaks)
    allTF <- names(singles)
    TFbinding.mat <- matrix(0, nrow = length(allTF), ncol = length(GRpeaks))
    rownames(TFbinding.mat) <- allTF
    TFbinding.mat[as.matrix(motifoverlap)] <- 1
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
#' # Load in the peak file and footprint calls from PIQ or CENTIPEDE algorithm
#' peak_file <- system.file("extdata", "input_peak_motif.Rdata",
#'  package = "BiFET")
#' load(peak_file)
#' # Generate a TF binding matrix M where i_th row represents
#' # i_th PWM and j_th column represetns j_th peak with M_i, j=1
#' # if the footprint of i_th PWM overlaps j_th peak and 0 otherwise
#' TFbinding.mat <- bindingTF_per_peak(GRpeaks, GRmotif)
#' # Other input parameter values for the function calculate_enrich_p
#' targetpeak <- which(GRpeaks$peaktype == "target")
#' backgroundpeak <- which(GRpeaks$peaktype == "background")
#' reads <- GRpeaks$reads
#' GCcontent <- GRpeaks$GC
#' result <- calculate_enrich_p(TFbinding.mat, reads,
#'  GCcontent, targetpeak, backgroundpeak)
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
              !(backgroundpeak %in% targetpeak)]], "greater")$p.value > 0.01

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
