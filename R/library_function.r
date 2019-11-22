#' Function to calculate p-value testing if footprints of a TF are
#'  over-represented in the target set of peaks compared to the background
#'   set of peaks correcting for the bias arising from the imbalance of
#'    GC-content and read counts between target and background set
#' @importFrom GenomicRanges GRangesList findOverlaps
#' @importFrom stats quantile t.test binomial glm optim phyper
#' @importFrom poibin ppoibin
#' @param GRpeaks ATAC-seq or DNase-seq peaks in GRanges class where each row
#' represents the location of each peak. In addition there must be
#' 3 metadata columns called "reads" (representing read counts in each peak),
#' "GC" (representing the GC content), and lastly "peaktype" which designates
#' each peak as either ("target","background","no").
#' @param GRmotif Footprint calls from a footprint algorithm in GRanges
#' class where each row represents the location of each PWM occurrence.
#' The footprint calls in the forward strand
#'  and those in the backward strand from the same PWM are not differentiated.
#'  The row names of GRmotif are the motif IDs (e.g. MA01371 STAT1).
#' @return Returns a list of parameter alpha_k, p values from BiFET algorithm
#'  and p values from the hypergeometric test
#' @author Ahrim Youn
#' @note The function calculate_enrich_p first generates a TF binding
#' matrix M where i_th row represents i_th PWM and j_th column represetns j_th
#' peak with M_i, j=1 if the footprint of i_th PWM overlaps j_th peak and
#' 0 otherwise. Finally, the function “calculate_enrich_p” returns a list of
#' parameter alpha_k, enrichment p values from BiFET algorithm
#' and enrichment p values from the hypergeometric test.
#' @details In this example, the file input_peak_motif.Rdata was obtained
#' as follows: we used ATAC-seq data obtained from five human PBMC
#' (Ucar, et al., 2017) and five human islet samples (Khetan, et al., 2017)
#' and called peaks using MACS version 2.1.0 (Zhang, et al., 2008) with
#' parameters "-nomodel -f BAMPE". The peak sets from all samples were merged
#' to generate one consensus peak set (N = 57,108 peaks) by using
#' package DiffBind_2.2.5. (Ross-Innes, et al., 2012), where only the peaks
#' present at least in any two samples were included in the analysis.
#' We used the **summits** option to re-center each peak around the point of
#' greatest enrichment and obtained consensus peaks of same width (200bp).
#'
#' Out of these consensus peaks, we defined regions that are specifically
#' accessible in PBMC samples as regions where at least 4 PBMC samples have
#'  a peak, whereas none of the islet samples have a peak
#'  (n=4106 peaks; these regions are used as target regions in this example).
#' Similarly, we defined islet-specific peaks as those that were called as
#' a peak in at least 4 islet samples but none in any of the
#' PBMC samples (n=12886 peaks). The rest of the peaks excluding the
#' PBMC/islet-specific peaks were used as the background
#' (i.e., non-specific) peaks in our analyses (n=40116 peaks). For each peak,
#' GC content was obtained using peak annotation program
#' annotatePeaks.pl from the HOMER software (Heinz. et al., 2010).
#'
#' In this example, TF footprints were called using PIQ algorithm
#' (Sherwood, et al., 2014) using the pooled islet samples and pooled PBMC
#' samples to increase the detection power for TF footprints. We used only the
#' TF footprints that have a purity score greater than 0.9. The example file
#' contains footprint calls for only five PWMs from the
#' JASPAR database to reduce computing time.
#' @examples
#' # Load in the peak file and footprint calls from a footprint algorithm
#' peak_file <- system.file("extdata", "input_peak_motif.Rdata",
#'  package = "BiFET")
#' load(peak_file)
#' result <- calculate_enrich_p(GRpeaks,GRmotif)
#' @export
calculate_enrich_p <- function(GRpeaks, GRmotif) {
  
  # produce error message if input peak file contains peaks with GC content = 0
  if (any(GRpeaks$GC == 0)) {
    stop("Error!: Input peak object contains entries with GC content equal to zero. Please remove these entries prior to using this function")
  } else {}
  
  TFbinding.mat <- bindingTF_per_peak(GRpeaks, GRmotif)
  targetpeak <- which(GRpeaks$peaktype == "target")
  backgroundpeak <- which(GRpeaks$peaktype == "background")
  reads <- GRpeaks$reads
  GCcontent <- GRpeaks$GC
  backgroundpeak <- c(targetpeak, backgroundpeak)
  if(min(reads)==0) reads=reads+1
  upperlimit <- quantile(reads, 0.999)
  reads[reads > upperlimit] <- upperlimit
  q <- rep(NA, nrow(TFbinding.mat))
  alpha <- matrix(NA, nrow = nrow(TFbinding.mat), ncol = 2)
  p.mat <- matrix(NA, nrow(TFbinding.mat), ncol = 2)
  rownames(p.mat) <- rownames(TFbinding.mat)
  for (i in (seq_along(TFbinding.mat[, 1]))) {
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

bindingTF_per_peak <- function(GRpeaks, GRmotif) {
  overlaps <- as.matrix(findOverlaps(GRmotif, GRpeaks, minoverlap=2L))  
  allTF <- unique(names(GRmotif))
  TFbinding.mat <- matrix(0, nrow = length(allTF), ncol = length(GRpeaks))
  rownames(TFbinding.mat) <- allTF
  TFbinding.mat[cbind(match(names(GRmotif)[overlaps[,1]],allTF),overlaps[,2])] <- 1
  return(TFbinding.mat)
}

calculate_enrich_p_per_TF <- function(TFbinding.mat, reads,
                                      GCcontent, targetpeak, backgroundpeak) {
  GCbias.cutoff <- 0.05
  
  model <- glm(TFbinding.mat[, backgroundpeak]~reads[backgroundpeak] +
                 GCcontent[backgroundpeak],
               family = binomial(link = "logit"))
  GCcoef <- summary(model)$coef[3, ]
  if ( GCcoef[1] > 0 & GCcoef[4] < GCbias.cutoff) {
    res <- main_optim(rbind(TFbinding.mat[, backgroundpeak]),
                      reads[backgroundpeak], GCcontent[backgroundpeak], c(300, 0.05))
    q <- res[[1]]
    alpha <- res[[2]]
    p.mat <- enrich_p(reads, GCcontent, TFbinding.mat,
                      q, alpha, targetpeak, backgroundpeak)

  } else {
    res <- main_optim(rbind(TFbinding.mat[, backgroundpeak]),
                      reads[backgroundpeak],
                      rep(Inf, length(backgroundpeak)), c(300, 0.05))
    q <- res[[1]]
    alpha <- res[[2]]
    p.mat <- enrich_p(reads, rep(Inf, length(GCcontent)),
                      TFbinding.mat, q, alpha, targetpeak, backgroundpeak)
    alpha[2] <- 10 ^ (-8)
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

    temp <- optim(init, fn = loglik_per_TF, gr = deriv_per_TF,
                  lower = 10 ^ (-16),
                  upper = 1 / max(f(reads, alpha, GCcontent)) - 10 ^ (-10),
                  method = "L-BFGS-B", reads = reads, alpha = alpha,
                  GCcontent = GCcontent, TFbinding.mat = TFbinding.mat,
                  control = list(fnscale = -1))

    q <- temp$par

    temp <- optim(par = alpha[1], fn = total_lik1, gr = total_deriv1,
                  lower = 10 ^ (-8) + max(0, (reads / (log(q * f0(GCcontent,
                                 alpha[2]) + 2) - log(abs(q * f0(GCcontent, alpha[2]) -
                                     2))))[q * f0(GCcontent, alpha[2]) > 2]),
                  upper = Inf, method = "L-BFGS-B", alpha2 = alpha[2], qi = q,
                  reads = reads, GCcontent = GCcontent,
                  TFbinding.mat = TFbinding.mat, control = list(fnscale = -1))
    alpha[1] <- temp$par

    if (GCcontent[1] != Inf) {
      temp <- optim(par = alpha[2],
                    fn = total_lik2, gr = total_deriv2,
                    lower = 10 ^ (-8),                    
                    upper = 1.4, method = "L-BFGS-B", alpha1 = alpha[1],
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
  for (ii in seq_along(TFbinding.mat[, 1])) {
    no.footprint <- TFbinding.mat[ii, ]
    pp <- q[ii] * f(reads[targetpeak], alpha, GCcontent[targetpeak])
    pvalue.bifet[ii] <- 1 -
      ppoibin(kk = sum(no.footprint[targetpeak]) - 1, pp = pp)

    pick <- sum(no.footprint[backgroundpeak])
    whitepick <- sum(no.footprint[targetpeak]) - 1
    white <- length(targetpeak)
    black <- length(backgroundpeak) - length(targetpeak)

    pvalue.hyper[ii] <- 1 - phyper(whitepick, white, black, pick)
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
