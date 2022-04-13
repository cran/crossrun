#' Joint Distribution for Crossings and Runs Using the Empirical Median.
#'
#' @description Joint probability distribution for the number of crossings
#' C and the longest run L in a sequence of n Bernoulli observations where
#' the number of successes is fixed at m, m between 0 and n. For fixed n,
#' the joint distribution is computed for all m, this makes the computation
#' demanding in terms of time and storage requirements. The joint distribution
#' is computed separately for sequences where the first observation is, or is
#' not, a success. The results are mainly intended for use when n is even and
#' m=n/2, but computation in this case requires that all distributions are
#' computed previously for all m, for all shorter sequences (lower n). In the
#' case of even n and m=n/2, the distributions for sequences starting or not
#' with a success are identical, and only the distribution among sequences
#' starting with a success is used. In that case, this may be interpreted as
#' the joint distribution for sequences around the empirical median.
#' @param nmax max sequence length.
#' @param prec mpft precision.
#' @param printn logical for progress output.
#' @return nfi, number of sequences with m successes, starting with a success, and
#' nfn, number of sequences with m successes, not starting with a success.
#' Three-dimensional Rmpfr arrays for each n up to nmax, with dimensions
#' n (C=0 to n-1), n (L=1 to n) and n+1 (m=0 to n). For n even and m=n/2,
#' only nfi, and the part corresponding to C=1 to n-1 and L=1 and m=n/2
#' is non-zero and should be used.
#' @examples
#' crem14 <- crossrunem(nmax=14, printn=TRUE)
#' Rmpfr::asNumeric(crem14$nfi[[14]][,,"m=7"]) # subsets of size 7=14/2
#' # restricted to possible values of C and L
#' Rmpfr::asNumeric(crem14$nfi[[14]][-1,1:7,"m=7"]) # same as stored data joint14em
#' Rmpfr::asNumeric(crem14$nfn[[14]][-1,1:7,"m=7"]) # the same
#'
#' # subsets of sizes different from 14/2
#' # size 4, first observation included
#' Rmpfr::asNumeric(crem14$nfi[[14]][,,"m=4"])
#' # size 14-4=10, first observation not included
#' Rmpfr::asNumeric(crem14$nfn[[14]][,,"m=10"]) # the same
#'
#' @export
crossrunem <- function(nmax = 100, prec = 120,
                       printn = FALSE) {
  # conditioning of S= first value included in the subset (S=1)
  # or not (S=0).
  # nfi: number of subsets, first value included,
  # nfn: number of subsets, first value not included.
  # separate code by brute force for low n (n=1,2,3,4):
  nfi <- list(Rmpfr::mpfrArray(0, prec, dim = c(1, 1, 2)))
  nfn <- list(Rmpfr::mpfrArray(0, prec, dim = c(1, 1, 2)))
  dimnames(nfi[[1]]) <- list("c=0","l=1",c("m=0","m=1"))
  dimnames(nfn[[1]]) <- list("c=0","l=1",c("m=0","m=1"))
  nfn[[1]][1,1,1] <- 1 # m=0
  nfi[[1]][1,1,2] <- 1 # m=1
  nfi[[2]] <- Rmpfr::mpfrArray(0, prec, dim = c(2, 2, 3))
  nfn[[2]] <- Rmpfr::mpfrArray(0, prec, dim = c(2, 2, 3))
  dimnames(nfi[[2]]) <- list(paste0("c=",0:1),paste0("l=",1:2),
                             paste0("m=",0:2))
  dimnames(nfn[[2]]) <- list(paste0("c=",0:1),paste0("l=",1:2),
                             paste0("m=",0:2))
  nfn[[2]][1,2,1] <- 1 # m=0
  nfi[[2]][2,1,2] <- 1 # m=1
  nfn[[2]][2,1,2] <- 1
  nfi[[2]][1,2,3] <- 1 # m=2
  nfi[[3]] <- Rmpfr::mpfrArray(0, prec, dim = c(3, 3, 4))
  nfn[[3]] <- Rmpfr::mpfrArray(0, prec, dim = c(3, 3, 4))
  dimnames(nfi[[3]]) <- list(paste0("c=",0:2),paste0("l=",1:3),
                             paste0("m=",0:3))
  dimnames(nfn[[3]]) <- list(paste0("c=",0:2),paste0("l=",1:3),
                             paste0("m=",0:3))
  nfn[[3]][1,3,1] <- 1 # m=0
  nfi[[3]][2,2,2] <- 1 # m=1
  nfn[[3]][2,2,2] <- 1
  nfn[[3]][3,1,2] <- 1
  nfi[[3]][2,2,3] <- 1 # m=2
  nfi[[3]][3,1,3] <- 1
  nfn[[3]][2,2,3] <- 1
  nfi[[3]][1,3,4] <- 1 # m=3
  nfi[[4]] <- Rmpfr::mpfrArray(0, prec, dim = c(4, 4, 5))
  nfn[[4]] <- Rmpfr::mpfrArray(0, prec, dim = c(4, 4, 5))
  dimnames(nfi[[4]]) <- list(paste0("c=",0:3),paste0("l=",1:4),
                             paste0("m=",0:4))
  dimnames(nfn[[4]]) <- list(paste0("c=",0:3),paste0("l=",1:4),
                             paste0("m=",0:4))
  nfn[[4]][1,4,1] <- 1 # m=0
  nfi[[4]][2,3,2] <- 1 # m=1
  nfn[[4]][2,3,2] <- 1
  nfn[[4]][3,2,2] <- 2
  nfi[[4]][2:3,2,3] <- 1 # m=2
  nfi[[4]][4,1,3] <- 1
  nfn[[4]][2:3,2,3] <- 1
  nfn[[4]][4,1,3] <- 1
  nfi[[4]][2,3,4] <- 1 # m=3
  nfi[[4]][3,2,4] <- 2
  nfn[[4]][2,3,4] <- 1
  nfi[[4]][1,4,5] <- 1 # m=4
  # iterative procedure for higher n (>=5):
  if (nmax>4) for (nn in 5:nmax) {
    nfi[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    nfn[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    dimnames(nfi[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    dimnames(nfn[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    # separate computation for m=0,n:
    nfn[[nn]][1,nn,1] <- 1 # m=0
    nfi[[nn]][1,nn,nn+1] <- 1 # m=n
    # iterative procedure nfi, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:mm) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + nfn[[1]][1,1,1]
        else
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),,mm-gg+1])[,nn-gg]
      }
      if (gg<nn-gg) {
        if (gg==1) {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1]
        }
        else {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1])[,gg]
        }
        nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfn[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm-gg+1]
      } # end low g
    } # end iterative procedure nfi
    # iterative procedure nfn, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:(nn-mm)) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + nfi[[1]][1,1,2]
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),,mm+1])[,nn-gg]
        }
      } # end high g
      if (gg<nn-gg) {
        if (gg==1) {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1]
        }
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1])[,gg]
        }
        nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfi[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm+1]
      } # end low g
    } # end iterative procedure nfi
    if (printn==TRUE) print(nn)
    if (printn==TRUE) print(Sys.time())
  } # end for nn
  return(list(nfi = nfi, nfn = nfn))
} # end function crossrunem

#' Continuation of an existing sequence of joint probabilities for
#' crossings and longest run, based on the empirical median.
#'
#' @description  Continuation of an existing sequence of the number of
#' crossings C and the longest run L in a sequence of n independent
#' continuous observations classified as above or below the empirical
#' median. To enhance precision, results are stored in mpfr arrays and
#' the probabilities are multiplied by \eqn{choose(n,m)/2} where m=n/2,
#' even n assumed. The probabilities are integers in this representation.
#'
#' @param emstart existing sequence
#' @param n1 sequence length for the first new case addedc
#' @param nmax max sequence length.
#' @param prec mpft precision.
#' @param printn logical for including progress output.
#' @return nfi, number of sequences with m successes, starting with a success, and
#' nfn, number of sequences with m successes, not starting with a success.
crossrunemcont <- function(emstart, n1=61, nmax = 100,
                           prec = 120, printn = FALSE) {
  nfi <- emstart$nfi
  nfn <- emstart$nfn
  for (nn in n1:nmax) {
    nfi[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    nfn[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    dimnames(nfi[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    dimnames(nfn[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    # separate computation for m=0,n:
    nfn[[nn]][1,nn,1] <- 1 # m=0
    nfi[[nn]][1,nn,nn+1] <- 1 # m=n
    # iterative procedure nfi, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:mm) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + nfn[[1]][1,1,1]
        else
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),,mm-gg+1])[,nn-gg]
      }
      if (gg<nn-gg) {
        if (gg==1) {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1]
        }
        else {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1])[,gg]
        }
        nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfn[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm-gg+1]
      } # end low g
    } # end iterative procedure nfi
    # iterative procedure nfn, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:(nn-mm)) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + nfi[[1]][1,1,2]
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),,mm+1])[,nn-gg]
        }
      } # end high g
      if (gg<nn-gg) {
        if (gg==1) {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1]
        }
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1])[,gg]
        }
        nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfi[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm+1]
      } # end low g
    } # end iterative procedure nfi
    if (printn==TRUE) print(nn)
    if (printn==TRUE) print(Sys.time())
  } # end for nn
  return(list(nfi = nfi, nfn = nfn))
} # end function crossrunemcont

#' Check of joint probabilities by simulations
#'
#' @description Simulation of a sequence of n=2m observations
#' around the median in the sequence. To be used for checking
#' the results of crossrunem.
#' @param m1, half the sequence length
#' @param nsim number of simulations
#' @return data frame with cs, number of crossings
#' and ls, longest run in the simulations.
#' @examples
#' simclem14 <- simclem(nsim=sum(joint14em))
#' print(table(simclem14)) # joint distributions in the simulations
#' print(joint14em) # for comparison
#'
#' @export
simclem <- function(m1=7, nsim = 100000) {
  n1 <- 2*m1
  cs <- rep(NA, nsim)
  ls <- rep(NA, nsim)
  for (sim in 1:nsim) {
    series <- stats::rnorm(n1)
    med <- stats::median(series)
    above <- as.numeric(series>med)
    rleabove <- rle(above)$lengths
    cs[sim] <- length(rleabove) - 1
    ls[sim] <- max(rleabove)
  } # end for sim
  return(data.frame(cs=cs,ls=ls))
} #end function simclem
