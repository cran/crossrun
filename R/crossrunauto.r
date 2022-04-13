#' Joint Distribution for Crossings and Runs, autocorrelated Sequence
#'
#' @description Joint probability distribution for the number of crossings
#' C and the longest run L in a sequence of n  autocorrelated Bernoulli
#' observations with success probability p. To enhance precision, results
#' are stored in mpfr arrays and the probabilities are multiplied by
#' \eqn{m^{n-1}} for a multiplier m.
#' @param nmax max sequence length.
#' @param prob success probability p.
#' @param changeprob unrestricted change probability. If \eqn{p \geq 0.5},
#'  probability of changing to success, if not probability of changing to
#'  failure.
#' @param mult multiplier for joint probabilities.
#' @param prec mpft precision.
#' @param printn logical for progress output.
#' @return list of joint probabilities.
#' @examples
#' # p=0.6, independence
#' cr10.6 <- crossrunbin(nmax=10, prob=0.6, printn=TRUE)
#' cra10.6 <- crossrunauto(nmax=10, prob=0.6, changeprob=.6, printn=TRUE)
#' Rmpfr::asNumeric(cr10.6$pt[[10]])
#' Rmpfr::asNumeric(cr10.6$pt[[10]])
#' Rmpfr::asNumeric(cr10.6$pt[[10]]) - Rmpfr::asNumeric(cra10.6$pt[[10]]) # equal
#'
#'
#' # p=0.6, some dependence
#' cr10.6 <- crossrunbin(nmax=10, prob=0.6, printn=TRUE)
#' cra10.6.u.5 <- crossrunauto(nmax=10, prob=0.6, changeprob=.5, printn=TRUE)
#' round(Rmpfr::asNumeric(cr10.6$pt[[10]]),1)
#' round(Rmpfr::asNumeric(cra10.6.u.5$pt[[10]]),1) # not the same
#' @export
crossrunauto <- function(nmax = 100, prob = 0.5, changeprob = 0.5,
                         mult = 2, prec = 120, printn = FALSE) {
    nill <- Rmpfr::mpfr(0, prec)
    one <- Rmpfr::mpfr(1, prec)
    half <- Rmpfr::mpfr(0.5, prec)
    multm <- Rmpfr::mpfr(mult, prec)
    pm <- Rmpfr::mpfr(prob, prec)
    qm <- one - pm
    changeprobm <- Rmpfr::mpfr(changeprob, prec)
    if (pm>=0.5) {
        upprobm <- changeprobm
        downprobm <- (one-pm)*upprobm/pm
    } else {
        downprobm <- changeprobm
        upprob <- pm*upprobm/(one-pm)
    }
    corrm <- one-upprobm/pm
    pmultm <- pm * multm
    qmultm <- qm * multm
    umultm <- upprobm * multm # multiplied probability for going up
    dmultm <- downprobm * multm  # multiplied probability for going down
    topmultm <- (one-downprobm) * mult  # multiplied probability for staying at top
    botmultm <- (one-upprobm) * mult  # multiplied probability for staying at bottom
    # conditioning of S= first value, pat: above 0, pbt: below 0 suffix
    # t: probabilities times multm^(n-1).  n=1:
    pat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    pbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    pt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    for (nn in 2:nmax) {
        pat[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
        pbt[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
        rownames(pat[[nn]]) <- c(0:(nn - 1))
        rownames(pbt[[nn]]) <- c(0:(nn - 1))
        colnames(pat[[nn]]) <- c(1:nn)
        colnames(pbt[[nn]]) <- c(1:nn)
        pat[[nn]][1, nn] <- (topmultm^(nn - 1))  # from cond on no crossing
        pbt[[nn]][1, nn] <- (botmultm^(nn - 1))  # from cond on no crossing
        for (ff in 2:nn) {
            # from cond on first crossing at ff if last part shortest:
            if (nn - ff + 1 <= ff - 1)
            {
                f1 <- ff  # unnecessary, but makes code checking easier
                pat[[nn]][2:(nn-f1+2), f1-1] <- pat[[nn]][2:(nn-f1+2),f1-1] +
                    (topmultm^(f1-2))*dmultm*qbt[[nn-f1+1]][1:(nn-f1+1), nn-f1+1]
                pbt[[nn]][2:(nn-f1+2), f1-1] <- pbt[[nn]][2:(nn-f1+2), f1-1] +
                    (botmultm^(f1-2))*umultm*qat[[nn-f1+1]][1:(nn-f1+1), nn-f1+1]
            }  # end if last part shortest
            if (nn - ff + 1 > ff - 1)
            {
                # if last part longest
                f2 <- ff  # unnecessary, but makes code checking easier
                pat[[nn]][2:(nn-f2+2), f2-1] <- pat[[nn]][2:(nn-f2+2), f2-1] +
                    (topmultm^(f2-2))*dmultm*qbt[[nn-f2+1]][1:(nn-f2+1), f2-1]
                pat[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] <- pat[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] +
                    (topmultm^(f2-2))*dmultm*pbt[[nn-f2+1]][1:(nn-f2+1), f2:(nn-f2+1)]
                pbt[[nn]][2:(nn-f2+2), f2-1] <- pbt[[nn]][2:(nn-f2+2), f2-1] +
                    (botmultm^(f2-2))*umultm*qat[[nn-f2+1]][1:(nn-f2+1), f2-1]
                pbt[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] <- pbt[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] +
                    (botmultm^(f2-2))*umultm*pat[[nn-f2+1]][1:(nn-f2+1), f2:(nn-f2+1)]
            }  # end if last part longest
        }  # end for ff
        pt[[nn]] <- pm * pat[[nn]] + qm * pbt[[nn]]
        qat[[nn]] <- cumsumm(pat[[nn]])
        qbt[[nn]] <- cumsumm(pbt[[nn]])
        qt[[nn]] <- pm * qat[[nn]] + qm * qbt[[nn]]
        rownames(pt[[nn]]) <- c(0:(nn - 1))
        colnames(pt[[nn]]) <- c(1:nn)
        rownames(qat[[nn]]) <- c(0:(nn - 1))
        colnames(qat[[nn]]) <- c(1:nn)
        rownames(qbt[[nn]]) <- c(0:(nn - 1))
        rownames(qat[[nn]]) <- c(0:(nn - 1))
        colnames(qt[[nn]]) <- c(1:nn)
        colnames(qt[[nn]]) <- c(1:nn)
        if (printn)
        {
            print(nn)
            print(Sys.time())
        }  # end optional timing information
    }  # end for nn
    names(pat) <- paste("pat", 1:nmax, sep = "")
    names(pbt) <- paste("pbt", 1:nmax, sep = "")
    names(pt) <- paste("pt", 1:nmax, sep = "")
    names(qat) <- paste("qat", 1:nmax, sep = "")
    names(qbt) <- paste("qbt", 1:nmax, sep = "")
    names(qt) <- paste("qt", 1:nmax, sep = "")
    return(list(pat = pat, pbt = pbt, pt = pt, qat = qat, qbt = qbt,
                qt = qt))
} # end function crossrunauto
