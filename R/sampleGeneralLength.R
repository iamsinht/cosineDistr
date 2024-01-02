#' Calculate the empirical Euclidean length of multivariate normal random vectors with eigenvalues sampled
#' from a gamma distribution and non-zero mean sampled from a normal distribution
#'
#' @inheritParams cosineDistCentered
#' @param distr Distribution to use - one of normal, laplace, mixture (of Gaussians), mixtureconst (of Gaussians with
#' the same covariance eigenvalue distribution)
#' @param meanSd Numeric, standard deviation of the normal distribution from which to sample means
#' @param nSamples Numeric, number of vectors to sample
#'
#' @returns Dataframe with columns: ndim - dimension of vector, meanLen - average length of vectors, sdLen - 
#' standard deviation of length, theoryLength - theoretical length predicted according to Central Limit Theorem, Smith et al.
#' @export
#' @importFrom stats sd
getGenLength <- function(ndim=seq(100, 500, 100), iter=10, distiter=10, sPar=c(1,2,1,2),
                               distr="normal", meanSd = 2, nSamples=2000){
  
  distr <- match.arg(distr, c("normal", "laplace", "mixture", "mixtureconst"))
  
  mylen <- length(ndim)*iter*distiter
  resdf <- data.frame(ndim=numeric(mylen), meanLen=numeric(mylen), sdLen=numeric(mylen), 
                      theoryLen=numeric(mylen), tLenNoMean=numeric(mylen), 
                      distribution=rep(distr, mylen))
  k <- 1
  mixcount <- 3
  
  for (mydim in ndim){
    print(sprintf("mydim = %d", mydim))
    
    for (ii in seq_len(distiter)){
      print(sprintf("ii = %d", ii))
      a <- stats::rgamma(1, shape=sPar[1], rate=sPar[2])
      b <- stats::rgamma(1, shape=sPar[3], rate=sPar[4])
      
      for (jj in seq_len(iter)){
        sigmat <- diag(rgamma(mydim, shape=a, rate=b) + 1e-6)
        lambda <- diag(sigmat)
        
        myMeans <- stats::rnorm(n = mydim, mean = 0, sd = meanSd) * lambda
        
        if (distr == "normal"){
          x <- t(MASS::mvrnorm(n=nSamples, mu=myMeans, Sigma=sigmat))
        } else if (distr == "laplace"){
          x <- t(LaplacesDemon::rmvl(n=nSamples, mu=myMeans, Sigma=sigmat))
          
          # Test assumption of eigenvalues
          xpr <- prcomp(t(x), center=TRUE, scale=FALSE)
          lambda2 <- xpr$sdev^2
          myMeans2 <- xpr$center %*% xpr$rotation
          
        } else if (distr == "mixture"){
          
          a <- stats::rgamma(mixcount, shape=sPar[1], rate=sPar[2])
          b <- stats::rgamma(mixcount, shape=sPar[3], rate=sPar[4])
          x <- c()
          
          for (kk in seq_len(mixcount)){
            sigmat <- diag(rgamma(mydim, shape=a[kk], rate=b[kk]) + 1e-6)
            lambda <- diag(sigmat)
            myMeans <- stats::rnorm(n = mydim, mean = 0, sd = meanSd) * lambda          
            
            x <- cbind(x, t(MASS::mvrnorm(n=nSamples, mu=myMeans, Sigma=sigmat)))
          }
          
          # The parameters used for the theory calculations need to be calculated for the mixture
          # Calculating lambda, the eigenvalues of the covariance matrix, as follows only works 
          # because we've chosen a basis with diagonal covariance. Technically, each mixture is
          # constructed to have diagonal covariance, which leads to loss of generality.
          
          # Update: calculating marginal variances is completely wrong. Instead, extract the eigenvalues
          # of the covariance matrix. You cannot assume the covariance matrix of the mixture is diagonal.
          
          #myMeans <- rowMeans(x)
          #lambda <- apply(x, 1, var)
          xpr <- prcomp(t(x), center=TRUE, scale=FALSE)
          lambda <- xpr$sdev^2
          myMeans <- xpr$center %*% xpr$rotation
        } else if (distr == "mixtureconst"){
          a <- stats::rgamma(1, shape=sPar[1], rate=sPar[2])
          b <- stats::rgamma(1, shape=sPar[3], rate=sPar[4])
          x <- c()
          
          for (kk in seq_len(mixcount)){
            sigmat <- diag(rgamma(mydim, shape=a, rate=b) + 1e-6)
            lambda <- diag(sigmat)
            myMeans <- stats::rnorm(n = mydim, mean = 0, sd = meanSd) * lambda          
            
            x <- cbind(x, t(MASS::mvrnorm(n=nSamples, mu=myMeans, Sigma=sigmat)))
          }
          
          xpr <- prcomp(t(x), center=TRUE, scale=FALSE)
          lambda <- xpr$sdev^2
          myMeans <- xpr$center %*% xpr$rotation
        }
        
        meanLen <- mean(sqrt(colSums(x^2)))
        sdLen <- stats::sd(sqrt(colSums(x^2)))
        
        theoryLen <- sqrt(sum(myMeans^2 + lambda))
        tLenNoMean <- sqrt(sum(lambda))
        
        resdf[k,] <- data.frame(ndim=mydim, 
                                meanLen=meanLen, 
                                sdLen=sdLen, 
                                theoryLen=theoryLen, 
                                tLenNoMean=tLenNoMean, 
                                distribution=distr)
        k <- k + 1
      }
    }
  }
  
  return (resdf)
  
}