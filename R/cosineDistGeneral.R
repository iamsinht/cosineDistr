#' Sample the variance of cosine similarity for vectors drawn from a generalized multivariate normal
#' distribution with non-zero mean and with eigenvalues drawn from a Gamma distribution. 
#' 
#' @param ndim numeric, a vector of dimensions to sample. Default seq(100, 1000, 100)
#' @param iter Number of iterations to use for a particular dimension and eigenvalue distribution
#' @param distiter Number of eigenvalue distributions to run for each dimension
#' @param sPar Numeric, vector of gamma distribution hyperparameters for the eigenvalues. The eigenvalues
#' are sampled from a gamma distribution with shape parameter gamma(shape=sPar[1], rate=sPar[2]). and rate
#' parameter gamma(shape=sPar[3], rate=sPar[4])
#' @param meanSd Numeric, standard deviation of the normal distribution from which to sample means
#' @param nSamples Numeric, number of vectors to sample
#' @param distr Distribution to use - one of normal, laplace, mixture (of Gaussians), mixtureconst (of Gaussians with
#' the same covariance eigenvalue distribution)
#' 
#' @returns Dataframe with columns: dimension, theoretical variance based on eigenvalues, theoretical variance 
#' on observed eigenvalues, empirical variance, and root sum squared component means/sd. 
#' @export
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rmvl
#' @import stats MASS
cosineDistGeneral <- function(ndim=seq(100, 1000, 100), iter=3, distiter=10, sPar=c(2, 2, 2, 2), 
                            meanSd=2, nSamples=2000, distr="normal"){
  
  distr <- match.arg(distr, c("normal", "laplace", "mixture", "mixtureconst"))
  
  mylen <- length(ndim)*iter*distiter
  resdf <- data.frame(ndim=numeric(mylen), obsMean=numeric(mylen), theoryMean=numeric(mylen),
                      obsVar=numeric(mylen), theoryVar=numeric(mylen), 
                      tVarNoMean=numeric(mylen), rssZMeans=numeric(mylen), 
                      approxMean=numeric(mylen), approxVar=numeric(mylen),
                      meanSd=numeric(mylen), nSamples=numeric(mylen),
                      distribution=rep(distr, mylen))
  
  k <- 1
  mixcount <- 3
  
  for (mydim in ndim){
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
        
        xcos <- cosine(x, x)
        
        ovar <- stats::var(xcos[upper.tri(xcos)])
        tvar <- sum(lambda*(lambda + 2*(myMeans **2)))/(sum(myMeans**2 + lambda) ** 2)
        tvarNoMean <- sum(lambda ** 2)/(sum(lambda) ** 2)
        
        omean <- mean(xcos[upper.tri(xcos)])
        tmean <- sum(myMeans ** 2)/sum(myMeans**2 + lambda)
        
        # Run the approximation, useful for checking if aberrent behavior is due to the 
        # approximation or due to a miscalculation in the math. 
        
        approxcos <- cosineApprox(x, x, myMeans, lambda)
        
        approxvar <- stats::var(approxcos[upper.tri(approxcos)])
        approxmean <- mean(approxcos[upper.tri(approxcos)])
        
        rssZ <- sqrt(sum((myMeans ** 2)/lambda))
        
        resdf[k,] <- c(ndim=mydim, obsMean=omean, theoryMean=tmean, 
                       obsVar=ovar, theoryVar=tvar, tVarNoMean=tvarNoMean, 
                       rssZMeans=rssZ, approxMean=approxmean, approxVar=approxvar,
                       meanSd=meanSd, nSamples=nSamples, distribution=distr) 
        k <- k + 1
      }
    }
  }

  return (resdf)
  
}


#' Calculate the variance of cosine specifically for the two-dimensional case
#' 
#' @param myn - Number of vectors to sample at each ratio
#' @param iter - Number of different values of u to sample on [0, 0.5]
#' @param meanvec - Mean values to use (default (0,0)).
#' 
#' @returns dataframe with mean, variance, and approximation variance
#' @export
#' @import MASS stats
cosineDist2D <- function(myn=1000, iter=100, meanvec=c(0,0)){
  uvals <- seq(0.5/iter, 0.5, 0.5/iter)
  
  resdf <- data.frame(u=uvals, meanCos=numeric(length(uvals)), varCos=numeric(length(uvals)), theoryVar=numeric(length(uvals)))
  
  # S1, S2 are the eigenvalues of the covariance matrix, i.e. sigma_i^{2}
  for (myu in uvals){
    print(sprintf("iter = %d/%d", match(myu, uvals), iter))
    s1 <- myu
    s2 <- 1 - myu
    
    x <- t(MASS::mvrnorm(n=myn, mu=meanvec, Sigma=diag(c(s1, s2))))
    
    xcos <- cosine(x,x)
    
    meancos <- mean(xcos[upper.tri(xcos)])
    varcos <- stats::var(xcos[upper.tri(xcos)])
    theoryvar <- (s1^2 + s2^2)/((s1+s2)^2)
    
    resdf[match(myu, uvals), 2:4] <- c(meancos, varcos, theoryvar)
  }
  
  params <- c(myn=myn, iter=iter, meanvec=meanvec)
  
  return(list(resdf=resdf, params=params))
}
