# Compute the distribution of cosine similarity for vectors of normally distributed random variables
# to illustrate how the distribution of cosine distance changes with increasing dimension

#' Simulate cosine distribution for spherically symmetric data 
#' 
#' `cosineDistSpherical` samples from spherically symmetric data (0 mean, identity covariance)
#' for varying dimension and computes the mean and variance of cosine similarity of standard 
#' iid normal vectors of specified dimension.
#' 
#' @param ndim numeric vector, length values to evaluate. Default: seq(5, 1000, 5)
#' @param N number of vectors to sample
#' 
#' @returns A dataframe with the dimension, mean, standard deviation, and variance of the cosine 
#' similarity of the sample, excluding the diagonal
#' 
#' @export
#' @import ggplot2 stats MASS
cosineDistSpherical <- function(ndim=seq(5, 1000, 5), N=1000){
  
  retdf <- data.frame(dim=numeric(), mean=numeric(), sd=numeric(), var=numeric())
  
  for (ii in ndim){
    print(ii)
    
    x <- t(MASS::mvrnorm(n=N, mu=numeric(ii), Sigma=diag(ii)))
    xcos <- cosine(x, x)
    
    retdf <- rbind(retdf, data.frame(dim=ii, mean=mean(xcos[upper.tri(xcos)]), 
                                     sd=stats::sd(xcos[upper.tri(xcos)]), 
                                     var=stats::var(xcos[upper.tri(xcos)])))
  }
  
  return(retdf)
}


#' Get cosine similarity distribution - for IID normals
#' 
#' @param mydim numeric, dimension of the vector space
#' @param N numeric, number of elements to sample from the vector space
#' @param mkfig logical, whether to create a figure. 
#' 
#' @returns the cosine similarity matrix of dimension N x N
#' @export
#' @import ggplot2 MASS
getCosineDist <- function(mydim=500, N=1000, mkfig=1){
  x <- t(MASS::mvrnorm(n=N, mu=numeric(mydim), Sigma=diag(mydim)))
  xcos <- cosine(x,x)
  
  ndist <- data.frame(x=seq(-1,1,0.002), y=0.5*stats::dbeta(seq(0,1,0.001), shape1=(mydim-1)/2, shape2=(mydim-1)/2))
  
  if (mkfig){
    xdf <- data.frame(val=xcos[upper.tri(xcos)])
    cols <- c("empirical"=ggplot2::alpha("black", 0.8), "theoretical"="orange")
    print(ggplot2::ggplot() + ggplot2::geom_density(data=xdf, aes(x=val, color="empirical", linetype="empirical"), size=2) + 
            geom_line(data=ndist, aes(x=x, y=y, color="theoretical", linetype="theoretical"), size=1) + xlab("Cosine similarity") + ylab("Density") + 
            ggtitle(sprintf("Cosine distribution of dimension %d, N=%d", mydim, N)) + theme_minimal() + 
            scale_color_manual("Distribution", breaks=names(cols), values=cols) + 
            scale_linetype_manual("Distribution", values=c("empirical"=1, "theoretical"=5)))  # +theme(legend.position="right"))
  }
  return(xcos)
}
