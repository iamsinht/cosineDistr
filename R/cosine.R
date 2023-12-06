#' cosine
#' 
#' This function computes the cosine similarity between all pairs of columns of two input matrices
#' @param mat1 Numeric, a matrix of size N x d1
#' @param mat2 Numeric, a matrix of size N x d2
#'
#' @return matrix
#'
#' @export
#' 

cosine <- function(mat1, mat2){
  a <- t(mat1) %*% mat2
  
  #mag1 <- sqrt(diag(t(mat1) %*% mat1))
  #mag2 <- sqrt(diag(t(mat2) %*% mat2))
  
  mag1 <- sqrt(colSums(mat1^2))
  mag2 <- sqrt(colSums(mat2^2))
  
  return(a/outer(mag1, mag2, "*"))
}



#' cosineApprox
#' 
#' This function computes the asymptotic approximation of cosine similarity used in the manuscript "On the 
#' distribution of cosine similarity." It is included here to simplify testing the properties of the
#' approximation. When computing similarities, always use cosine, not the asymptotic approximation. 
#' 
#' The calculation assumes that both matrices are samples from the same probability distribution over R^n
#' and that the column vectors are independently and identically distributed. In practice, when the population
#' parameters are not known, computing the sample mean and covariance eigenvalues from the data would be 
#' sensible. 
#' 
#' @param mat1 Numeric, a matrix of size N x d1.
#' @param mat2 Numeric, a matrix of size N x d2.
#' @param means Numeric, a vector of size N of population means in the basis with diagonal covariance.
#' @param lambda Numeric, the vector of eigenvalues of the covariance matrix. 
#'
#' @return matrix
#'
#' @export

cosineApprox <- function(mat1, mat2, means, lambda){
  
  a <- t(mat1) %*% mat2
  
  # lambda is the eigenvalues of the covariance matrix, i.e. lambda_i = sigma_i^2 in the basis with a 
  # diagonal covariance. 
  magSq <- sum(means**2 + lambda)
  
  return(a/magSq)
  
}