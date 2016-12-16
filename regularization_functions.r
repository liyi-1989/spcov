require(lorec)

#' regularization_by_thresholding
#' 
#' Regularizes a correlation matrix by thresholding
#' 
#' @param cor_matrix A correlation matrix
#' 
#' @param n The number of independent observations that were used to construct 
#' the 
#' 
#' @return The regularized correlation matrix
#' 
#' @references
#' Bickel, P. J. & Levina, E. Covariance regularization by thresholding. 
#' Ann. Statist. 36, 2577–2604 (2008). http://dx.doi.org/10.1214/08-AOS600 
#' 
regularization_by_thresholding <- function ( cor_matrix, MP=1 ) {
	# Regularization by thresholding
	# Bickel, P. J. & Levina, E. Covariance regularization by thresholding. 
	# Ann. Statist. 36, 2577–2604 (2008). http://dx.doi.org/10.1214/08-AOS600
	
	p <- ncol( cor_matrix )
	n=nrow(cor_matrix)
	regularized <- cor_matrix * ( abs(cor_matrix) > MP*sqrt(log(p)/n) )
	
	return( regularized )
}


#' regularization_by_convex_minimization
#' 
#' Regularizes a correlation matrix by convex minimization
#' 
#' @param cor_matrix A correlation matrix
#' 
#' @param n The number of independent observations that were used to construct 
#' the 
#' 
#' @return The regularized correlation matrix
#' 
#' @references
#' Luo, X. High Dimensional Low Rank and Sparse Covariance Matrix Estimation via 
#' Convex Minimization. arXiv.org (2011). http://arxiv.org/abs/1111.1133
#' 
regularization_by_convex_minimization <- function ( cor_matrix, n ) {

	p <- ncol( cor_matrix )
	
	re.lorec <- lorec( cor_matrix, diag(1, 100), diag(1, 100), sqrt(p/n), 
		sqrt(log(p)/n) )
	re.lorec.eig <- eigen( re.lorec$L )
	# threholding both
	sel <- re.lorec.eig$values > sqrt( p/n )
	if ( sum(sel) > 0 ) {
		V <- re.lorec.eig$vectors
		V <- V * ( abs(V)>sqrt(1/p) )
		V <- V[ ,1:sum(sel) ]
		Lhat <- V%*%diag( re.lorec.eig$values[1:sum(sel)] )%*%t(V)
	} else {
		Lhat <- matrix(0, p, p)
	}
	Shat <- re.lorec$S
	Shat <- Shat*(abs(Shat)>sqrt(1/p))
	regularized <- Lhat + Shat
	return( regularized )
}
