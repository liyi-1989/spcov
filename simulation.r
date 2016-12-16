source("regularization_functions.r")
require(picante)
require(ape) 
require(geiger)
require(fields)

set.seed(123456)


# A couple functions to facilitate simulation

test_cov_matrix <- function() {
	# Builds a simple character covariance matrix for simulating data
	
	G <- 100
	trueCovariance = matrix(0,G,G)
	trueCovariance[1:10,1:10] = 0.95
	trueCovariance[11:80,11:80] = 0.3
	trueCovariance[81:100,81:100] = 0.7
	diag(trueCovariance) = 1

	return( trueCovariance ) 
}

plot_matrix <- function(m, ... ) {
	# Basic plotting
	nr <- nrow(m)
	nc <- ncol(m)
	image(1:nc, 1:nr, t(m[nr:1, ]), axes=F,xlab="", ylab="", ... )
}


# Build the true matrix
trueCovariance = test_cov_matrix()

# Run the phylogenetic analyses
tree_text <- "(((Species_A:10.0,Species_B:10.0)I:40.0,(Species_C:40.0,Species_D:40.0)J:10.0)K:10.0,(Species_E:45.0,(Species_F:30.0,(Species_G:15.0,Species_H:15.0)O:15.0)N:15.0)M:15.0)L:0.0;"
phy <- read.tree( text=tree_text )
Z <- sim.char(phy, trueCovariance, nsim = 1, model = "BM", root = 1)
W <- Z[,,1]	# Grab one simulation
contrasts <- apply( W, 2, function(a) pic(a, phy) )


# The covariance matrix estimated directly from the contrasts

contrastcor <- cor.table(contrasts)$r

n <- nrow( contrasts )

# Regularization by thresholding

bickel <- regularization_by_thresholding( contrastcor, n )

# Regularization by convex Minimization

luo <- regularization_by_convex_minimization( contrastcor, n )


# Set up plotting parameters
# Color plots
# mycols <- tim.colors(100)
# Greyscale plots
mycols <- gray(1:100 / 100)

mybreaks <- seq(-1, 1, length.out=101)


# Plot the results

pdf(file="regularization.pdf", width=7, height=3.5)
layout(rbind(c(1,2,3,4),c(0,5,5,0)), heights=c(.75,.5))

m <- 1
par( mar=c(m,m,m,m) )
par( oma=c(0,0,0,0) )

plot_matrix( trueCovariance, main="(a) True", breaks=mybreaks, col=mycols )
plot_matrix( contrastcor, main="(b) Contrasts", breaks=mybreaks, col=mycols )
plot_matrix( luo, main="(c) Minimization", breaks=mybreaks, col=mycols )
plot_matrix( bickel, main="(d) Thresholding", breaks=mybreaks, col=mycols )

frame()
par( mar=c(0,0,0,0))
image.plot( legend.only=TRUE, zlim=range(-1:1), col=mycols, legend.width=3, 
	horizontal=TRUE)

dev.off()


# Write sessionInfo() results to file to record version numbers.

si <- toLatex( sessionInfo() )

sink( "sessionInfo.txt" )
cat( si )
sink( )
