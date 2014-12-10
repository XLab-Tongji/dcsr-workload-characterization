## MG_CLUSTER_TENDENCY
##
##
## Copyright (C) 2009  Distributed Computing System (DCS) Group, Computer
## Science Department - University of Piemonte Orientale, Alessandria (Italy).
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##
## Authors:
##   Marco Guazzone (marco.guazzone@gmail.com)
##

library("cluster");


## MG_CLUSTER.TENDENCY.HOPKINS
##
## Compute the Hopkin's statistic for the data matrix {\code x}.
##
## \param X Numeric matrix or data frame, where columns are variables and rows
##  are observations.
## \param metric Character string specifying the distance metrix to be used.
##  The currently availbale options are "euclidean" (the default).
## \param standardize Logical flag: if \c TRUE, then the measurements in \a X
##  are standardized before calculating the statistics.
##  Measurements are standardized for each variable (column), by subtracting the
##  variable's mean value and dividing by the variable's mean absolute
##  deviation.
## Calculates and compares the Euclidean distances between
## - randomly selected experimental objects and their nearest neighbours (W)
## - and artificially generated objects and their nearest experimental
##   neighbours (U). 
## .
## When U/(U+W) exceeds 0.75,  then the data set is considered to be
## significantly clustered.                              
## Forina's modification can be also applied with this algorithm.
##   H*= (Nreal*sum(U)-Npseudo*sum(W)) / (Nreal*sum(U)+Npseudo*sum(W))
## One important step in this modification is to select correctly 
## Npseudo, it is recommended to take a large value for it (here 80%). 
##
mg_cluster.tendency.hopkins <- function(X, win.size=0.1*nrow(X), metric=c("euclidean"), standardize=FALSE, n.iter=1, ...)
{
	# Gets the distance name and function
	metric <- match.arg(metric);
	if ( is.null(metric) )
	{
		stop( "[MG::CLUSTER::TENDENCY::HOPKINS] 'metric' parameter must be a distance measure name or a function." );
	}
	if (metric == "euclidean")
	{
		metric.fun <- function(x,y) { return(sqrt(sum((x-y)^2))); };
	}

	if (!standardize)
	{
		X <- scale(X, center=TRUE, scale=FALSE);
		sX <- colMeans(abs(X), na.rm=TRUE);
		if (0 %in% sX)
		{
			warning("'X' has constant columns; these are standardized to 0");
			sX[sX == 0] <- 1;
		}
		X <- scale(X, center=FALSE, scale=sX);
	}

	d <- dim(X);
	ranges <- apply(X, 2, range); # ==> 1st row: min; 2nd row: max
	# Generate a (win.size x d[2]) matrix of artifical points in space R^d[2]
	# (actually not properly in R^d[2] but in a smaller space that takes into
	# account the min and max ranges of the data
	Y <- apply(
			t(1:d[2]),
			2,
			function(ix) {
				# min and max are "enlarged" to allow runif to possibly
				# generate real extremes.
				cbind(runif(win.size, ranges[1,ix]-1e-07, ranges[2,ix]+1e-07));
			}
	);
	H.min <- +Inf;
	H.max <- -Inf;
	H.mean <- 0;
	for (i in 1:n.iter)
	{
#       S <- apply(
#               t(1:d[2]),
#               2,
#               function(ix) {
#                   sample(X[,ix], win.size, replace=FALSE);
#               }
#       );
		# Sample from original data
		ix <- sample(1:d[1], win.size, replace=FALSE);
		S <- X[ix,];

		# Compute distance between artificial and original points
		# and keep the minimum distance
		U <- apply(
			t(1:nrow(Y)),
			2,
			function(r1) {
				min(
					apply(
						t(1:nrow(X)),
						2,
						function(r2) {
							dist(
								rbind(Y[r1,],X[r2,])
							);
						}
					)
				);
			}
		);
		#W <- U <- rep(1,d[2]);
		#U <- X[S,S];
		#W <- X[T,S];
		# Compute distance between sample and original points
		# and keep the minimum distance
		W <- apply(
			t(1:nrow(S)),
			2,
			function(r1) {
				min(
					apply(
						t(1:nrow(X)),
						2,
						function(r2) {
							dist(
								rbind(S[r1,],X[r2,])
							);
						}
					)
				);
			}
		);
		Us <- sum(U^(d[2]));
		Ws <- sum(W^(d[2]));
		H <- Us/(Us+Ws);
		if (H < H.min)
		{
			H.min <- H;
		}
		if (H > H.max)
		{
			H.max <- H;
		}
		H.mean <- H.mean + H;
	}
	H.mean <- H.mean / n.iter;

	return(list(mean=H.mean, min=H.min, max=H.max));
}


## MG_CLUSTER.TENDENCY.VAT
##
## \param R A numeric matrix. If \a metric is NA or NULL, \a R is interpreted as
##  a dissimilarities (or distance) matrix; otherwise, \a R is a data matrix
##  and pair-wise distance are computed with respect to the distance metric
##  \a metric.
## \param metrix A character string representing the distance metric.
##  May be \c NA or \c NULL.
## \param stand Logical flag. If TRUE then the measurements in \a R are
##  standardized before computing the dissimilarities.
## \param plot Logical flag. A TRUE value causes an image to be plotted.
## \param ... Additional graphical parameters.
##
## References:
## - J. C. Bezdek and R. J. Hathaway.
##   "VAT: a tool for visual assessment of (cluster) tendency",
##   In Proc. of the 2002 International Joint Conference on Neural Networks,
##   2002.
## .
##
mg_cluster.tendency.vat <- function(R, metric="euclidean", stand=FALSE, plot=FALSE, ...)
{
	#R<-rbind(c(0,0.73,0.19,0.71,0.16),c(0.73, 0,0.59,0.12,0.78),c(0.19,0.59,0,0.55,0.19),c(0.71,0.12,0.55,0,0.74),c(0.16,0.78,0.19,0.74,0));
	#G.max <- 255;
	#G <- (0:G.max);
	#GG <- G/G.max;
	#image(t(R),col=grey(GG));
	#heatmap(R,col=grey(GG),Colv=NA,Rowv=NA,distfun=NA,symm=T)
	if (!(is.na(metric) || is.null(metric)))
	{
		# Compute pair-wise distances
		R <- as.matrix(daisy(R), metric=metric, stand=stand);
		# Transform distances to dissimilarities between 0 and 1
		R <- (R-min(R))/(max(R)-min(R));
		# Convert a distance object to a numeric matrix
		R <- as.matrix(R);
	}
	n <- nrow(R);
	R.max <- max(R);
	I <- J <- c();
	P <- vector(mode="integer", length=n);
	for (i in 1:n)
	{
		if (R.max %in% R[i,])
		{
			I <- c(I, i);
			P[1] <- i;
			break;
		}
	}
	J <- setdiff(1:n, I);
	for (r in 2:n)
	{
	#cat("r =>",r,"\n");
	#cat("I =>",paste(I,sep=","),"\n");
	#cat("J =>",paste(J,sep=","),"\n");
	#cat("P =>",paste(P,sep=","),"\n");
	#cat("R.min =>",R.min,"\n");
		R.min <- min(R[I,J]);
		for (j in 1:length(J))
		{
	#cat("j =>",j,"\n");
			if (R.min %in% R[I,j])
			{
	#cat("j =>",j,"Found MIN","\n");
				P[r] <- j;
				I <- c(I, j);
				J <- setdiff(J, j);
				break;
			}
		}
	}
	if (length(J) > 0)
	{
		P[J] <- J;
	}
	R.tilde <- R[P,P];

	if (plot)
	{
		image(t(R.tilde), ...);
	}

	return(invisible(R.tilde));
}

mg_cluster.tendency.revat <- function(X, plot=FALSE, ...)
{
	# Compute pair-wise distances
	R <- as.matrix(daisy(X));
	# Normalize distances to dissimilarities between 0 and 1
	R.max <- max(R);
	R.min <- min(R);
	R <- (R-R.min)/(R.max-R.min);
	n <- nrow(R);

	# Step.1
	k <- 0
	I <- 1:n;
	J <- c();
	S <- c();
	P <- matrix(NA,nrow=n,ncol=n);
	repeat
	{
		# Step.2
		i <- sample(I, 1);
		S <- union(S, i); # Note: instead of storing data we store positions
		# Step.3
		P[i,] <- 1-(R[i,]/R.max);
		#P[i,] <- 1-R[i,];
		p.mean <-mean(P[i,]);
		delta <-(1+p.mean)/2;
		# Step.4
		j <- which(P[i,] > delta); # the i-th "Tall Set"
		# Step.5
		if (length(j) > 0)
		{
			I <- setdiff(I, j);
			J <- union(J, j);
		}
		# Step.6
		if (length(I) == 0)
		{
			break;
		}
		k <- k+1;
	}


	P <- na.omit(P);

	if (plot)
	{
		rows <- nrow(P);
		par.old <- par(mfrow=n2mfrow(rows+1));
		sapply(1:rows, function(ix){ hist(P[ix,],xlab=paste(c("p[",ix,"]"),collapse=""),...); });
		image(t(R[S,S]), ...);
		par(mfrow=par.old);
	}

#	return(invisible(R.tilde));
	return(invisible(list(S=S,P=P)));
}
