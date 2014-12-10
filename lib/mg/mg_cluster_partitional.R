## MG_CLUSTER_PARTITIONAL
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

library('cluster');


## MG_CLUSTER.PARTITIONAL.KMEANS
##
mg_cluster.partitional.kmeans <- function(x, centers, iter.max=10, nstart=1, algorithm=c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
{
	res <- cluster::kmeans(
			x = X,
			centers = centers,
			iter.max = iter.max,
			nstart = nstart,
			algorithm = algorithm
	);

	return(res);
}


## MG_CLUSTER.PARTITIONAL.CLARA
##
## \param x data matrix or data frame, each row corresponds to an observation,
##  and each column corresponds to a variable.
##  All variables must be numeric.
##  Missing values (NAs) are allowed.
## \param k the number of clusters.
##  It is required that 0 < k < n where n is the number of observations
##  (i.e., n = nrow(x)).
## \param metric character string specifying the metric to be used for
##  calculating dissimilarities between observations.
##  The currently available options are "euclidean" and "manhattan".
##  Euclidean distances are root sum-of-squares of differences, and manhattan
##  distances are the sum of absolute differences.
## \param standardized logical flag, indicating if the measurements in x are
##  standardized before calculating the dissimilarities.
##  Measurements are standardized for each variable (column), by subtracting
##  the variable's mean value and dividing by the variable's mean absolute
##  deviation.
## \param samples.num  number of samples to be drawn from the dataset.
## \param sample.size number of observations in each sample.
##  Should be higher than the number of clusters (k) and at most the number of
##  observations (n = nrow(x)).
##
## RETURN:
##   An object of class "clara@ representing the clustering.
##
## SEE ALSO:
##   cluster
##
mg_cluster.partitional.clara <- function(x, k, metric="euclidean", standardized=FALSE, samples.num=5, sample.size=min(max(dim(x)),40+2*k))
{
	if (is.vector(x))
	{
		x <- as.matrix(x);
	} else {
		x <- data.matrix(x);
	}

	res <- cluster::clara(
		x,
		k,
		metric=metric,
		stand=standardized,
		samples=samples.num,
		sampsize=sample.size,
		rngR = TRUE # use the R's rng (making the results reproducible)
	);

	return(res);
}

