## MG_CLUSTER_DENSITY
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

library("fpc");

mg_cluster.density.dbscan <- function(X, min.points=5, eps=NA, method=c("hybrid", "raw", "dist"))
{
	if (is.na(eps) || is.null(eps))
	{
		# Analytical way of estimating neighborhood radius for DBSCAN
		d <- dim(X);
		X.maxs <- sapply(1:d[2], function(ix){ max(X[,ix]); });
		X.mins <- sapply(1:d[2], function(ix){ min(X[,ix]); });
		eps <- (
			(
				prod(X.maxs - X.mins)*min.points*gamma(0.5*d[2]+1)
			) / (
				d[1]*sqrt(pi^d[2])
			)
		)^(1/d[2]);
	}

	res <- fpc::dbscan(X, eps, MinPts=min.points);

	return(res);
}

###DEPRECATED: use fpc::dbscan
#
#mg_cluster.density.distvector <- function(x,data){
#  ddata <- t(data)-x
#  dv <- apply(ddata^2,2,sum)
#}
#
## data may be nxp or distance matrix
## eps is the dbscan distance cutoff parameter
## MinPts is the minimum size of a cluster
## scale: Should the data be scaled?
## distances: has to be TRUE if data is a distance matrix
## showplot: Should the computation process be visualized?
## countmode: dbscan gives messages when processing point no. (countmode)
#mg_cluster.density.dbscan <- function(data,eps,MinPts=5, scale=FALSE, distances=FALSE, showplot=FALSE, countmode=c(1,2,3,5,10,100,1000,5000,10000,50000)){
#
#  data <- as.matrix(data)
#  n <- nrow(data)
#  if (scale) data <- scale(data)
#  unregpoints <- rep(0,n)
#
#  e2 <- eps^2
#  cv <- rep(0,n)
#  cn <- 0
#
#
#  i <- 1
#  for (i in 1:n){
#    if (i %in% countmode) cat("Processing point ", i," of ",n, ".\n")     unclass <- cv<1
#    if (cv[i]==0){
#
#      if (distances) seeds <- data[i,]<=eps
#      else{
#        seeds <- rep(FALSE,n)
#        seeds[unclass] <- distvector(data[i,],data[unclass,])<=e2
#      }
#      if (sum(seeds)+unregpoints[i]<MinPts) cv[i] <- (-1)
#      else{
#        cn <- cn+1
#        cv[i] <- cn
#        seeds[i] <- unclass[i] <- FALSE
#        unregpoints[seeds] <- unregpoints[seeds]+1
#        while (sum(seeds)>0){
#          if (showplot) plot(data,col=1+cv)
#          unclass[seeds] <- FALSE
#          cv[seeds] <- cn
#          ap <- (1:n)[seeds]
#
#
##          print(ap)
#
#
#          seeds <- rep(FALSE,n)          
#          for (j in ap){
#
#
##            if (showplot) plot(data,col=1+cv)
#
#
#            jseeds <- rep(FALSE,n)          
#            if (distances) jseeds[unclass] <- data[j,unclass]<=eps
#            else{
#              jseeds[unclass] <- distvector(data[j,],data[unclass,])<=e2
#            }
#            unregpoints[jseeds] <- unregpoints[jseeds]+1
#
#
##            if (cn==1)
#
##              cat(j," sum seeds=",sum(seeds)," unreg=",unregpoints[j],
#
#
##                  " newseeds=",sum(cv[jseeds]==0),"\n")
#
#
#            if (sum(jseeds)+unregpoints[j]>=MinPts){              
#              seeds[jseeds] <- cv[jseeds]==0
#              cv[jseeds & cv<0] <- cn
#            }
#          } # for j
#        } # while sum seeds>0
#      } # else (sum seeds + ... >= MinPts)
#
#    } # if cv==0
#  } # for i
#  if (sum(cv==(-1))>0){
#    noisenumber <- cn+1
#    cv[cv==(-1)] <- noisenumber
#  }
#  else
#    noisenumber <- FALSE
#  out <- list(classification=cv, noisenumber=noisenumber,
#
#              eps=eps, MinPts=MinPts, unregpoints=unregpoints)   out
#} # dbscan
## classification: classification vector
## noisenumber: number in the classification vector indicating noise points
## unregpoints: ignore... 
