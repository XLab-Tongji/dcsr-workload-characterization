## MG_CLUSTER_VALIDITY
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


## MG_CLUSTER_VALIDITY.SILHOUETTE
##
## \param x: a cluster object or an integer vector with k different integer
##  cluster codes or a list with such an "x$clustering" component.
##  Note that silhouette statistics are only defined if 2 <= k <= n-1.
##
mg_cluster.valid.silhouette <- function(x, ...)
{
	res <- cluster::silhouette(x, ...);

	return(res);
}
