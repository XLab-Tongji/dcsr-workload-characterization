## MG_CLUSTER_MBC
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

library("mclust");

mg_cluster.mbc.gaussian <- function(x, k=2)
{
	res <- mclust::Mclust(x, G=k);

	return(res);
}
