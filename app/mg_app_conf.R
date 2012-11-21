##
## Copyright (C) 2008       Marco Guazzone
##                          [Distributed Computing System (DCS) Group,
##                           Computer Science Institute,
##                           Department of Science and Technological Innovation,
##                           University of Piemonte Orientale,
##                           Alessandria (Italy)]
##
## This file is part of dcsr-workload-characterization
##
## dcsr-workload-characterization is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## dcsr-workload-characterization is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with dcsr-workload-characterization  If not, see <http://www.gnu.org/licenses/>.
##

## MG_APP_CONF
##
## SUMMARY
##   Global configurations.
##
## DESCRIPTION
##   Contains a set of global variables used as configuration for all "mg_app-*"
##   applications.
##
## AUTHOR
##   Marco Guazzone (marco.guazzone@gmail.com)
##

MG_APP_CONF_BASEPATH = "..";

# LCG
MG_APP_CONF_LCG_LOGPATH = paste( MG_APP_CONF_BASEPATH, "/traces/LCG", sep = "" );
MG_APP_CONF_LCG_LOGFULLNAME = paste( MG_APP_CONF_LCG_LOGPATH, "/LCG-2005-0", sep = "" );
MG_APP_CONF_LCG_VIEWPATH = paste( MG_APP_CONF_LCG_LOGPATH, "/views", sep = "" );
MG_APP_CONF_LCG_IMGPATH = paste( MG_APP_CONF_LCG_VIEWPATH, "/images", sep = "" );
MG_APP_CONF_LCG_STATSPATH = paste( MG_APP_CONF_LCG_VIEWPATH, "/stats", sep = "" );

# TeraGrid
MG_APP_CONF_TERAGRID_LOGPATH = paste( MG_APP_CONF_BASEPATH, "/traces/TeraGrid", sep = "" );
MG_APP_CONF_TERAGRID_LOGFULLNAME = paste( MG_APP_CONF_TERAGRID_LOGPATH, "/mout" , sep = "" );
MG_APP_CONF_TERAGRID_VIEWPATH = paste( MG_APP_CONF_TERAGRID_LOGPATH, "/views", sep = "" );
MG_APP_CONF_TERAGRID_IMGPATH = paste( MG_APP_CONF_TERAGRID_VIEWPATH, "/images", sep = "" );
MG_APP_CONF_TERAGRID_STATSPATH = paste( MG_APP_CONF_TERAGRID_VIEWPATH, "/stats", sep = "" );
