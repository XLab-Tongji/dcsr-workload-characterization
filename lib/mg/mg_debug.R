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

##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## Class representation
setClass("file")
setClass("connection")
setClass("terminal")
setIs("file", "connection")
setIs("terminal", "connection")

setClass(
        "mg_debug",
        representation( active="logical", con="connection" )
);

## Constructor
setMethod(
        "initialize",
        "mg_debug",
        function( .Object, activate = TRUE, con = stderr() )
        {
                .Object <- callNextMethod()
                .Object@active <- activate | !is.null(con);
		.Object@con <- con;
                return( .Object )
        }
);

####
## Generic methods
####

setGeneric(
	"mg_debug.message",
	function( o, ..., domain = NULL, appendLF = TRUE ) { standardGeneric( "mg_debug.message" ) }
);
setMethod(
	"mg_debug.message",
	"mg_debug",
	function( o, ..., domain = NULL, appendLF = TRUE )
	{
		if ( !o@active )
		{
			return( invisible() );
		}

		args <- list(...);

		cond <- if (length(args) == 1 && inherits(args[[1]], "condition")) {
			if (nargs() > 1) 
			{
				warning("[MG::DEBUG::MESSAGE] Additional arguments ignored in message()");
			}
			args[[1]]
		} else {
			msg <- .makeMessage(..., domain = domain, appendLF = appendLF)
			call <- sys.call()
			simpleMessage(msg, call)
		};

		defaultHandler <- function(c) {
			cat(conditionMessage(c), file = o@con, sep = "")
		};

		withRestarts({
			signalCondition(cond)
			defaultHandler(cond)
		}, muffleMessage = function() NULL);

		return( invisible() );
	}
);

setGeneric(
	"mg_debug.active",
	function( o, ... ) { standardGeneric( "mg_debug.active" ) }
);
setMethod(
	"mg_debug.active",
	"mg_debug",
	function( o, ... )
	{
		dots <- list(...);

		if ( length( dots ) == 0 )
		{
			return( o@active );
		}

		return( invisible( o@active <- as.logical( dots[[1]] ) ) );
	}
);
