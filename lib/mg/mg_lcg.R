## vim: set tabstop=4 expandtab shiftwidth=4 softtabstop=4:

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

## MG_LCG
##
## SUMMARY
##  LCG class.
##
## DESCRIPTION
##  Class for LCG log.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_plot.R" );
source( "lib/mg/mg_fit.R" );

## Class representation
mg_lcg <- setClass(
	"mg_lcg",
	slots = c( data="data.frame" )
);

## Constructor
setMethod(
    "initialize",
    "mg_lcg",
    function( .Object, file )
    {
        #.Object <- callNextMethod()
        .Object@data <- read.table( file, header = FALSE )
        # timestamp: Unix timestamp of submit time
        # vo: group name
        # uid: user ID
        # node: compute element name
        # runtime: job runtime
        names(.Object@data) <- c("timestamp", "vo", "uid", "node", "runtime");
        return( .Object )
    }
);

####
## Geneic methods
####

setGeneric(
	"mg_lcg_voNames",
	function( o, sorted = FALSE ) { standardGeneric( "mg_lcg_voNames" ) }
);
setMethod(
	"mg_lcg_voNames",
	"mg_lcg",
	function( o, sorted = FALSE )
	{
		vos <- levels( o@data$vo );
		if ( sorted )
		{
			vos - sort( vos );
		}
		return( vos )
	}
);

setGeneric(
	"mg_lcg_voFrequency",
	function( o, whatVO ) { standardGeneric( "mg_lcg_voFrequency" ) }
);
setMethod(
	"mg_lcg_voFrequency",
	"mg_lcg",
	function( o, whatVO )
	{
		return( length( o@data$vo[ o@data$vo == whatVO ] ) );
	}
);

####
## Number of Jobs
####

## Overall number of jobs by date (no time)
setGeneric(
	"mg_lcg_overallNJobsByDate",
	function( o ) { standardGeneric( "mg_lcg_overallNJobsByDate" ) }
);
setMethod(
	"mg_lcg_overallNJobsByDate",
	"mg_lcg",
	function( o )
	{
		jan1970 <- ISOdatetime( 1970, 1, 1, 0, 0, 0,"GMT" );

		x <- aggregate(
			o@data$timestamp,
			by = list( as.Date( jan1970 + o@data$timestamp ) ),
			length
		)

		return( list( date = levels( x$Group ), njobs = x$x ) )
	}
);

## Number of jobs per VO.
setGeneric(
	"mg_lcg_njobsByVO",
	function( o ) { standardGeneric( "mg_lcg_njobsByVO" ) }
);
setMethod(
	"mg_lcg_njobsByVO",
	"mg_lcg",
	function( o )
	{
		x <- aggregate(
			o@data$vo,
			by = list( o@data$vo ),
			length
		)

		return( list( vo = levels( x$Group ), njobs = x$x ) )
	}
);

## Number of jobs per User
setGeneric(
	"mg_lcg_njobsByUser",
	function( o ) { standardGeneric( "mg_lcg_njobsByUser" ) }
);
setMethod(
	"mg_lcg_njobsByUser",
	"mg_lcg",
	function( o )
	{
		x <- aggregate(
			o@data$uid,
			by = list( o@data$uid ),
			length
		)

		return( list( uid = levels( x$Group ), njobs = x$x ) )
	}
);

# VO Number of jobs by Date
setGeneric(
	"mg_lcg_voNJobsByDate",
	function( o, whatVO ) { standardGeneric( "mg_lcg_voNJobsByDate" ) }
);
setMethod(
	"mg_lcg_voNJobsByDate",
	"mg_lcg",
	function( o, whatVO )
	{
		jan1970 <- ISOdatetime( 1970, 1, 1, 0, 0, 0,"GMT" );

		t <- o@data$timestamp[ o@data$vo == whatVO ];
		x <- aggregate(
			t,
			by = list( as.Date( jan1970 + t ) ),
			length
		)

		return( list( date = levels( x$Group ), njobs = x$x ) )
	}
);

## Overall Number of Jobs (all days) histogram
setGeneric(
	"mg_lcg_overallNJobsByDateHistPlot",
	function( o ) { standardGeneric( "mg_lcg_overallNJobsByDateHistPlot" ) }
);
setMethod(
	"mg_lcg_overallNJobsByDateHistPlot",
	"mg_lcg",
	function( o )
	{
		x <- mg_lcg_overallNJobsByDate(o);
		res <- mg_plot_bar(
			x$njobs,
			bar.names = format.Date(as.Date(x$date),format="%b-%d"),
			xlab = "day",
			ylab = "# jobs"
		);
		return( invisible(res) );
	}
);

## Number of Jobs by VO histogram
setGeneric(
	"mg_lcg_njobsByVOHistPlot",
	function( o ) { standardGeneric( "mg_lcg_njobsByVOHistPlot" ) }
);
setMethod(
	"mg_lcg_njobsByVOHistPlot",
	"mg_lcg",
	function( o )
	{
		x <- mg_lcg_njobsByVO(o);
		op <- par( las = 1 );
		ord <- order( x$vo, decreasing=TRUE );
		res <- mg_plot_bar(
			x$njobs[ord],
			x$vo[ord],
			xlab = "# jobs",
			ylab = NULL,
			show.polygon = FALSE,
			horiz = TRUE
		);
		par( op  );
		return( invisible(res) );
	}
);

## Number of Jobs by User histogram
setGeneric(
	"mg_lcg_njobsByUserHistPlot",
	function( o ) { standardGeneric( "mg_lcg_njobsByUserHistPlot" ) }
);
setMethod(
	"mg_lcg_njobsByUserHistPlot",
	"mg_lcg",
	function( o )
	{
		x <- mg_lcg_njobsByUser(o);
		op <- par( las = 1 );
		ord <- order( as.numeric(x$uid) );
		res <- mg_plot_bar(
			x$njobs[ord],
			x$uid[ord],
			xlab = "user id",
			ylab = "# jobs",
			show.polygon = FALSE,
			horiz = FALSE
		);
		par( op  );
		return( invisible(res) );
	}
);

## VO Number of jobs (by date) Histogramm
setGeneric(
	"mg_lcg_voNJobsByDatePlot",
	function( o, whatVO ) { standardGeneric( "mg_lcg_voNJobsByDatePlot" ) }
);
setMethod(
	"mg_lcg_voNJobsByDatePlot",
	"mg_lcg",
	function( o, whatVO )
	{
		x <- mg_lcg_voNJobsByDate( o, whatVO );
		res <- mg_plot_bar(
			x$njobs,
			bar.names = format.Date(as.Date(x$date),format="%b-%d"),
			xlab = "day",
			ylab = "# jobs",
			show.polygon = TRUE
		);
		return( invisible(res) );
	}
);


####
## Interarrival time
####

## Overall Interarrival Times
setGeneric(
	"mg_lcg_overallInterArrTimes",
	function( o ) { standardGeneric( "mg_lcg_overallInterArrTimes" ) }
);
setMethod(
	"mg_lcg_overallInterArrTimes",
	"mg_lcg",
	function( o )
	{
		t <- sort.int( o@data$timestamp )

		for ( i in seq( length(t), 2 ) )
		{
			t[i] <- t[i] - t[i-1]
		}
		t[1] <- 0

		return ( t )
	}
);

## VO Interarrival Times
setGeneric(
	"mg_lcg_voInterArrTimes",
	function( o, whatVO ) { standardGeneric( "mg_lcg_voInterArrTimes" ) }
);
setMethod(
	"mg_lcg_voInterArrTimes",
	"mg_lcg",
	function( o, whatVO )
	{
		##t <- subset( o@data, vo == whatVO, select = timestamp, drop = TRUE )
		t <- subset( o@data, vo == whatVO, select = timestamp );
		t <- sort.int( t$timestamp );

		#firstArrTs <- min( t$timestamp );
		##firstArrTs <- min( t );
		##t <- transform( t, timestamp = timestamp - firstArrTs );
		#f <- function(y) { return( y - firstArrTs ) }
		#t <- lapply( t, f );
		#t <- t$timestamp; # to vector
		#
		#return( t )

		n <- length( t );

		if ( n > 1 )
		{
			for ( i in seq( length(t), 2 ) )
			{
				t[i] <- t[i] - t[i-1];
			}
		}
		t[1] <- 0;

		return( t );
	}
);

####
## Runtime
####

## Overall Runtimes
setGeneric(
	"mg_lcg_overallRuntimes",
	function( o ) { standardGeneric( "mg_lcg_overallRuntimes" ) }
);
setMethod(
	"mg_lcg_overallRuntimes",
	"mg_lcg",
	function( o )
	{
		return( o@data$runtime );
	}
);

## VO Runtimes
setGeneric(
	"mg_lcg_voRuntimes",
	function( o, whatVO ) { standardGeneric( "mg_lcg_voRuntimes" ) }
);
setMethod(
	"mg_lcg_voRuntimes",
	"mg_lcg",
	function( o, whatVO )
	{
		return( o@data$runtime[ o@data$vo == whatVO ] );
	}
);
