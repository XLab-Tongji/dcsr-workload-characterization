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

## MG_TERAGRID
##
## SUMMARY
##   TeraGrid class.
##
## DESCRIPTION
##   Class for TeraGrid log.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( 'lib/mg/mg_plot.R' )

## Class representation
mg_teragrid <- setClass(
	'mg_teragrid',
	slots = ( data='data.frame' )
);

## Constructor
setMethod(
	'initialize',
	'mg_teragrid',
	function( .Object, file, omitAgtS = FALSE, omitSeqE = FALSE, omitSgtE = FALSE, omitNle0 = FALSE, fixAgtS = FALSE )
	{
		#.Object <- callNextMethod()
		.Object@data <- read.table( file, header = TRUE )

		if ( omitAgtS )
		{
			.Object@data <- subset( .Object@data, timestamp <= starttime );
		} else if ( fixAgtS ) {
			d1 <- subset( .Object@data, timestamp <= starttime );
			d2 <- subset( .Object@data, timestamp > starttime );
			d2[ c("starttime","endtime") ] <- c( d2$timestamp, d2$timestamp + ( d2$endtime - d2$starttime ) );
			.Object@data <- rbind( d1, d2 );

			d1 <- d2 <- NULL; # clean-up
		}
		if ( omitSeqE )
		{
			.Object@data <- subset( .Object@data, starttime != endtime );
		}
		if ( omitSgtE )
		{
			.Object@data <- subset( .Object@data, starttime <= endtime );
		}
		if ( omitNle0 )
		{
			.Object@data <- subset( .Object@data, nodes > 0 );
		}

		# sort the data.frame by timestamp
		.Object@data <- .Object@data[ order( .Object@data$timestamp ), ];

		return( .Object )
	}
);

####
## Number of Jobs
####

## Overall number of jobs by date (no time)
setGeneric(
	"mg_teragrid_overallNJobsByDate",
	function( o ) { standardGeneric( "mg_teragrid_overallNJobsByDate" ) }
);
setMethod(
	"mg_teragrid_overallNJobsByDate",
	"mg_teragrid",
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

## Number of jobs per nodes
setGeneric(
	"mg_teragrid_njobsByNodes", 
	function( o ) { standardGeneric( "mg_teragrid_njobsByNodes" ) }
);
setMethod(
	"mg_teragrid_njobsByNodes",
	"mg_teragrid", 
	function( o )
	{
		x <- aggregate(
			o@data$nodes,
			by = list( o@data$nodes ),
			length
		)

		return( list( nodes = levels( x$Group ), njobs = x$x ) )
	}
);

## Overall Number of Jobs (all days) histogram
setGeneric(
	"mg_teragrid_overallNJobsByDateHistPlot",
	function( o ) { standardGeneric( "mg_teragrid_overallNJobsByDateHistPlot" ) }
);
setMethod(
	"mg_teragrid_overallNJobsByDateHistPlot",
	"mg_teragrid",
	function( o )
	{
		x <- mg_teragrid_overallNJobsByDate(o);
		res <- mg_plot_bar(
			x$njobs,
			bar.names = format.Date(as.Date(x$date),format="%b-%d"),
			xlab = "day",
			ylab = "# jobs"
		);
		return( invisible(res) );
	}
);

## Number of Jobs by Nodes histogram
setGeneric(
	"mg_teragrid_njobsByNodesHistPlot",
	function( o ) { standardGeneric( "mg_teragrid_njobsByNodesHistPlot" ) }
);
setMethod(
	"mg_teragrid_njobsByNodesHistPlot",
	"mg_teragrid",
	function( o )
	{
		x <- mg_teragrid_njobsByNodes(o);
		op <- par( las = 1 );
		ord <- order( as.numeric(x$nodes) );
		res <- mg_plot_bar(
			x$njobs[ord],
			x$nodes[ord],
			xlab = "nodes",
			ylab = "# jobs",
			show.polygon = FALSE,
			horiz = FALSE
		);
		par( op  );
		return( invisible(res) );
	}
);

####
## Interarrival time
####

## Overall Interarrival Times
setGeneric(
	"mg_teragrid_overallInterArrTimes",
	function( o ) { standardGeneric( "mg_teragrid_overallInterArrTimes" ) }
);
setMethod(
	"mg_teragrid_overallInterArrTimes",
	"mg_teragrid",
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

####
## Runtime
####

## Overall Runtimes
setGeneric(
	"mg_teragrid_overallRuntimes",
	function( o ) { standardGeneric( "mg_teragrid_overallRuntimes" ) }
);
setMethod(
	"mg_teragrid_overallRuntimes",
	"mg_teragrid",
	function( o )
	{
		return( o@data$endtime - o@data$starttime );
	}
);

##### BEGIN DEPRECATED METHODS #####

####
## Number of Nodes
####

#### Node Trend plot.
#setGeneric(
#	'mg_nodesTrendPlot',
#	function( o ) { standardGeneric( 'mg_nodesTrendPlot' ) }
#);
#setMethod(
#	'mg_nodesTrendPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_trend(
#			o@data$timestamp,
#			o@data$nodes,
#			title = 'TeraGrid Jobs Size - Trend',
#			xlabel = 'Time',
#			ylabel = 'Node',
#		);
#	}
#);
#
### Node CDF plot.
#setGeneric(
#	'mg_nodesCDFPlot',
#	function( o ) { standardGeneric( 'mg_nodesCDFPlot' ) }
#);
#setMethod(
#	'mg_nodesCDFPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_ecdf(
#			o@data$nodes,
#			title = 'TeraGrid Jobs Size - CDF',
#			xlabel = 'x [# nodes]',
#			ylabel = 'Fn(x)'
#		)
#	}
#);
#
### Node CCDF plot.
#setGeneric(
#	'mg_nodesCCDFPlot',
#	function( o ) { standardGeneric( 'mg_nodesCCDFPlot' ) }
#);
#setMethod(
#	'mg_nodesCCDFPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_ccdf(
#			o@data$nodes,
#			title = 'TeraGrid Jobs Size - CCDF',
#			xlabel = 'x [# nodes]',
#			ylabel = '1-Fn(x)'
#		)
#	}
#);
#
### Node LLCD plot.
#setGeneric(
#	'mg_nodesLLCDPlot',
#	function( o ) { standardGeneric( 'mg_nodesLLCDPlot' ) }
#);
#setMethod(
#	'mg_nodesLLCDPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_llcd(
#			o@data$nodes,
#			title = 'TeraGrid Jobs Size - LLCD',
#			xlabel = 'x [# nodes]',
#			ylabel = 'Fn(x)'
#		)
#	}
#);
#
### Nodes Histogram plot.
#setGeneric(
#	'mg_nodesHistPlot',
#	function( o ) { standardGeneric( 'mg_nodesHistPlot' ) }
#);
#setMethod(
#	'mg_nodesHistPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_hist(
#			o@data$nodes,
#			title = 'TeraGrid Jobs Size - Histogram',
#			xlabel = 'x [# nodes]',
#			ylabel = 'Frequency'
#		)
#	}
#);

### Runtime Trend plot.
#setGeneric(
#	'mg_runtimeTrendPlot',
#	function( o ) { standardGeneric( 'mg_runtimeTrendPlot' ) }
#);
#setMethod(
#	'mg_runtimeTrendPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_trend(
#			o@data$timestamp,
#			o@data$runtime,
#			title = 'TeraGrid Jobs Runtime - Trend',
#			xlabel = 'Time',
#			ylabel = 'Runtime [s]'
#		);
#	}
#);
#
### Runtime CDF plot.
#setGeneric(
#	'mg_runtimeCDFPlot',
#	function( o ) { standardGeneric( 'mg_runtimeCDFPlot' ) }
#);
#setMethod(
#	'mg_runtimeCDFPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_ecdf(
#			o@data$runtime,
#			title = 'TeraGrid Jobs Runtime - CDF',
#			xlabel = 'x [s]',
#			ylabel = 'Fn(x)'
#		)
#	}
#);
#
### Runtime CCDF plot.
#setGeneric(
#	'mg_runtimeCCDFPlot',
#	function( o ) { standardGeneric( 'mg_runtimeCCDFPlot' ) }
#);
#setMethod(
#	'mg_runtimeCCDFPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_ecdf(
#			o@data$runtime,
#			title = 'TeraGrid Jobs Runtime - CCDF',
#			xlabel = 'x [s]',
#			ylabel = '1-Fn(x)'
#		)
#	}
#);
#
### Runtime LLCD plot.
#setGeneric(
#	'mg_runtimeLLCDPlot',
#	function( o ) { standardGeneric( 'mg_runtimeLLCDPlot' ) }
#);
#setMethod(
#	'mg_runtimeLLCDPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_llcd(
#			o@data$runtime,
#			title = 'TeraGrid Jobs Runtime - LLCD',
#			xlabel = 'x [s]',
#			ylabel = '1-Fn(x)'
#		)
#	}
#);
#
### Runtime Histogram plot.
#setGeneric(
#	'mg_runtimeHistPlot',
#	function( o ) { standardGeneric( 'mg_runtimeHistPlot' ) }
#);
#setMethod(
#	'mg_runtimeHistPlot',
#	'mg_teragrid',
#	function( o )
#	{
#		mg_plot_hist(
#			o@data$runtime,
#			title = 'TeraGrid Jobs Runtime - Histogram',
#			xlabel = 'x [s]',
#			ylabel = 'Frequency'
#		)
#	}
#);

##### END DEPRECATED METHODS #####
