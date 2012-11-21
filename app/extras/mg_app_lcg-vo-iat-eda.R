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

## MG_APP_LCG-VO-IAT-ALL
##
## SUMMARY
##  LCG - Overall interarrival times.
##
## DESCRIPTION
##  Writes:
##  * a data-file containing the LCG interarrival times for each LCG VO.
##  * a stat (ascii) file containing varoius statistics.
##  Plots:
##  * a frequency histogram;
##  * a density histogram;
##  * a time-series plot;
##  * an LLCD plot with several fitting distributions;
##  * a QQ-plot;
##  * an ACF plot;
##  * a Hill estimator plot.
##  * ...
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

oldwarn <- getOption( "warn" );
options( warn = 1 );

source( "app/mg_app_conf.R" );
source( "lib/mg/mg_debug.R" );
source( "lib/mg/mg_dists.R" );
source( "lib/mg/mg_eda.R" );
source( "lib/mg/mg_fit.R" );
source( "lib/mg/mg_gof.R" );
source( "lib/mg/mg_lcg.R" );
source( "lib/mg/mg_lrd.R" );
source( "lib/mg/mg_moments.R" );
source( "lib/mg/mg_plot.R" );
source( "lib/mg/mg_statsutils.R" );

set.seed(123456); # used to eventually repeat the same experiment

dbg <- new( "mg_debug", activate = TRUE );

lcg <- new( "mg_lcg", MG_APP_CONF_LCG_LOGFULLNAME );

#voNames <- mg_lcg_voNames( lcg );
#voNames <- c( "alice" );#, "atlas", "babar" );#, "cms", "lhcb" );#XXX
voNames <- c( "alice", "lhcb", "dteam", "cms", "atlas", "babar", "zeus" );

for ( vo in voNames )
{
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Analyzing VO: '", vo, "'" );

	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Reading IAT file" );

	iat <- mg_lcg_voInterArrTimes( lcg, vo );

	sink( paste( MG_APP_CONF_LCG_STATSPATH, "/lcg-vo_", vo, "-iat-stats-eda.out", sep = "" ) );

	cat( "\n" );
	cat( paste( rep( "#", 80 ), collapse = "" ), "\n" );
	cat( paste( ">>> VO: '", vo, "' <<<", sep = "" ), "\n" );
	cat( paste( rep( "#", 80 ), collapse = "" ), "\n" );
	cat( "\n" );

	##@{ EDA shape plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Plotting EDA Shape Plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-eda-shape.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_eda.shape.plot( iat )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-eda-shape.eps", sep = "" ) );
		warning( "Error while drawing EDA shape plot for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
	}
	##@} EDA shape plot

	# BEGIN Writes samples statistics to data file
	ci.mean.lev <- ci.var.lev <- ci.median.lev <- 0.95;
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Writing summary statistics" );
	cat( "## Summary statistics ##", "\n" );
	cat( paste( "# Observations: ", length( iat ), sep = "" ), "\n" );
	#res <- summary( iat );
	qq <- quantile(iat, probs = c( 0, 0.25, 0.75, 1 ) ); # get min, 1st Quartile, 3rd Quartile, max
	cat( paste( "Min:", qq[1] ), "\n" );
	cat( paste( "1st. Quartile:", qq[2] ), "\n" );
	ci <- mg_statsutils_confint.median( iat, level = ci.median.lev );
	cat( paste( "Median: ", ci[["median"]], " with (", ci[["lower"]], ",", ci[["upper"]], ") @ ", (100*ci.median.lev), "% level", sep = "" ), "\n" );
	cat( paste( "3rd. Quartile:", qq[3] ), "\n" );
	cat( paste( "Max:", qq[4] ), "\n" );
	cat( paste( "IQR:", IQR(iat) ), "\n" );
	cat( paste( "MAD:", mad(iat) ), "\n" );
	cat( paste( "Skewness:", mg_moments_skewness(iat) ), "\n" );
	cat( paste( "Quartile coefficient of Skewness:", mg_moments_skewness.quartile(iat) ), "\n" );
	cat( paste( "Kurtosis:", mg_moments_kurtosis(iat) ), "\n" );
	cat( paste( "Kurtosis Excess:", mg_moments_kurtosis.excess(iat) ), "\n" );
	ci <- mg_statsutils_confint.mean( iat, level = ci.mean.lev );
	cat( paste( "Mean: ", ci[["mean"]], " with (", ci[["lower"]], ",", ci[["upper"]], ") @ ", (100*ci.mean.lev), "% level", sep = "" ), "\n" );
	ci <- mg_statsutils_confint.var( iat, level = ci.var.lev );
	cat( paste( "Std Err: ", sqrt(ci[["var"]]), " with (", sqrt(ci[["lower"]]), ",", sqrt(ci[["upper"]]), ") @ ", (100*ci.var.lev), "% level", sep = "" ), "\n" );
	cat( paste( "CV:", mg_statsutils_cv( iat ) ), "\n" );

	cat( "\n" );
	# END Writes samples statistics to data file

	# BEGIN Plot trend
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Plotting trend" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-trend.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_trend( iat )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-trend.eps", sep = "" ) );
		warning( "Error while drawing trend plot for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
	}
	# END Plot trend

	# BEGIN Plot frequency histogram
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Plotting frequency histogram" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-fhist.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_freqHist( iat, title = NULL )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-fhist.eps", sep = "" ) );
		warning( "Error while drawing frequency histogram plot for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
	}
	# END Plot frequency histogram

	# BEGIN Plot density histogram
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Plotting density histogram" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-dhist.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_densHist( iat, title = NULL )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-dhist.eps", sep = "" ) );
		warning( "Error while drawing density histogram plot for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
	}
	# END Plot density histogram

	# BEGIN Plot lag-1 plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Plotting Lag-1 plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo-iat-lag1.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_lag( iat, lag = 1 )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo-iat-lag1.eps", sep = "" ) );
		warning( "Error while drawing Lag-1 plot." );
		warning( "Error message is: ", geterrmessage() );
	}
	# END Plot lag-1 plot

	# BEGIN Plot and print ACF
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Drawing ACF Plot" );
	cat( "## ACF statistics and plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-acf.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		acf( iat )
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		print( res );
	} else {
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-acf.eps", sep = "" ) );
		warning( "Error while drawing ACF plot for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) ACF failed.", "\n" );
	}
	# END Plot and print ACF

	##@{ EDA autocorrelation plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Plotting EDA Autocorrelation Plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-eda-acor.eps", sep = "" ), horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_eda.acor.plot( iat )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-eda-acor.eps", sep = "" ) );
		warning( "Error while drawing EDA Autocorrelation plot for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
	}
	##@} EDA autocorrelation plot

	##@{ Plot and print Long Range Dependence indicators

	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Drawing Long Range Dependence indicators Plot" );
	cat( "## Long Range Dependence indicators statistics and plots:", "\n" );

	##@{ Aggregated trend plot
	cat( "## LRD Aggregated Trends plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-aggtrend.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		#mg_plot_aggregation.trend( iat )
		mg_plot_aggregation.trend(
			iat,
			scale.start = 1,
			scale.end = ceiling( log10( length(iat) ) ) * 10 + 1,
			scale.incr = 0,
			scale.mult = 10,
			fit.method = "enlarge"
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		print( res );
	} else {
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-aggtrend.eps", sep = "" ) );
		warning( "Error while plotting Aggregated Trends for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Aggregated Trends plot failed.", "\n" );
	}
	##@} Aggregated trend plot

	##@{ Aggregated Log-Log CDF plot
	cat( "## LRD Aggregated Log-Log Empirical CDF plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-aggllecdf.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_aggregation.ecdf(
			iat,
			scale.end = ceiling( log2( length(iat) ) ) * 2 + 1,
			scale.incr = 0,
			scale.mult = 2,
			fit.method = "enlarge",
			log = "xy"
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		print( res );
	} else {
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-aggllecdf.eps", sep = "" ) );
		warning( "Error while plotting Aggregated Log-Log Empirical CDF for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Aggregated Log-Log Empirical CDF plot failed.", "\n" );
	}
	##@} Aggregated Log-Log CDF plot

	##@{ Aggregated Log-Log Empirical CCDF plot
	cat( "## LRD Aggregated Log-Log Empirical CCDF plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-agglleccdf.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_aggregation.eccdf(
			iat,
			scale.end = ceiling( log2( length(iat) ) ) * 2 + 1,
			scale.incr = 0,
			scale.mult = 2,
			fit.method = "enlarge",
			log = "xy"
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		print( res );
	} else {
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-agglleccdf.eps", sep = "" ) );
		warning( "Error while plotting Aggregated Log-Log CCDF for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Aggregated Log-Log Empirical CCDF plot failed.", "\n" );
	}
	##@} Aggregated Log-Log CDF plot

	##@{ Aggregate of Variance method (for Hurst exponent)
	cat( "## LRD Aggregate of Variance method:", "\n" );
	res <- try(
		mg_lrd_aggvar( iat, ret.log.m = TRUE, ret.log.vars = TRUE )
	);
	if ( !is( res, "try-error" ) )
	{
		print( res );

		postscript(
			file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-av.eps", sep = "" ),
			horizontal = FALSE,
			onefile = FALSE
		);
		res <- try(
			mg_plot_lrd.aggvar( res, plot.fitline = TRUE )
		);
		dev.off();
		if ( is( res, "try-error" ) )
		{
			unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-av.eps", sep = "" ) );
			warning( "Error while plotting Aggregate of Variance method for VO '", vo, "'." );
			warning( "Error message is: ", geterrmessage() );
			cat( "(ERROR) Aggregate of Variance method failed.", "\n" );
		}
	} else {
		warning( "Error while performing Aggregate of Variance method for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Aggregate of Variance method failed.", "\n" );
	}
	##@} Aggregate of Variance method (for Hurst exponent)

	##@{ R/S method (for Hurst exponent)
	cat( "## LRD R/S method:", "\n" );
	res <- try(
		mg_lrd_rs( iat, ret.log.m = TRUE, ret.log.rs = TRUE )
	);
	if ( !is( res, "try-error" ) )
	{
		print( res );

		postscript(
			file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-rs.eps", sep = "" ),
			horizontal = FALSE,
			onefile = FALSE
		);
		res <- try(
			mg_plot_lrd.rs( res, plot.fitline = TRUE )
		);
		dev.off();
		if ( is( res, "try-error" ) )
		{
			unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-rs.eps", sep = "" ) );
			warning( "Error while plotting R/S method for VO '", vo, "'." );
			warning( "Error message is: ", geterrmessage() );
			cat( "(ERROR) R/S method failed.", "\n" );
		}
	} else {
		warning( "Error while performing R/S method for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) R/S method failed.", "\n" );
	}
	##@} R/S method (for Hurst exponent)

	##@{ Periodogram method (for Hurst exponent)
	cat( "## LRD Periodogram method:", "\n" );
	res <- try(
		mg_lrd_periodogram( iat, ret.log.freqs = TRUE, ret.log.pgrams = TRUE )
	);
	if ( !is( res, "try-error" ) )
	{
		print( res );

		postscript(
			file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-pgram.eps", sep = "" ),
			horizontal = FALSE,
			onefile = FALSE
		);
		res <- try(
			mg_plot_lrd.periodogram( res, plot.fitline = TRUE )
		);
		dev.off();
		if ( is( res, "try-error" ) )
		{
			unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-pgram.eps", sep = "" ) );
			warning( "Error while plotting Periodogram method for VO '", vo, "'." );
			warning( "Error message is: ", geterrmessage() );
			cat( "(ERROR) Periodogram method failed.", "\n" );
		}
	} else {
		warning( "Error while performing Periodogram method for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Periodogram method failed.", "\n" );
	}
	##@} Periodogram method (for Hurst exponent)

	##@{ Cumulative Periodogram method (for Hurst exponent)
	cat( "## LRD Cumulative Periodogram method:", "\n" );
	res <- try(
		mg_lrd_periodogram.cum( iat, ret.log.freqs = TRUE, ret.log.pgrams = TRUE )
	);
	if ( !is( res, "try-error" ) )
	{
		print( res );

		postscript(
			file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-cumpgram.eps", sep = "" ),
			horizontal = FALSE,
			onefile = FALSE
		);
		res <- try(
			mg_plot_lrd.periodogram( res, plot.fitline = TRUE )
		);
		dev.off();
		if ( is( res, "try-error" ) )
		{
			unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-lrd-cumpgram.eps", sep = "" ) );
			warning( "Error while plotting Cumulative Periodogram method for VO '", vo, "'." );
			warning( "Error message is: ", geterrmessage() );
			cat( "(ERROR) Cumulative Periodogram method failed.", "\n" );
		}
	} else {
		warning( "Error while performing Cumulative Periodogram method for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Cumulative Periodogram method failed.", "\n" );
	}
	##@} Cumulative Periodogram method (for Hurst exponent)

	##@} Plot and print Long Range Dependence indicators

	##@{ Plot Hill estimators
	mg_debug.message( dbg, "[MG::APP::VO::LCG::IAT>> Drawing Hill Plot" );
	cat( "## Hill estimator:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-hill.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_hill( iat, conf.level = 0.95 )
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		#cat( "Hill estimators: " );
		#print( res );
		#cat( "Last Hill estimator: " );
		#print( res$alpha[ length( res$alpha ) ] );
		cat( paste( "Last Hill estimator: ",  res$statistic[ length( res$statistic ) ], " with (",  res$conf.int$lower[ length( res$conf.int$lower ) ], ", ", res$conf.int$upper[ length( res$conf.int$upper ) ], ") @ 95% level", sep = "" ), "\n" );
	} else {
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-iat-hill.eps", sep = "" ) );
		warning( "Error while drawing Hill plot for VO '", vo, "'." );
		warning( "Error message is: ", geterrmessage(), sep = "" )
		cat( "(ERROR) Hill estimator failed.", "\n" );
	}
	##@} Plot Hill estimators

	cat( paste( "<<< VO: '", vo, "' >>>", sep = "" ), "\n" );

	sink( file = NULL );
}

options( warn = oldwarn );
