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

## MG_APP_LCG-VO-RT-ALL
##
## SUMMARY
##  LCG - Overall runtimes.
##
## DESCRIPTION
##  Writes:
##  * a data-file containing the LCG runtimes for each LCG VO.
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
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Analyzing VO: '", vo, "'" );

	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Reading RT file" );

	rt <- mg_lcg_voRuntimes( lcg, vo );

	##@{ Distribution fitting (parameters estimated from data set)
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Fitting distributions" );
	fitDistrs <- mg_fit_methodsNameForDistrs(
		c(
#			MG_CONSTS_BETA_DIST,
			MG_CONSTS_CAUCHY_DIST,
#			MG_CONSTS_CHISQUARED_DIST,
			MG_CONSTS_EXPONENTIAL_DIST,
#			MG_CONSTS_FISHERF_DIST,
#			MG_CONSTS_FRECHET_DIST,
			MG_CONSTS_GAMMA_DIST,
			MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST,
			MG_CONSTS_GENERALIZEDPARETO_DIST,
#			#MG_CONSTS_GUMBEL_DIST,
#			#MG_CONSTS_HYPEREXPONENTIAL_DIST,
			MG_CONSTS_LOGISTIC_DIST,
			MG_CONSTS_LOGNORMAL_DIST,
#			#MG_CONSTS_MMPP_DIST,
			MG_CONSTS_NORMAL_DIST,
			MG_CONSTS_PARETO_DIST,
			MG_CONSTS_CONTINUOUSPHASETYPE_DIST,
#			MG_CONSTS_REVERSEWEIBULL_DIST,
#			MG_CONSTS_SKEWEDSTABLE_DIST,
#			MG_CONSTS_STUDENTT_DIST,
#			MG_CONSTS_SYMSTABLE_DIST,
			MG_CONSTS_WEIBULL_DIST
		)
	);
	if ( length( fitDistrs ) == 0 )
	{
		warning( "Cannot perform distributions fitting for VO: '", vo, "'." );
		cat( paste( "(WARNING) No method available for performing distributions fitting for VO '", vo, "'.", sep = "" ), "\n" );
		next;
	}
#	postscript(
#		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-llcd.eps", sep = "" ),
#		#title = paste( "LCG - VO '", vo, "' Interarrival Times", sep = "" ),
#		horizontal = FALSE,
#		onefile = FALSE
#	);
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Distribution multi fitting" );
	fit <- try(
		mg_fit_multifit(
			rt,
			fitDistrs
			#plotType = MG_CONSTS_LLCD_PLOT
		)
	);
#	dev.off();
	if ( is( fit, "try-error" ) )
	{
		#unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-llcd.eps", sep = "" ) );
		warning( "Error while performing distributions fitting for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		next;
	}
	cat( paste( "## Distributions Fitting for VO '", vo, "':", sep = "" ), "\n" );
	print( fit );

	fit.est <- NULL; # clean-up mem

	##@{ Draw some plots

	# Rearranges the data structure of estimated parameters
	fit.est <- list();
	for ( d in names(fit) )
	{
		if ( !is.null( fit[[ d ]] ) && length( fit[[ d ]] ) != 0 )
		{
			fit.est[[ d ]] <- fit[[ d ]]$estimate;
		}
	}

	# Draw a CCDF plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Drawing a CCDF plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-ccdf.eps", sep = "" ),
		#title = paste( "LCG - VO '", vo, "' Interarrival Times", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_multiplot(
			rt,
			fit.est,
			plot.type = MG_CONSTS_CCDF_PLOT
		)
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-ccdf.eps", sep = "" ) );
		warning( "Error while plotting CCDF distributions fitting for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		next;
	}

	# Draw a LLCD (LLCDF) plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Drawing a LLCCDF plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-llccdf.eps", sep = "" ),
		#title = paste( "LCG - VO '", vo, "' Interarrival Times", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_multiplot(
			rt,
			fit.est,
			plot.type = MG_CONSTS_LLCD_PLOT
		)
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-llccdf.eps", sep = "" ) );
		warning( "Error while plotting LLCD distributions fitting for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		next;
	}

	# Draw a CDF plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Drawing a CDF plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-cdf.eps", sep = "" ),
		#title = paste( "LCG - VO '", vo, "' Interarrival Times", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_multiplot(
			rt,
			fit.est,
			plot.type = MG_CONSTS_CDF_PLOT
		)
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-cdf.eps", sep = "" ) );
		warning( "Error while plotting CDF distributions fitting for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		next;
	}

	# Draw a LLCDF plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Drawing a LLCDF plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-llcdf.eps", sep = "" ),
		#title = paste( "LCG - VO '", vo, "' Interarrival Times", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_multiplot(
			rt,
			fit.est,
			plot.type = MG_CONSTS_LLCDF_PLOT
		)
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-llcdf.eps", sep = "" ) );
		warning( "Error while plotting LLCDF distributions fitting for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		next;
	}

	# Draw a PDF plot
	mg_debug.message( dbg, "[MG::APP::VO::LCG::RT>> Drawing a PDF plot" );
	postscript(
		file = paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-pdf.eps", sep = "" ),
		#title = paste( "LCG - VO '", vo, "' Interarrival Times", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_multiplot(
			rt,
			fit.est,
			plot.type = MG_CONSTS_PDF_PLOT
		)
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_LCG_IMGPATH, "/lcg-vo_", vo, "-rt-fit-pdf.eps", sep = "" ) );
		warning( "Error while plotting PDF distributions fitting for VO: '", vo, "'." );
		warning( "Error message is: ", geterrmessage() );
		next;
	}

	##@} Draw some plots

	##@} Distribution fitting (parameters estimated from data set)

	cat( paste( "<<< VO: '", vo, "' >>>", sep = "" ), "\n" );
}

options( warn = oldwarn );
