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

## MG_APP_TERAGRID-OVERALL-IAT-ALL
##
## SUMMARY
##  TeraGrid - Overall interarrival times.
##
## DESCRIPTION
##  Writes:
##  * a data-file containing the overall interarrival times for TERAGRID log;
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
source( "lib/mg/mg_teragrid.R" );
source( "lib/mg/mg_lrd.R" );
source( "lib/mg/mg_moments.R" );
source( "lib/mg/mg_plot.R" );
source( "lib/mg/mg_statsutils.R" );

set.seed(123456); # used to eventually repeat the same experiment

dbg <- new( "mg_debug", activate = TRUE );

#teragrid <- new( "mg_teragrid", MG_APP_CONF_TERAGRID_LOGFULLNAME );
teragrid <- new(
	"mg_teragrid",
	MG_APP_CONF_TERAGRID_LOGFULLNAME,
	omitSeqE = TRUE,
	omitSgtE = TRUE,
	omitNle0 = TRUE,
	fixAgtS = TRUE
);

mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Reading IAT file" );

iat <- mg_teragrid_overallInterArrTimes( teragrid );

# BEGIN Plot fitted distributions and writes estimated parameters to data file.
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Fitting distributions" );
fitDistrs <- mg_fit_methodsNameForDistrs(
	c(
#		MG_CONSTS_BETA_DIST,
		MG_CONSTS_CAUCHY_DIST,
#		MG_CONSTS_CHISQUARED_DIST,
		MG_CONSTS_EXPONENTIAL_DIST,
#		MG_CONSTS_FISHERF_DIST,
#		#MG_CONSTS_FRECHET_DIST,
		MG_CONSTS_GAMMA_DIST,
		MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST,
		MG_CONSTS_GENERALIZEDPARETO_DIST,
#		#MG_CONSTS_GUMBEL_DIST,
#		#MG_CONSTS_HYPEREXPONENTIAL_DIST,
		MG_CONSTS_LOGISTIC_DIST,
		MG_CONSTS_LOGNORMAL_DIST,
#		#MG_CONSTS_MMPP_DIST,
		MG_CONSTS_NORMAL_DIST,
		MG_CONSTS_PARETO_DIST,
		MG_CONSTS_CONTINUOUSPHASETYPE_DIST,
#		#MG_CONSTS_REVERSEWEIBULL_DIST,
#		MG_CONSTS_SKEWEDSTABLE_DIST,
#		MG_CONSTS_STUDENTT_DIST,
#		MG_CONSTS_SYMSTABLE_DIST,
		MG_CONSTS_WEIBULL_DIST
	)
);
if ( length( fitDistrs ) == 0 )
{
	warning( "Cannot perform distributions fitting." );
	cat( "(WARNING) No method available for performing distributions fitting.", "\n" );
	quit();
}
#postscript(
#	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llcd.eps", sep = "" ),
#	#title = "TERAGRID - Overall Interarrival Times",
#	horizontal = FALSE,
#	onefile = FALSE
#);
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Distribution multi fitting" );
fit <- try(
	mg_fit_multifit(
		iat,
		fitDistrs
		#plotType = MG_CONSTS_LLCD_PLOT
	)
);
#dev.off();
if ( is( fit, "try-error" ) )
{
	#unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llcd.eps", sep = "" ) );
	warning( "Error while performing distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}
cat( "## Distributions Fitting:", "\n" );
print( fit );

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
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing a CCDF plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-ccdf.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_multiplot(
		iat,
		fit.est,
		plot.type = MG_CONSTS_CCDF_PLOT
	)
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-ccdf.eps", sep = "" ) );
	warning( "Error while plotting CCDF distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}

# Draw a LLCD (LLCDF) plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing a LLCCDF plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llccdf.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_multiplot(
		iat,
		fit.est,
		plot.type = MG_CONSTS_LLCD_PLOT
	)
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llccdf.eps", sep = "" ) );
	warning( "Error while plotting LLCD distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}

# Draw a CDF plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing a CDF plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-cdf.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_multiplot(
		iat,
		fit.est,
		plot.type = MG_CONSTS_CDF_PLOT
	)
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-cdf.eps", sep = "" ) );
	warning( "Error while plotting CDF distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}

# Draw a LLCDF plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing a LLCDF plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llcdf.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_multiplot(
		iat,
		fit.est,
		plot.type = MG_CONSTS_LLCDF_PLOT
	)
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llcdf.eps", sep = "" ) );
	warning( "Error while plotting LLCDF distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}

# Draw a PDF plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing a PDF plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-pdf.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_multiplot(
		iat,
		fit.est,
		plot.type = MG_CONSTS_PDF_PLOT
	)
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-pdf.eps", sep = "" ) );
	warning( "Error while plotting PDF distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}

##@} Draw some plots

options( warn = oldwarn );
