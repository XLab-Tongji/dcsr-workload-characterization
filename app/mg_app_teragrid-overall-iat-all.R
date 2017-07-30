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

olderror <- getOption( "error" );
options( error = quote( dump.frames( "mg_app_teragrid-overall-iat-all", TRUE ) ) ); #DEBUG

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

dbg <- mg_debug( activate = TRUE );

#teragrid <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME );
teragrid <- mg_teragrid(
	MG_APP_CONF_TERAGRID_LOGFULLNAME,
	omitSeqE = TRUE,
	omitSgtE = TRUE,
	omitNle0 = TRUE,
	fixAgtS = TRUE
);

mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Reading IAT file" );

iat <- mg_teragrid_overallInterArrTimes( teragrid );

# BEGIN Writes samples to data file
write(
	iat,
	paste( MG_APP_CONF_TERAGRID_VIEWPATH, "/teragrid-overall-iat.data", sep = "" ),
	ncolumns = 1
);
# END Writes samples to data file

sink( file = paste( MG_APP_CONF_TERAGRID_STATSPATH, "/teragrid-overall-iat-stats.out", sep = "" ) );

##@{ EDA shape plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Plotting EDA Shape Plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-eda-shape.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_eda.shape.plot( iat )
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-eda-shape.eps", sep = "" ) );
	warning( "Error while drawing EDA shape plot." );
	warning( "Error message is: ", geterrmessage() );
}
##@} EDA shape plot

# BEGIN Writes samples statistics to data file
ci.mean.lev <- ci.var.lev <- ci.median.lev <- 0.95;
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Writing summary statistics" );
cat( paste( "# Observations: ", length( iat ), sep = "" ), "\n" );
#cat( paste( summary( iat ), "\n" );
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
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Plotting trend" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-trend.eps", sep = "" ),
	#title = "TERAGRID - Overall Interarrival Times",
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_trend( iat )
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-trend.eps", sep = "" ) );
	warning( "Error while drawing trend plot." );
	warning( "Error message is: ", geterrmessage() );
}
# END Plot trend

# BEGIN Plot frequency histogram
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Plotting frequency histogram" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fhist.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_freqHist( iat, title = NULL )
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fhist.eps", sep = "" ) );
	warning( "Error while drawing frequency histogram plot." );
	warning( "Error message is: ", geterrmessage() );
}
# END Plot frequency histogram

# BEGIN Plot density histogram
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Plotting density histogram" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-dhist.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_densHist( iat, title = NULL )
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-dhist.eps", sep = "" ) );
	warning( "Error while drawing density histogram plot." );
	warning( "Error message is: ", geterrmessage() );
}
# END Plot density histogram

# BEGIN Plot lag-1 plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Plotting Lag-1 plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lag1.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_plot_lag( iat, lag = 1 )
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lag1.eps", sep = "" ) );
	warning( "Error while drawing Lag-1 plot." );
	warning( "Error message is: ", geterrmessage() );
}
# END Plot lag-1 plot

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

##@{ Performs GoF tests

cat( "## Goodness-of-Fit tests:", "\n" );
fitDistrs <- names( fit );
for ( fitDistr in fitDistrs )
{
	cat( paste( "GoF tests for distribution: ", fitDistr, sep = "" ), "\n" );
	if ( length( fit[[ fitDistr ]] ) == 0 )
	{
		cat( "(WARNING) Fitting informations not available.", "\n" );
		next;
	}

	# BEGIN GoF tests
#BEGIN don't work
#	res <- try(
#		mg_gof_multigof.test(
#			iat,
#			dist = fitDistr,
#			tests = c(
#				MG_CONSTS_ANDERSONDARLING_GOF,
#				MG_CONSTS_KOLMOGOROVSMIRNOV_GOF,
#				MG_CONSTS_CHISQUARED_GOF
#			),
#			fit[[ fitDistr ]]
#		)
#	);
#END don't work
	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Doing Goodness-of-Fit tests for '", fitDistr, "' distribution" );
	res <- try(
		do.call(
			"mg_gof_multigof.test",
			c(
				list(iat),
				list(dist = fitDistr),
				fit[[ fitDistr ]]$estimate,
				list(
					tests = c(
#						MG_CONSTS_ANDERSONDARLING_GOF,
#						MG_CONSTS_ANDERSONDARLINGBOOTPAR_GOF,
#						MG_CONSTS_KOLMOGOROVSMIRNOV_GOF,
#						MG_CONSTS_KOLMOGOROVSMIRNOVBOOTPAR_GOF,
#						#MG_CONSTS_KOLMOGOROVSMIRNOVBOOTNONPAR_GOF,
						#MG_CONSTS_KOLMOGOROVSMIRNOVSIM_GOF,
						#MG_CONSTS_KOLMOGOROVSMIRNOV2_GOF,
						#MG_CONSTS_KOLMOGOROVSMIRNOV2BOOT_GOF,
						MG_CONSTS_CHISQUAREDPAR_GOF
#						MG_CONSTS_QQ_GOF
					),
					dist.nparams = length( fit[[ fitDistr ]]$estimate ),
					boot.n = 999,
					boot.fit.methods = fit[[ fitDistr ]]$methods
				)
			)
		)
	);
	if ( !is( res, "try-error" ) )
	{
		print( res );
	} else {
		warning( "Error while doing GoF tests for distribution '", fitDistr, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) GoF tests failed.", "\n" );
	}
	# END GoF tests

	# BEGIN Q-Q Plot
	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing Q-Q Plot for '", fitDistr, "' distribution" );
	cat( "## Q-Q Plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-qqplot.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		do.call(
			"mg_plot_qq",
			c(
				list(iat),
				list(dist = fitDistr),
				fit[[ fitDistr ]]$estimate
			)
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		cat( "\t", paste( "Slope: ", res$slope, " - Intercept: ", res$int ), "\n" );
		cat( "\t", "Correlation Coeff.: ", paste( "Pearson: ", res$cor.pearson$estimate, " @ [", res$cor.pearson$conf.int[1], ",", res$cor.pearson$conf.int[2], "] (p-value: ", res$cor.pearson$p.value, ")", sep = "" ), "\n" );
		cat( "\t", "Area Difference: Absolute: ", res$abs.area.diff, ", Relative: ", res$rel.area.diff, "\n" );
	} else {
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-qqplot.eps", sep = "" ) );
		warning( "Error while drawing Q-Q Plot distribution: '", fitDistr, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Q-Q Plot failed.", "\n" );
	}
	# END Q-Q Plot

#	# BEGIN Log-Log Q-Q Plot
#	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing Log-Log Q-Q Plot for '", fitDistr, "' distribution" );
#	cat( "## Log-Log Q-Q Plot:", "\n" );
#	postscript(
#		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-llqqplot.eps", sep = "" ),
#		horizontal = FALSE,
#		onefile = FALSE
#	);
#	res <- try(
#		do.call(
#			"mg_plot_qq",
#			c(
#				list(iat),
#				list(dist = fitDistr),
#				fit[[ fitDistr ]]$estimate,
#				log = "xy"
#			)
#		)
#	);
#	dev.off();
#	if ( !is( res, "try-error" ) )
#	{
#		cat( "\t", paste( "Slope: ", res$slope, " - Intercept: ", res$int ), "\n" );
#		cat( "\t", "Correlation Coeff.: ", paste( "Pearson: ", res$cor.pearson$estimate, " @ [", res$cor.pearson$conf.int[1], ",", res$cor.pearson$conf.int[2], "] (p-value: ", res$cor.pearson$p.value, ")", sep = "" ), "\n" );
#		cat( "\t", "Area Difference: Absolute: ", res$abs.area.diff, ", Relative: ", res$rel.area.diff, "\n" );
#	} else {
#		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-llqqplot.eps", sep = "" ) );
#		warning( "Error while drawing Log-Log Q-Q Plot distribution: '", fitDistr, "'." );
#		warning( "Error message is: ", geterrmessage() );
#		cat( "(ERROR) Log-Log Q-Q Plot failed.", "\n" );
#	}
#	# END Log-Log Q-Q Plot

	# BEGIN P-P Plot
	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing P-P Plot for '", fitDistr, "' distribution" );
	cat( "## P-P Plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-ppplot.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		do.call(
			"mg_plot_pp",
			c(
				list(iat),
				list(fitDistr),
				fit[[ fitDistr ]]$estimate
			)
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		cat( "\t", paste( "Slope: ", res$slope, " - Intercept: ", res$int ), "\n" );
		cat( "\t", "Correlation Coeff.: ", paste( "Pearson: ", res$cor.pearson$estimate, " @ [", res$cor.pearson$conf.int[1], ",", res$cor.pearson$conf.int[2], "] (p-value: ", res$cor.pearson$p.value, ")", sep = "" ), "\n" );
		cat( "\t", "Area Difference: Absolute: ", res$abs.area.diff, ", Relative: ", res$rel.area.diff, "\n" );
	} else {
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-ppplot.eps", sep = "" ) );
		warning( "Error while drawing P-P Plot distribution: '", fitDistr, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) P-P Plot failed.", "\n" );
	}
	# END P-P Plot
}

##@} Performs GoF tests

##@} Distribution fitting (parameters estimated from data set)

# BEGIN Plot and print ACF
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing ACF Plot" );
cat( "## ACF statistics and plot:", "\n" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-acf.eps", sep = "" ),
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
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-acf.eps", sep = "" ) );
	warning( "Error while drawing ACF plot." );
	warning( "Error message is: ", geterrmessage() );
	cat( "(ERROR) ACF failed.", "\n" );
}
# END Plot and print ACF

##@{ EDA autocorrelation plot
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Plotting EDA Autocorrelation Plot" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-eda-acor.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
	mg_eda.acor.plot( iat )
);
dev.off();
if ( is( res, "try-error" ) )
{
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-eda-acor.eps", sep = "" ) );
	warning( "Error while drawing EDA Autocorrelation plot." );
	warning( "Error message is: ", geterrmessage() );
}
##@} EDA autocorrelation plot

##@{ Plot and print Long Range Dependence indicators

mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing Long Range Dependence indicators Plot" );
cat( "## Long Range Dependence indicators statistics and plots:", "\n" );

##@{ Aggregated trend plot
cat( "## LRD Aggregated Trends plot:", "\n" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-aggtrend.eps", sep = "" ),
	horizontal = FALSE,
	onefile = FALSE
);
res <- try(
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
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-aggtrend.eps", sep = "" ) );
	warning( "Error while plotting Aggregated Trends." );
	warning( "Error message is: ", geterrmessage() );
	cat( "(ERROR) Aggregated Trends plot failed.", "\n" );
}
##@} Aggregated trend plot

##@{ Aggregated Log-Log CDF plot
cat( "## LRD Aggregated Log-Log Empirical CDF plot:", "\n" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-aggllecdf.eps", sep = "" ),
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
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-aggllecdf.eps", sep = "" ) );
	warning( "Error while plotting Aggregated Log-Log Empirical CDF." );
	warning( "Error message is: ", geterrmessage() );
	cat( "(ERROR) Aggregated Log-Log Empirical CDF plot failed.", "\n" );
}
##@} Aggregated Log-Log CDF plot

##@{ Aggregated Log-Log Empirical CCDF plot
cat( "## LRD Aggregated Log-Log Empirical CCDF plot:", "\n" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-agglleccdf.eps", sep = "" ),
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
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-agglleccdf.eps", sep = "" ) );
	warning( "Error while plotting Aggregated Log-Log CCDF." );
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
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-av.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_lrd.aggvar( res, plot.fitline = TRUE )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-av.eps", sep = "" ) );
		warning( "Error while plotting Aggregate of Variance method." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Aggregate of Variance method failed.", "\n" );
	}
} else {
	warning( "Error while performing Aggregate of Variance method." );
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
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-rs.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_lrd.rs( res, plot.fitline = TRUE )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-rs.eps", sep = "" ) );
		warning( "Error while plotting R/S method." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) R/S method failed.", "\n" );
	}
} else {
	warning( "Error while performing R/S method." );
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
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-pgram.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_lrd.periodogram( res, plot.fitline = TRUE )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-pgram.eps", sep = "" ) );
		warning( "Error while plotting Periodogram method." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Periodogram method failed.", "\n" );
	}
} else {
	warning( "Error while performing Periodogram method." );
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
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-cumpgram.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		mg_plot_lrd.periodogram( res, plot.fitline = TRUE )
	);
	dev.off();
	if ( is( res, "try-error" ) )
	{
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-lrd-cumpgram.eps", sep = "" ) );
		warning( "Error while plotting Cumulative Periodogram method." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Cumulative Periodogram method failed.", "\n" );
	}
} else {
	warning( "Error while performing Cumulative Periodogram method." );
	warning( "Error message is: ", geterrmessage() );
	cat( "(ERROR) Cumulative Periodogram method failed.", "\n" );
}
##@} Cumulative Periodogram method (for Hurst exponent)

##@} Plot and print Long Range Dependence indicators

##@{ Plot Hill estimators
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing Hill Plot" );
cat( "## Hill estimator:", "\n" );
postscript(
	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-hill.eps", sep = "" ),
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
	cat( paste( "Last Hill estimator: ",  res$statistic[ length( res$statistic ) ], " with (",  res$conf.int$lower[ length( res$conf.int$lower ) ], ", ", res$conf.int$upper[ length( res$conf.int$upper ) ], ") @ 95% level", sep = "" ), "\n" );
} else {
	unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-hill.eps", sep = "" ) );
	warning( "Error while drawing Hill plot." );
	warning( "Error message is: ", geterrmessage() ); 
	cat( "(ERROR) Hill estimator failed.", "\n" );
}
##@} Plot Hill estimators

sink( file = NULL );

options( error = olderror );
options( warn = oldwarn );
