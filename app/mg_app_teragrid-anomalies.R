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

## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

oldwarn <- getOption( "warn" );
options( warn = 1 );

olderror <- getOption( "error" );
options( error = quote( dump.frames( "mg_app_teragrid-anomalies", TRUE ) ) ); #DEBUG

source( "app/mg_app_conf.R" );
source( "lib/mg/mg_dists.R" );
source( "lib/mg/mg_fit.R" );
source( "lib/mg/mg_gof.R" );
source( "lib/mg/mg_teragrid.R" );
source( "lib/mg/mg_plot.R" );
source( "lib/mg/mg_statsutils.R" );

# no omit && no fix
tgFull <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME );

# omit{timestamp > startime}
tgNoAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE );
# omit{starttime == endtime}
tgNoSeqE <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE );
# omit{starttime > endtime}
tgNoSgtE <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSgtE = TRUE );
# omit{nodes <= 0}
tgNoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitNle0 = TRUE );
# fix{timestamp > starttime}
tgFixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, fixAgtS = TRUE );

# omit{timestamp > starttime} && omit{starttime == endtime}
tgNoAgtSNoSeqE <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitSeqE = TRUE );
# omit{timestamp > starttime} && omit{starttime > endtime}
tgNoAgtSNoSgtE <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitSgtE = TRUE );
# omit{timestamp > starttime} && omit{nodes <= 0}
tgNoAgtSNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitNle0 = TRUE );

# omit{starttime == endtime} && omit{starttime > endtime}
tgNoSeqENoSgtE <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE, omitSgtE = TRUE );
# omit{nodes == endtime} && omit{nodes <= 0}
tgNoSeqENoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE, omitNle0 = TRUE );
# omit{nodes == endtime} && fix{timestamp >  starttime}
tgNoSeqEFixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE, fixAgtS = TRUE );

# omit{starttime > endtime} && omit{nodes <= 0}
tgNoSgtENoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSgtE = TRUE, omitNle0 = TRUE );
# omit{starttime > endtime} && fix{timestamp > starttime}
tgNoSgtEFixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSgtE = TRUE, fixAgtS = TRUE );

# omit{nodes <= 0} && fix{timestamp > starttime}
tgNoNle0FixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitNle0 = TRUE, fixAgtS = TRUE );

# omit{timestamp > startime} && omit{starttime == endtime} && omit{starttime > endtime}
tgNoAgtSNoSeqENoSgtE <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitSeqE = TRUE, omitSgtE = TRUE );
# omit{timestamp > startime} && omit{starttime == endtime} && omit{nodes <= 0}
tgNoAgtSNoSeqENoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitSeqE = TRUE, omitNle0 = TRUE );
# omit{timestamp > startime} && omit{starttime > endtime} && omit{nodes <= 0}
tgNoAgtSNoSgtENoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitSgtE = TRUE, omitNle0 = TRUE );

# omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0}
tgNoSeqENoSgtENoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE, omitSgtE = TRUE, omitNle0 = TRUE );
# omit{starttime == endtime} && omit{starttime > endtime} && fix{timestamp gt starttime}
tgNoSeqENoSgtEFixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE, omitSgtE = TRUE, fixAgtS = TRUE );

# omit{starttime > endtime} && omit{nodes <= 0} && fix{timestamp gt starttime}
tgNoSgtENoNle0FixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSgtE = TRUE, omitNle0 = TRUE, fixAgtS = TRUE );

# omit{timestamp > starttime> && omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0}
tgNoAgtSNoSeqENoSgtENoNle0 <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitAgtS = TRUE , omitSeqE = TRUE , omitSgtE = TRUE, omitNle0 = TRUE );

# omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0} && fix{timestamp > starttime}
tgNoSeqENoSgtENoNle0FixAgtS <- mg_teragrid( MG_APP_CONF_TERAGRID_LOGFULLNAME, omitSeqE = TRUE , omitSgtE = TRUE, omitNle0 = TRUE, fixAgtS = TRUE );

# BEGIN print summary informations
colors <- rainbow(23);
tgMethods <- c(
	"iat" = "mg_teragrid_overallInterArrTimes",
	"rt" = "mg_teragrid_overallRuntimes"
);
for ( info in names( tgMethods ) )
{
	print( "" );
	print( paste( "### ", info, " ###", sep = "" ) );
	print( "" );

	m <- tgMethods[ info ];

	# no omit && no fix
	print( "TeraGrid - no omit && no fix" );
	dataFull <- do.call( m, list( tgFull ) );
	print( paste( "# Observations: ", length( dataFull ), sep = "" ) );
	print( summary( dataFull ) );
	print( paste( "CV: ", mg_statsutils_cv( dataFull ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataFull ), sep = "" ) );

	# omit{timestamp > startime}
	print( "TeraGrid - omit{timestamp > startime}" );
	dataNoAgtS <- do.call( m, list( tgNoAgtS ) );
	print( paste( "# Observations: ", length( dataNoAgtS ), sep = "" ) );
	print( summary( dataNoAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtS ), sep = "" ) );

	# omit{starttime == endtime}
	print( "TeraGrid - omit{starttime == endtime}" );
	dataNoSeqE <- do.call( m, list( tgNoSeqE ) );
	print( paste( "# Observations: ", length( dataNoSeqE ), sep = "" ) );
	print( summary( dataNoSeqE ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqE ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqE ), sep = "" ) );

	# omit{starttime > endtime}
	print( "TeraGrid - omit{starttime > endtime}" );
	dataNoSgtE <- do.call( m, list( tgNoSgtE ) );
	print( paste( "# Observations: ", length( dataNoSgtE ), sep = "" ) );
	print( summary( dataNoSgtE ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSgtE ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSgtE ), sep = "" ) );

	# omit{nodes <= 0}
	print( "TeraGrid - omit{nodes <= 0}" );
	dataNoNle0 <- do.call( m, list( tgNoNle0 ) );
	print( paste( "# Observations: ", length( dataNoNle0 ), sep = "" ) );
	print( summary( dataNoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoNle0 ), sep = "" ) );

	# fix{timestamp > starttime}
	print( "TeraGrid - fix{timestamp > starttime}" );
	dataFixAgtS <- do.call( m, list( tgFixAgtS ) );
	print( paste( "# Observations: ", length( dataFixAgtS ), sep = "" ) );
	print( summary( dataFixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataFixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataFixAgtS ), sep = "" ) );

	# omit{timestamp > starttime} && omit{starttime == endtime}
	print( "TeraGrid - omit{timestamp > starttime} && omit{starttime == endtime}" );
	dataNoAgtSNoSeqE <- do.call( m, list( tgNoAgtSNoSeqE ) );
	print( paste( "# Observations: ", length( dataNoAgtSNoSeqE ), sep = "" ) );
	print( summary( dataNoAgtSNoSeqE ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNoSeqE ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNoSeqE ), sep = "" ) );

	# omit{timestamp > starttime} && omit{starttime > endtime}
	print( "TeraGrid - omit{timestamp > starttime} && omit{starttime > endtime}" );
	dataNoAgtSNoSgtE <- do.call( m, list( tgNoAgtSNoSgtE ) );
	print( paste( "# Observations: ", length( dataNoAgtSNoSgtE ), sep = "" ) );
	print( summary( dataNoAgtSNoSgtE ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNoSgtE ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNoSgtE ), sep = "" ) );

	# omit{timestamp > starttime} && omit{nodes <= 0}
	print( "TeraGrid - omit{timestamp > starttime} && omit{nodes <= 0}" );
	dataNoAgtSNle0 <- do.call( m, list( tgNoAgtSNle0 ) );
	print( paste( "# Observations: ", length( dataNoAgtSNle0 ), sep = "" ) );
	print( summary( dataNoAgtSNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNle0 ), sep = "" ) );

	# omit{starttime == endtime} && omit{starttime > endtime}
	print( "TeraGrid - omit{starttime == endtime} && omit{starttime > endtime}" );
	dataNoSeqENoSgtE <- do.call( m, list( tgNoSeqENoSgtE ) );
	print( paste( "# Observations: ", length( dataNoSeqENoSgtE ), sep = "" ) );
	print( summary( dataNoSeqENoSgtE ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqENoSgtE ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqENoSgtE ), sep = "" ) );

	# omit{nodes == endtime} && omit{nodes <= 0}
	print( "TeraGrid - omit{nodes == endtime} && omit{nodes <= 0}" );
	dataNoSeqENoNle0 <- do.call( m, list( tgNoSeqENoNle0 ) );
	print( paste( "# Observations: ", length( dataNoSeqENoNle0 ), sep = "" ) );
	print( summary( dataNoSeqENoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqENoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqENoNle0 ), sep = "" ) );

	# omit{nodes == endtime} && fix{timestamp >  starttime}
	print( "TeraGrid - omit{nodes == endtime} && fix{timestamp >  starttime}" );
	dataNoSeqEFixAgtS <- do.call( m, list( tgNoSeqEFixAgtS ) );
	print( paste( "# Observations: ", length( dataNoSeqEFixAgtS ), sep = "" ) );
	print( summary( dataNoSeqEFixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqEFixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqEFixAgtS ), sep = "" ) );

	# omit{starttime > endtime} && omit{nodes <= 0}
	print( "TeraGrid - omit{starttime > endtime} && omit{nodes <= 0}" );
	dataNoSgtENoNle0 <- do.call( m, list( tgNoSgtENoNle0 ) );
	print( paste( "# Observations: ", length( dataNoSgtENoNle0 ), sep = "" ) );
	print( summary( dataNoSgtENoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSgtENoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSgtENoNle0 ), sep = "" ) );

	# omit{starttime > endtime} && fix{timestamp > starttime}
	print( "TeraGrid - omit{starttime > endtime} && fix{timestamp > starttime}" );
	dataNoSgtEFixAgtS <- do.call( m, list( tgNoSgtEFixAgtS ) );
	print( paste( "# Observations: ", length( dataNoSgtEFixAgtS ), sep = "" ) );
	print( summary( dataNoSgtEFixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSgtEFixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSgtEFixAgtS ), sep = "" ) );

	# omit{nodes <= 0} && fix{timestamp > starttime}
	print( "TeraGrid - omit{nodes <= 0} && fix{timestamp > starttime}" );
	dataNoNle0FixAgtS <- do.call( m, list( tgNoNle0FixAgtS ) );
	print( paste( "# Observations: ", length( dataNoNle0FixAgtS ), sep = "" ) );
	print( summary( dataNoNle0FixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoNle0FixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoNle0FixAgtS ), sep = "" ) );

	# omit{timestamp > startime} && omit{starttime == endtime} && omit{starttime > endtime}
	print( "TeraGrid - omit{timestamp > startime} && omit{starttime == endtime} && omit{starttime > endtime}" );
	dataNoAgtSNoSeqENoSgtE <- do.call( m, list( tgNoAgtSNoSeqENoSgtE ) );
	print( paste( "# Observations: ", length( dataNoAgtSNoSeqENoSgtE ), sep = "" ) );
	print( summary( dataNoAgtSNoSeqENoSgtE ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNoSeqENoSgtE ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNoSeqENoSgtE ), sep = "" ) );

	# omit{timestamp > startime} && omit{starttime == endtime} && omit{nodes <= 0}
	print( "TeraGrid - omit{timestamp > startime} && omit{starttime == endtime} && omit{nodes <= 0}" );
	dataNoAgtSNoSeqENoNle0 <- do.call( m, list( tgNoAgtSNoSeqENoNle0 ) );
	print( paste( "# Observations: ", length( dataNoAgtSNoSeqENoNle0 ), sep = "" ) );
	print( summary( dataNoAgtSNoSeqENoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNoSeqENoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNoSeqENoNle0 ), sep = "" ) );

	# omit{timestamp > startime} && omit{starttime > endtime} && omit{nodes <= 0}
	print( "TeraGrid - omit{timestamp > startime} && omit{starttime > endtime} && omit{nodes <= 0}" );
	dataNoAgtSNoSgtENoNle0 <- do.call( m, list( tgNoAgtSNoSgtENoNle0 ) );
	print( paste( "# Observations: ", length( dataNoAgtSNoSgtENoNle0 ), sep = "" ) );
	print( summary( dataNoAgtSNoSgtENoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNoSgtENoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNoSgtENoNle0 ), sep = "" ) );

	# omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0}
	print( "TeraGrid - omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0}" );
	dataNoSeqENoSgtENoNle0 <- do.call( m, list( tgNoSeqENoSgtENoNle0 ) );
	print( paste( "# Observations: ", length( dataNoSeqENoSgtENoNle0 ), sep = "" ) );
	print( summary( dataNoSeqENoSgtENoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqENoSgtENoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqENoSgtENoNle0 ), sep = "" ) );

	# omit{starttime == endtime} && omit{starttime > endtime} && fix{timestamp gt starttime}
	print( "TeraGrid - omit{starttime == endtime} && omit{starttime > endtime} && fix{timestamp gt starttime}" );
	dataNoSeqENoSgtEFixAgtS <- do.call( m, list( tgNoSeqENoSgtEFixAgtS ) );
	print( paste( "# Observations: ", length( dataNoSeqENoSgtEFixAgtS ), sep = "" ) );
	print( summary( dataNoSeqENoSgtEFixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqENoSgtEFixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqENoSgtEFixAgtS ), sep = "" ) );

	# omit{starttime > endtime} && omit{nodes <= 0} && fix{timestamp gt starttime}
	print( "TeraGrid - omit{starttime > endtime} && omit{nodes <= 0} && fix{timestamp gt starttime}" );
	dataNoSgtENoNle0FixAgtS <- do.call( m, list( tgNoSgtENoNle0FixAgtS ) );
	print( paste( "# Observations: ", length( dataNoSgtENoNle0FixAgtS ), sep = "" ) );
	print( summary( dataNoSgtENoNle0FixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSgtENoNle0FixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSgtENoNle0FixAgtS ), sep = "" ) );

	# omit{timestamp > starttime} && omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0}
	print( "TeraGrid - omit{timestamp > starttime} && omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0}" );
	dataNoAgtSNoSeqENoSgtENoNle0 <- do.call( m, list( tgNoAgtSNoSeqENoSgtENoNle0 ) );
	print( paste( "# Observations: ", length( dataNoAgtSNoSeqENoSgtENoNle0 ), sep = "" ) );
	print( summary( dataNoAgtSNoSeqENoSgtENoNle0 ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoAgtSNoSeqENoSgtENoNle0 ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoAgtSNoSeqENoSgtENoNle0 ), sep = "" ) );

	# omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0} && fix{timestamp > starttime}
	print( "TeraGrid - omit{starttime == endtime} && omit{starttime > endtime} && omit{nodes <= 0} && fix{timestamp > starttime}" );
	dataNoSeqENoSgtENoNle0FixAgtS <- do.call( m, list( tgNoSeqENoSgtENoNle0FixAgtS ) );
	print( paste( "# Observations: ", length( dataNoSeqENoSgtENoNle0FixAgtS ), sep = "" ) );
	print( summary( dataNoSeqENoSgtENoNle0FixAgtS ) );
	print( paste( "CV: ", mg_statsutils_cv( dataNoSeqENoSgtENoNle0FixAgtS ), sep = "" ) );
	print( paste( "Std Err: ", sd( dataNoSeqENoSgtENoNle0FixAgtS ), sep = "" ) );

	postscript(
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-", info, "-boxplot.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	mg_plot_boxplot(
		dataFull,
#		dataNoAgtS,
#		dataNoSeqE,
#		dataNoSgtE,
#		dataNoNle0,
#		dataFixAgtS,
#		dataNoAgtSNoSeqE,
#		dataNoAgtSNoSgtE,
#		dataNoAgtSNle0,
#		dataNoSeqENoSgtE,
#		dataNoSeqENoNle0,
#		dataNoSeqEFixAgtS,
#		dataNoSgtENoNle0,
#		dataNoSgtEFixAgtS,
#		dataNoNle0FixAgtS,
#		dataNoAgtSNoSeqENoSgtE,
#		dataNoAgtSNoSeqENoNle0,
#		dataNoAgtSNoSgtENoNle0,
#		dataNoSeqENoSgtENoNle0,
#		dataNoSeqENoSgtEFixAgtS,
#		dataNoSgtENoNle0FixAgtS,
#		dataNoAgtSNoSeqENoSgtENoNle0,
		dataNoSeqENoSgtENoNle0FixAgtS,
		#border = colors,
		names = c(
			"full",
#			"! A>S",
#			"! S=E",
#			"! S>E",
#			"! N<=0",
#			"Fix A>S",
#			"(! A>S) & (! S=E)",
#			"(! A>S) & (! S>E)",
#			"(! A>S) & (! N<=0)",
#			"(! S=E) & (! S>E)",
#			"(! S=E) & (! N<=0)",
#			"(! S=E) & (Fix A>S=",
#			"(! S>E) & (! N<=0)",
#			"(! S>E) & (Fix A>S)",
#			"(! N<=0) & (Fix A>S)",
#			"(! A>S) & (! S=E) & (! S>E)",
#			"(! A>S) & (! S=E) & (! N<=0)",
#			"(! A>S) & (! S>E) & (! N<=0)",
#			"(! S=E) & (! S>E) & (! N<=0)",
#			"(! S=E) & (! S>E) & (Fix A>S)",
#			"(! S>E) & (! N<=0) & (Fix A>S)",
#			"(! A>S) & (! S=E) & (! S>E) & (! N<=0)",
#			"(! S=E) & (! S>E) & (! N<=0) & (Fix A>S)"
			"all fixed"
		)
	);
#	legend(
#		"topright",
#		legend = c(
#			"full",
#			"! A>S",
#			"! S=E",
#			"! S>E",
#			"! N<=0",
#			"Fix A>S",
#			"(! A>S) & (! S=E)",
#			"(! A>S) & (! S>E)",
#			"(! A>S) & (! N<=0)",
#			"(! S=E) & (! S>E)",
#			"(! S=E) & (! N<=0)",
#			"(! S=E) & (Fix A>S=",
#			"(! S>E) & (! N<=0)",
#			"(! S>E) & (Fix A>S)",
#			"(! N<=0) & (Fix A>S)",
#			"(! A>S) & (! S=E) & (! S>E)",
#			"(! A>S) & (! S=E) & (! N<=0)",
#			"(! A>S) & (! S>E) & (! N<=0)",
#			"(! S=E) & (! S>E) & (! N<=0)",
#			"(! S=E) & (! S>E) & (Fix A>S)",
#			"(! S>E) & (! N<=0) & (Fix A>S)",
#			"(! A>S) & (! S=E) & (! S>E) & (! N<=0)",
#			"(! S=E) & (! S>E) & (! N<=0) & (Fix A>S)"
#		),
#		col = colors
#	);
#	grid();
	dev.off();
}
# END print summary informations

options( warn = oldwarn );
options( error = olderror );
