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

## MG_GOF
##
## SUMMARY
##  Goodness-of-fit tests.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_consts.R" );
source( "lib/mg/mg_gof_ad.R" );
source( "lib/mg/mg_gof_chisq.R" );
source( "lib/mg/mg_gof_ks.R" );
source( "lib/mg/mg_gof_qq.R" );

mg_gof_multigof.test <- function( x, dist, ..., tests, dist.nparams = NULL, boot.n = 999, boot.fit.methods = NULL, boot.rep = FALSE, sim.n = 1000, bin.method = c("moore", "sturges"), alternative = c("two.sided", "less", "greater"), exact = NULL )
{
	if ( !is.character( dist ) )
	{
		stop( "[MG::GOF::MULTIGOF::test] Expected distribution name." );
	}

	res <- list();

	for ( t in tests )
	{
		if ( t == MG_CONSTS_KOLMOGOROVSMIRNOV_GOF )
		{
			res$ks <- try( mg_gof_ks.test( x, dist, ..., alternative = alternative, exact = exact ) );
		} else if ( t == MG_CONSTS_KOLMOGOROVSMIRNOVBOOTPAR_GOF ) {
			res$ks.boot.par <- try( mg_gof_ks.test.boot.par( x, dist, ..., n.boot = boot.n, fit.methods = boot.fit.methods, alternative = alternative, exact = exact ) );
		} else if ( t == MG_CONSTS_KOLMOGOROVSMIRNOVBOOTNONPAR_GOF ) {
			res$ks.boot.nonpar <- try( mg_gof_ks.test.boot.nonpar( x, dist, ..., n.boot = boot.n, alternative = alternative, exact = exact ) );
		} else if ( t == MG_CONSTS_KOLMOGOROVSMIRNOVSIM_GOF ) {
			res$ks.sim <- try( mg_gof_ks.test.sim( x, dist, ..., n.sim = sim.n, alternative = alternative, exact = exact ) );
		} else if ( t == MG_CONSTS_KOLMOGOROVSMIRNOV2_GOF ) {
			# Generate random sample from 'dist'
			y <- mg_dists_rvg( dist, length(x), ... );

			res$ks2 <- try( mg_gof_ks.test2( x, y, alternative = alternative, exact = exact ) );
		} else if ( t == MG_CONSTS_KOLMOGOROVSMIRNOV2BOOT_GOF ) {
			# Generate random sample from 'dist'
			y <- mg_dists_rvg( dist, length(x), ... );

			res$ks2.boot <- try( mg_gof_ks.test2.boot( x, y, n.boot = boot.n, replace = boot.rep, alternative = alternative, exact = exact ) );
		} else if ( t == MG_CONSTS_ANDERSONDARLING_GOF ) {
			res$ad <- try( mg_gof_ad.test( x, dist, ... ) );
		} else if ( t == MG_CONSTS_ANDERSONDARLINGBOOTPAR_GOF ) {
			res$ad.boot.par <- try( mg_gof_ad.test.boot.par( x, dist, ..., n.boot = boot.n, fit.methods = boot.fit.methods, alternative = alternative ) );
		} else if ( t == MG_CONSTS_CHISQUAREDPAR_GOF ) {
			# Parametric Chi-Square test
			res$chisq.par <- try( mg_gof_chisq.test.par( x, dist, ..., n.params = ifelse( is.null(dist.nparams), length(list(...)), dist.nparams ), bin.method = bin.method ) );
		} else if ( t == MG_CONSTS_CHISQUAREDNONPAR_GOF ) {
			# Non-Parametric Chi-Square test
			res$chisq.nonpar <- try( mg_gof_chisq.test.nonpar( x, dist, ..., bin.method = bin.method ) );
		} else if ( t == MG_CONSTS_QQ_GOF ) {
			res$qq <- try( mg_gof_qq.test( x, dist, ... ) );
		}
	}
	return( res );
}
