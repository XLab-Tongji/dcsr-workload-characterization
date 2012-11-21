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

# Distributions Identifiers
MG_CONSTS_BETA_DIST <- "beta";
MG_CONSTS_BINOMIAL_DIST <- "binomial";
MG_CONSTS_CAUCHY_DIST <- "cauchy";
MG_CONSTS_CHISQUARED_DIST <- "chisq";
MG_CONSTS_CONTINUOUSPHASETYPE_DIST <- "cph";
MG_CONSTS_EXPONENTIAL_DIST <- "exp";
MG_CONSTS_FISHERF_DIST <- "f";
MG_CONSTS_FRECHET_DIST <- "frechet";
MG_CONSTS_GAMMA_DIST <- "gamma";
MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST <- "gev";
MG_CONSTS_GENERALIZEDPARETO_DIST <- "gpd";
MG_CONSTS_GUMBEL_DIST <- "gumbel";
MG_CONSTS_HYPEREXPONENTIAL_DIST <- "hyperexp";
MG_CONSTS_LOGISTIC_DIST <- "logistic";
MG_CONSTS_LOGNORMAL_DIST <- "lognormal";
MG_CONSTS_MMPP_DIST <- "mmpp";
MG_CONSTS_NORMAL_DIST <- "normal";
MG_CONSTS_PARETO_DIST <- "pareto";
MG_CONSTS_POISSON_DIST <- "poisson";
MG_CONSTS_REVERSEWEIBULL_DIST <- "rweibull";
MG_CONSTS_SKEWEDSTABLE_DIST <- "skewedstable";
MG_CONSTS_STABLE_DIST <- "stable";
MG_CONSTS_STUDENTT_DIST <- "t";
MG_CONSTS_SYMSTABLE_DIST <- "symstable";
MG_CONSTS_UNIFORM_DIST <- "uniform";
MG_CONSTS_WEIBULL_DIST <- "weibull";
MG_CONSTS_WEIBULL3_DIST <- "weibull3";

# Distribution Fitting Identifiers
MG_CONSTS_EM_FIT <- "em";
MG_CONSTS_HILL_FIT <- "hill";
MG_CONSTS_LEASTSQUARED_FIT <- "lsq";
MG_CONSTS_MLE_FIT <- "mle";
#MG_CONSTS_MOMENTMATCHING_FIT <- "mmatch"; # see MG_CONSTS_METHODOFMOMENTS_FIT
MG_CONSTS_METHODOFMOMENTS_FIT <- "mom";
MG_CONSTS_QUANTILEMETHOD_FIT <- "quantile";
MG_CONSTS_YULE_FIT <- "yule";

# Goodness-of-Fit Identifiers
MG_CONSTS_ANDERSONDARLING_GOF <- "ad"; # Anderson-Darling GoF test
MG_CONSTS_ANDERSONDARLINGBOOTPAR_GOF <- "ad.boot.par"; # Anderson-Darling GoF test with parametric bootstraping
MG_CONSTS_CHISQUAREDPAR_GOF <- "chisq.par"; # Chi-Square GoF parametric test
MG_CONSTS_CHISQUAREDNONPAR_GOF <- "chisq.nonpar"; # Chi-Square GoF non-parametric test
MG_CONSTS_KOLMOGOROVSMIRNOV_GOF <- "ks"; # Kolmogorov-Smirnov 1-sample GoF test
MG_CONSTS_KOLMOGOROVSMIRNOVBOOTPAR_GOF <- "ks.boot.par"; # Kolmogorov-Smirnov 1-sample GoF test with parametric bootstraping
MG_CONSTS_KOLMOGOROVSMIRNOVBOOTNONPAR_GOF <- "ks.boot.nonpar"; # Kolmogorov-Smirnov 1-sample GoF test with non-parametric bootstraping
MG_CONSTS_KOLMOGOROVSMIRNOVSIM_GOF <- "ks.sim"; # Kolmogorov-Smirnov 1-sample GoF test with simulation
MG_CONSTS_KOLMOGOROVSMIRNOV2_GOF <- "ks2"; # Kolmogorov-Smirnov 2-sample GoF test
MG_CONSTS_KOLMOGOROVSMIRNOV2BOOT_GOF <- "ks2.boot"; # Kolmogorov-Smirnov 2-sample GoF test with bootstraping
MG_CONSTS_QQ_GOF <- "qq"; # QQ-plot based GoF test

# Plot Type Identifiers
MG_CONSTS_CDF_PLOT <- "cdf";
MG_CONSTS_LLCDF_PLOT <- "llcdf";
MG_CONSTS_CCDF_PLOT <- "ccdf";
MG_CONSTS_LLCD_PLOT <- "llcd";
MG_CONSTS_PDF_PLOT <- "pdf";
#MG_CONSTS_LLPDF_PLOT <- "llpdf";
