/*
 * Copyright (C) 2008       Marco Guazzone
 *                          [Distributed Computing System (DCS) Group,
 *                           Computer Science Institute,
 *                           Department of Science and Technological Innovation,
 *                           University of Piemonte Orientale,
 *                           Alessandria (Italy)]
 *
 * This file is part of dcsr-workload-characterization.
 *
 * dcsr-workload-characterization is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * dcsr-workload-characterization is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with dcsr-workload-characterization. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MG_MATRIX_H_

#define _MG_MATRIX_H_

#include <R_ext/Lapack.h>
#include "mg_utils.h"

SEXP mg_matrix_exp(SEXP m);

void F77_NAME(dgesdd)(const char *jobz,
                      const int *m, const int *n,
                      double *a, const int *lda, double *s,
                      double *u, const int *ldu,
                      double *vt, const int *ldvt,
                      double *work, const int *lwork, int *iwork, int *info);

#endif /* _MG_MATRIX_H_ */
