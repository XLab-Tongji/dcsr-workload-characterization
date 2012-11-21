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

#ifndef _MG_UTILS_H_
#define _MG_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <Rdefines.h> /* Rinternals.h + GET_SLOT etc */
#include <R.h>  /* includes Rconfig.h */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("mg_atrix", String)
#else
#define _(String) (String)
#endif

#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/* enum constants from cblas.h and some short forms */
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
#define RMJ CblasRowMajor
#define CMJ CblasColMajor
#define NTR CblasNoTrans
#define TRN CblasTrans
#define CTR CblasConjTrans
#define UPP CblasUpper
#define LOW CblasLower
#define NUN CblasNonUnit
#define UNT CblasUnit
#define LFT CblasLeft
#define RGT CblasRight

#define PACKED_TO_FULL(TYPE)                                            \
TYPE *packed_to_full_ ## TYPE(TYPE *dest, const TYPE *src,              \
                             int n, enum CBLAS_UPLO uplo)
PACKED_TO_FULL(double);
PACKED_TO_FULL(int);
#undef PACKED_TO_FULL

#define FULL_TO_PACKED(TYPE)                                            \
TYPE *full_to_packed_ ## TYPE(TYPE *dest, const TYPE *src, int n,       \
                              enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag)
FULL_TO_PACKED(double);
FULL_TO_PACKED(int);
#undef FULL_TO_PACKED

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

#ifdef __cplusplus
}
#endif

#endif /* _MG_UTILS_H_ */
