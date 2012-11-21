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

#include "mg_matrix.h"
/*
#include "mg_utils.h"
#include <R_ext/Lapack.h>
*/

const static double padec [] = /* for matrix exponential calculation. */
{
  5.0000000000000000e-1,
  1.1666666666666667e-1,
  1.6666666666666667e-2,
  1.6025641025641026e-3,
  1.0683760683760684e-4,
  4.8562548562548563e-6,
  1.3875013875013875e-7,
  1.9270852604185938e-9,
};

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])

/**
 * Matrix exponential - based on the code for Octave's expm function.
 *
 * @param x real square matrix to exponentiate
 *
 * @return matrix exponential of x
 */
SEXP mg_matrix_exp(SEXP m)
{
	SEXP val = PROTECT(duplicate(m));
	//int *Dims = INTEGER(GET_SLOT(m, Matrix_DimSym));/*FIXME*/
	//int *Dims = INTEGER(getAttrib(m, R_DimSymbol));
	//SEXP Dims = PROTECT(allocVector(INTSXP,2));
	int i;
	int ilo;
	int ilos;
	int ihi;
	int ihis;
	int j;
	int sqpow;
	int nrow = INTEGER(GET_DIM(m))[0];
	int ncol = INTEGER(GET_DIM(m))[1];
	int ncp1 = ncol + 1, ncsqr = ncol * ncol;
	//int *pivot = Calloc(ncol, int);
	int *pivot = Alloca(ncol, int);
	//int *iperm = Calloc(ncol, int);
	int *iperm = Alloca(ncol, int);
	//double *dpp = Calloc(ncsqr, double); /* denominator power Pade' */
	double *dpp = Alloca(ncsqr, double); /* denominator power Pade' */
	//double *npp = Calloc(ncsqr, double); /* numerator power Pade' */
	double *npp = Alloca(ncsqr, double); /* numerator power Pade' */
	//double *perm = Calloc(ncol, double);
	double *perm = Alloca(ncol, double);
	//double *scale = Calloc(ncol, double);
	double *scale = Alloca(ncol, double);
	//double *v = REAL(GET_SLOT(val, Matrix_xSym));/*FIXME*/
	double *v = REAL(val);
	//double *work = Calloc(ncsqr, double);
	double *work = Alloca(ncsqr, double);
	double inf_norm;
	double m1_j; /* (-1)^j */
	double one = 1.0;
	double trshift;
	double zero = 0.0;

	R_CheckStack();

	if (ncol < 1 || nrow != ncol)
		error(_("Matrix exponential requires square, non-null matrix"));

	/* FIXME: Add special treatment for ncol == 1 */

	/* Preconditioning 1.  Shift diagonal by average diagonal if positive. */
	trshift = 0;		/* determine average diagonal element */
	for (i = 0; i < ncol; i++)
	{
		trshift += v[i * ncp1];
	}
	trshift /= ncol;
	if (trshift > 0.)
	{
		/* shift diagonal by -trshift */
		for (i = 0; i < ncol; i++)
		{
			v[i * ncp1] -= trshift;
		}
	}

	/* Preconditioning 2. Balancing with dgebal. */
	F77_CALL(dgebal)( "P", &ncol, v, &ncol, &ilo, &ihi, perm, &j );
	if (j)
	{
		error(_("mg_matrix_exp: LAPACK routine dgebal returned %d"), j);
	}
	F77_CALL(dgebal)( "S", &ncol, v, &ncol, &ilos, &ihis, scale, &j );
	if (j)
	{
		error(_("mg_matrix_exp: LAPACK routine dgebal returned %d"), j);
	}

	/* Preconditioning 3. Scaling according to infinity norm */
	inf_norm = F77_CALL(dlange)( "I", &ncol, &ncol, v, &ncol, work );
	sqpow = (inf_norm > 0) ? (int) (1 + log(inf_norm)/log(2.)) : 0;
	if (sqpow < 0)
	{
		sqpow = 0;
	}
	if (sqpow > 0)
	{
		double scale_factor = 1.0;
		for (i = 0; i < sqpow; i++)
		{
			scale_factor *= 2.;
		}
		for (i = 0; i < ncsqr; i++)
		{
			v[i] /= scale_factor;
		}
	}

	/* Pade' approximation. Powers v^8, v^7, ..., v^1 */
	AZERO(npp, ncsqr);
	AZERO(dpp, ncsqr);
	m1_j = -1;
	for (j = 7; j >=0; j--)
	{
		double mult = padec[j];
		/* npp = m * npp + padec[j] *m */
		F77_CALL(dgemm)("N", "N", &ncol, &ncol, &ncol, &one, v, &ncol, npp, &ncol,
			&zero, work, &ncol);
		for (i = 0; i < ncsqr; i++) npp[i] = work[i] + mult * v[i];
		/* dpp = m * dpp * (m1_j * padec[j]) * m */
		mult *= m1_j;
		F77_CALL(dgemm)("N", "N", &ncol, &ncol, &ncol, &one, v, &ncol, dpp, &ncol,
			&zero, work, &ncol);
		for (i = 0; i < ncsqr; i++) dpp[i] = work[i] + mult * v[i];
		m1_j *= -1;
	}
	/* Zero power */
	for (i = 0; i < ncsqr; i++)
	{
		dpp[i] *= -1.;
	}
	for (j = 0; j < ncol; j++)
	{
		npp[j * ncp1] += 1.;
		dpp[j * ncp1] += 1.;
	}

	/* Pade' approximation is solve(dpp, npp) */
	F77_CALL(dgetrf)(&ncol, &ncol, dpp, &ncol, pivot, &j);
	if (j)
	{
		error(_("mg_matrix_exp: dgetrf returned error code %d"), j);
	}
	F77_CALL(dgetrs)("N", &ncol, &ncol, dpp, &ncol, pivot, npp, &ncol, &j);
	if (j)
	{
		error(_("mg_matrix_exp: dgetrs returned error code %d"), j);
	}
	Memcpy(v, npp, ncsqr);

	/* Now undo all of the preconditioning */
	/* Preconditioning 3: square the result for every power of 2 */
	while (sqpow--)
	{
		F77_CALL(dgemm)("N", "N", &ncol, &ncol, &ncol, &one, v, &ncol, v, &ncol, &zero, work, &ncol);
		Memcpy(v, work, ncsqr);
	}
	/* Preconditioning 2: apply inverse scaling */
	for (j = 0; j < ncol; j++)
	{
		for (i = 0; i < ncol; i++)
		{
			v[i + j * ncol] *= scale[i]/scale[j];
		}
	}
	/* Construct balancing permutation vector */
	for (i = 0; i < ncol; i++)
	{
		iperm[i] = i; /* identity permutation */
	}
	/* Leading permutations applied in forward order */
	for (i = 0; i < (ilo - 1); i++)
	{
		int swapidx = (int) (perm[i]) - 1;
		int tmp = iperm[i];
		iperm[i] = iperm[swapidx];
		iperm[swapidx] = tmp;
	}
	/* Trailing permutations applied in reverse order */
	for (i = ncol - 1; i >= ihi; i--)
	{
		int swapidx = (int) (perm[i]) - 1;
		int tmp = iperm[i];
		iperm[i] = iperm[swapidx];
		iperm[swapidx] = tmp;
	}
	/* Construct inverse balancing permutation vector */
	Memcpy(pivot, iperm, ncol);
	for (i = 0; i < ncol; i++)
	{
		iperm[pivot[i]] = i;
	}
	/* Apply inverse permutation */
	Memcpy(work, v, ncsqr);
	for (j = 0; j < ncol; j++)
	{
		for (i = 0; i < ncol; i++)
		{
			v[i + j * ncol] = work[iperm[i] + iperm[j] * ncol];
		}
	}

	/* Preconditioning 1: Trace normalization */
	if (trshift > 0.)
	{
		double mult = exp(trshift);
		for (i = 0; i < ncsqr; i++)
		{
			v[i] *= mult;
		}
	}

	/* Clean up */
	//Free(dpp); Free(npp); Free(perm); Free(iperm); Free(pivot); Free(scale); Free(work);

	UNPROTECT(1);
//	SEXP ans;
//    SEXP Dims = PROTECT(allocVector(INTSXP,2));
//double *v = REAL(val);
//    INTEGER(Dims)[0] = 1;
//    INTEGER(Dims)[1] = 10;
//fprintf(stderr, "Ciao %d\n", INTEGER(Dims)[0] );
//fprintf(stderr, "Ciao %e\n", v[0] );
////SEXP ans =  PROTECT(allocMatrix(REALSXP, 2, 2));
// //   ans = Memcpy(REAL(ans), REAL(m), 2 * 2);
//
//    UNPROTECT(1);
//    UNPROTECT(1);
//    return val;
//	int nrow=INTEGER(GET_DIM(m))[0];
//	int ncol=INTEGER(GET_DIM(m))[1];
//	fprintf(stderr, "Rows %d\n", nrow );
//	fprintf(stderr, "Cols %d\n", ncol );
	//fprintf(stderr, "First %g\n", RMATRIX(m, 0, 1) );
	return( val );
}
