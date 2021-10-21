#ifndef __LAPACK_H
#define __LAPACK_H

typedef long int integer;
typedef double doublereal;

#if defined __cplusplus
extern "C" {
#endif

int dspev_(const char *jobz, const char *uplo, const integer *n, doublereal *ap,
	doublereal *w, doublereal *z, const integer *ldz, doublereal *work,
	integer *info);
int dgesv_(const integer *n, const integer *nrhs, doublereal *a, const integer *lda,
	integer *ipiv, doublereal *b, const integer *ldb, integer *info);

int dsyev_(const char *jobz, const char *uplo, const integer *n, doublereal *a, const integer *lda, 
	doublereal *w, doublereal *work, const integer *lwork, integer *info);

int dspsv_(const char *uplo, const integer *n, const integer *nrhs, doublereal *ap,
	integer *ipiv, doublereal *b, const integer *ldb, integer *info);

int dppsv_(const char *uplo, const integer *n, const integer *nrhs, doublereal *ap,
	doublereal *b, const integer *ldb, integer *info);

int dgetrf_(const integer *m, const integer *n, doublereal *a, const integer * lda,
	integer *ipiv, integer *info);

int dgetri_(const integer *n, doublereal *a, const integer *lda, const integer *ipiv,
	doublereal *work, const integer *lwork, integer *info);

int dgetrs_(const char *trans, const integer *n, const integer *nrhs, const doublereal *a, 
	const integer *lda, const integer *ipiv, doublereal *b, const integer * ldb,
	integer *info);

int dpptrf_(const char *uplo, const integer *n, doublereal *ap, integer *info);

int dpptri_(const char *uplo, const integer *n, doublereal *ap, integer *info);

int dpptrs_(const char *uplo, const integer *n, const integer *nrhs, const doublereal *ap,
	doublereal *b, const integer *ldb, integer *info);

int dsptrf_(const char *uplo, const integer *n, doublereal *ap, integer *ipiv, integer *info);

int dsptri_(const char *uplo, const integer *n, doublereal *ap, const integer *ipiv,
	doublereal *work, integer *info);

int dsptrs_(const char *uplo, const integer *n, const integer *nrhs, const doublereal *ap,
	const integer *ipiv, doublereal *b, const integer *ldb, integer *info);

#if defined __cplusplus
}
#endif

#endif /* __LAPACK_H */

