/* Experimental (dev/drop-in-replacement.md, B6): stand-in for libopenblas64_ in the
   juliac bundle. Packmol makes no BLAS calls, but LinearAlgebra loads the library at
   startup, and libblastrampoline probes these symbols to detect its interface.
   Build: gcc -O2 -shared -fPIC -o libopenblas64_.so blas_stub.c -lm */
#include <stdint.h>
#include <math.h>
int64_t isamax_64_(const int64_t *n, const float *x, const int64_t *incx) {
    if (*n <= 0 || *incx <= 0) return 0;
    int64_t imax = 1; float vmax = fabsf(x[0]);
    for (int64_t i = 1; i < *n; i++) { float v = fabsf(x[i * *incx]); if (v > vmax) { vmax = v; imax = i + 1; } }
    return imax;
}
#include <complex.h>
double complex zdotc_64_(const int64_t *n, const double complex *x, const int64_t *incx,
                         const double complex *y, const int64_t *incy) {
    double complex s = 0;
    for (int64_t i = 0; i < *n; i++) s += conj(x[i * *incx]) * y[i * *incy];
    return s;
}
float complex cdotc_64_(const int64_t *n, const float complex *x, const int64_t *incx,
                        const float complex *y, const int64_t *incy) {
    float complex s = 0;
    for (int64_t i = 0; i < *n; i++) s += conjf(x[i * *incx]) * y[i * *incy];
    return s;
}
float sdot_64_(const int64_t *n, const float *x, const int64_t *incx, const float *y, const int64_t *incy) {
    float s = 0;
    for (int64_t i = 0; i < *n; i++) s += x[i * *incx] * y[i * *incy];
    return s;
}
