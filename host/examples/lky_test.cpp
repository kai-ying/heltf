#include <iostream>
#include <fftw3.h>

int main()
{
    int N = 100;
    double *in = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(p);

    // Do some processing on the frequency domain data here

    fftw_plan q = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    fftw_execute(q);

    // Do some processing on the time domain data here

    fftw_destroy_plan(p);
    fftw_destroy_plan(q);
    fftw_free(in);
    fftw_free(out);

    return 0;
}