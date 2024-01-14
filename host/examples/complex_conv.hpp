#include <iostream>
#include <vector>
#include <complex>

// Function to perform complex signal convolution
std::vector<std::complex<float>> complexConvolution(std::vector<std::complex<float>>& x,
                                                const std::vector<std::complex<float>>& h) {
    int N = x.size();
    int M = h.size();
    int YSize = N + M - 1;

    // Initialize the result vector with zeros
    std::vector<std::complex<float>> y(YSize, 0.0f);

    // Perform convolution
    for (int n = 0; n < YSize; ++n) {
        for (int k = 0; k < N; ++k) {
            if (n - k >= 0 && n - k < M) {
                y[n] += x[k] * h[n - k];
            }
        }
    }

    return y;
}