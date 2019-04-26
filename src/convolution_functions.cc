#include "convolution_functions.h"

void factorize(const int n, int* n_factors, int factors[],
               int* implemented_factors) {
    int nf = 0;
    int ntest = n;
    int factor;
    int i = 0;

    if (n == 1) {
        factors[0] = 1;
        *n_factors = 1;
        return;
    }

    while (implemented_factors[i] && ntest != 1) {
        factor = implemented_factors[i];
        while ((ntest % factor) == 0) {
            ntest = ntest / factor;
            factors[nf] = factor;
            nf++;
        }
        i++;
    }

    if (ntest != 1) {
        factors[nf] = ntest;
        nf++;
    }

    {
        int product = 1;

        for (i = 0; i < nf; i++) {
            product *= factors[i];
        }
    }

    *n_factors = nf;
}

bool is_optimal(int n, int* implemented_factors) {
    // We check that n is not a multiple of 4*4*4*2
    if (n % 4 * 4 * 4 * 2 == 0) return false;

    int nf = 0;
    int factors[64];
    int i = 0;
    factorize(n, &nf, factors, implemented_factors);

    while (implemented_factors[i]) {
        if (factors[nf - 1] == implemented_factors[i]) return true;
        ++i;
    }
    return false;
}

int find_closest_factor(int n) {
    // FFTW optimizes for these factors
    int implemented_factors[] = {13, 11, 7, 5, 3, 2, 0};
    int j;
    if (is_optimal(n, implemented_factors)) return n;
    j = n + 1;
    while (!is_optimal(j, implemented_factors)) ++j;
    return j;
}
