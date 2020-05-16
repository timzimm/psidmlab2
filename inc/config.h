#ifndef PSIDMLAB_CONFIG
#define PSIDMLAB_CONFIG
#include <limits>

#define PSIDMLAB_MINIMUM_INTEGRATION_STEP \
    100 * std::numeric_limits<double>::epsilon()
#define PSIDMLAB_RENORMALIZATION_THRESHOLD 1000
#define PSIDMLAB_CONVOLUTION_THRESHOLD 20

#endif
