#ifndef __SCHROEDINGER_INTERFACE__
#define __SCHROEDINGER_INTERFACE__
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <complex>
#include <memory>
#include <utility>
#include "common.h"
#include "cosmology.h"
#include "fftw3.h"

namespace Schroedinger {

// This defines the common Algorithm interface all Schroedinger integrators will
// share
class Algorithm {
   public:
    virtual void operator()(SimState& state) = 0;
    virtual ~Algorithm() = default;
};

#include "schroedinger/uso_dkd.h"
#include "schroedinger/uso_kdk.h"

// Add new schroedinger algorithms here
enum class AlgorithmType { USO_DKD, USO_KDK };

// Potential solver factory. Add new solvers here.
struct Factory {
    static std::unique_ptr<Algorithm> create(const AlgorithmType& t,
                                             const Parameters& p) {
        switch (t) {
            case AlgorithmType::USO_DKD:
                return std::make_unique<USO_DKD>(p);
            case AlgorithmType::USO_KDK:
                return std::make_unique<USO_KDK>(p);
            default:
                return nullptr;
        }
    }
};
}  // namespace Schroedinger

#endif
