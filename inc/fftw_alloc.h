#ifndef __FFTW_ALLOC__
#define __FFTW_ALLOC__

#include <fftw3.h>
#include <memory>

// A custom allocator meant to make fftw usable with std::vector while still
// allowing for proper SIMD allignment
template <class T>
class FFTWAlloc : std::allocator<T> {
   public:
    using value_type = T;
    value_type* allocate(std::size_t n) {
        return static_cast<value_type*>(fftw_malloc(n * sizeof(value_type)));
    }

    void deallocate(value_type* p, std::size_t) noexcept { fftw_free(p); }
};

#endif
