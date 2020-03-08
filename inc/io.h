#ifndef __IO__
#define __IO__

#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <complex>

template <typename Head, typename... Tail>
struct Last {
    using T = typename Last<Tail...>::T;
};

template <typename Head>
struct Last<Head> {
    using T = Head;
};

// Reads columns of stream s into the supplied vectors.
// Passed in vectors already have to be large enough to hold the columns
// of s. If Container holds complex elmenents, even numbered (0-based) columns
// are interpreted as real part and odd numbered columns as imaginary parts.
template <typename... Ts>
void fill_from_file(std::istream& s, Ts&... vecs) {
    using ElementType = blaze::RemoveReference_t<decltype(
        std::declval<typename Last<Ts...>::T>()[0])>;

    const int N = std::max({std::size(vecs)...});
    if constexpr (blaze::IsComplex_v<ElementType>) {
        auto real = [](std::complex<double>& z) -> double& {
            return reinterpret_cast<double(&)[2]>(z)[0];
        };

        auto imag = [](std::complex<double>& z) -> double& {
            return reinterpret_cast<double(&)[2]>(z)[1];
        };
        for (int i = 0; i < N; ++i) {
            double drop;
            size_t pos = s.tellg();
            (void(s >> real(vecs[i]) >> drop), ...);
            s.seekg(pos);
            (void(s >> drop >> imag(vecs[i])), ...);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            (s >> ... >> vecs[i]);
        }
    }
}

#endif
