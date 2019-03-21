#ifndef __DEBUGGING__
#define __DEBUGGING__
#include <iomanip>
#include <ostream>

template <typename Iter>
void print(std::ostream& stream, const Iter& begin, const Iter& end) {
    stream << std::setprecision(6);
    for (auto i = begin; i != end; ++i) stream << *i << std::endl;
}

#endif
