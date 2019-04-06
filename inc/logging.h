#ifndef __DEBUGGING__
#define __DEBUGGING__
#include <iomanip>
#include <ostream>

template <typename Iter>
void print(std::ostream& stream, const Iter& begin, const Iter& end) {
    stream << std::setprecision(6);
    for (auto i = begin; i != end; ++i) stream << *i << std::endl;
}

// Some convenience macros and typedefs no one wants to write out
#define INFOTAG(message) "\033[1;34m[INFO]\033[0m " << message
#define WARNINGTAG(message) "\033[1;33m[WARNING]\033[0m " << message
#define ERRORTAG(message) "\033[1;33m[ERROR]\033[0m " << message

#endif
