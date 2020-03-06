#ifndef __FACTORY__
#define __FACTORY__

#include "logging.h"

#include <boost/type_index.hpp>
#include <exception>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

// CRTP. Inject Base class dependent Factory interface into Base class (i.e. the
// ABC) by inheriting via
//                     class Base: public Factory<Base, ...>
// In short: each ABC has (static) Factory functionality right away. This is
// different than the standard self-registering factory pattern where the ABC
// and the factory constitute two seperate identities.
template <typename... Args>
struct pack {};

template <class Base, class... Args>
class Factory {
   public:
    using func_t = std::function<std::unique_ptr<Base>(Args...)>;
    using CtorArgs = pack<Args...>;

    template <class... T>
    static std::unique_ptr<Base> make(const std::string& s, T&&... args) {
        try {
            return data().at(s)(std::forward<T>(args)...);
        } catch (std::out_of_range) {
            std::cerr << ERRORTAG(s + " not registered") << std::endl;
            exit(1);
        }
    }
    static bool add(const std::string& name, const func_t& creator) {
        if (auto it = data().find(name); it == data().end()) {
            Factory::data()[name] = creator;
            return true;
        }
        return false;
    }

    // Avoids class Base1 : Factory<Base2, ...>. See also private c'tor
    friend Base;

   private:
    Factory() = default;

    static auto& data() {
        static std::unordered_map<std::string, func_t> s;
        return s;
    }
};

#define REGISTER(IMPLEMENTATION)                                             \
    inline static const std::string name =                                   \
        boost::typeindex::type_id<IMPLEMENTATION>().pretty_name();           \
    template <typename... Args>                                              \
    static std::function<                                                    \
        std::unique_ptr<typename IMPLEMENTATION::InterfaceType>(Args...)>    \
    create_creator(pack<Args...>) {                                          \
        return                                                               \
            [](Args... args)                                                 \
                -> std::unique_ptr<typename IMPLEMENTATION::InterfaceType> { \
            return std::make_unique<IMPLEMENTATION>(                         \
                std::forward<Args>(args)...);                                \
        };                                                                   \
    };                                                                       \
    inline static bool registered = IMPLEMENTATION::add(                     \
        name,                                                                \
        create_creator(typename IMPLEMENTATION::InterfaceType::CtorArgs()));

#endif
