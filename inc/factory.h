#ifndef __FACTORY__
#define __FACTORY__

#include "logging.h"

#include <boost/core/typeinfo.hpp>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <utility>

// CRTP. Inject Base class dependent Factory interface into Base class (i.e. the
// ABC) by inheriting via
//                     class Base: public Factory<Base, ...>
// In short: each ABC has (static) Factory functionality right away. This is
// different than the standard self-registering factory pattern where the ABC
// and the factory constitute two seperate identities.
template <class Base, class... Args>
class Factory {
   public:
    template <class... T>
    static std::unique_ptr<Base> make(const std::string& s, T&&... args) {
        try {
            return data().at(s)(std::forward<T>(args)...);
        } catch (std::out_of_range) {
            std::cout << ERRORTAG(s + " not registered") << std::endl;
            exit(1);
        }
    }

    // Avoids class Base1 : Factory<Base2, ...>. See also private c'tor
    friend Base;

    // Self register functionality injected into Derived_from_Base via
    // class Derived_from_Base: public Base::Registrar<Derived_from_Base>
    template <class Derived_from_Base>
    struct Registrar : Base {
        friend Derived_from_Base;

        static bool add() {
            // RTTI must be enabled. Use demangled class name as map key.
            const auto name =
                boost::core::demangled_name(typeid(Derived_from_Base));
            Factory::data()[name] = [](Args... args) -> std::unique_ptr<Base> {
                return std::make_unique<Derived_from_Base>(
                    std::forward<Args>(args)...);
            };
            return true;
        }
        static bool registered;

       private:
        // Force instantiation of static member
        Registrar() { (void)registered; }
    };

   private:
    using func_t = std::unique_ptr<Base> (*)(Args...);
    Factory() = default;

    static auto& data() {
        static std::unordered_map<std::string, func_t> s;
        return s;
    }
};

// Registering all derived classes.
template <class Base, class... Args>
template <class Derived_from_Base>
bool Factory<Base, Args...>::Registrar<Derived_from_Base>::registered =
    Factory<Base, Args...>::Registrar<Derived_from_Base>::add();

#endif
