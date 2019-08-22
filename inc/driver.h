#ifndef __DRIVER__
#define __DRIVER__
#include "interfaces.h"

class SimpleDriver : public Driver::Registrar<SimpleDriver> {
    const double dt;

   public:
    SimpleDriver(const Parameters& p);
    void integrate(std::unique_ptr<Stepper>& stepper, SimState& state,
                   const double t_final);
};

class StableDriver : public Driver::Registrar<StableDriver> {
   public:
    StableDriver(const Parameters& p);
    void integrate(std::unique_ptr<Stepper>& stepper, SimState& state,
                   const double t_final);
};

#endif
