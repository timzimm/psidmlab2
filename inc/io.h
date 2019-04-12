#ifndef __IO__
#define __IO__
#include <string>
#include "blaze/math/DynamicVector.h"
#include "hdf5.h"

// Forward Declare
class Parameters;
class SimState;

// TODO Switch to H5CPP to allow for blaze vectors
class OutputFile {
   private:
    hid_t file;
    hid_t complex_type;

   public:
    OutputFile(const Parameters& p);
    void write(const SimState& state, const Parameters& p);
    void read(const SimState& state, const Parameters& p) const;

    // Used for debugging
    void write(const blaze::DynamicVector<double>& data,
               const std::string& dset_name);
    ~OutputFile();
};
#endif
