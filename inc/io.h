#ifndef __IO__
#define __IO__
#include "common.h"
#include "hdf5.h"

// TODO Switch to H5CPP to allow for blaze vectors
class OutputFile {
   private:
    hid_t file;
    hid_t complex_type;

   public:
    OutputFile(const Parameters& params);
    void write(const SimState& state, const Parameters& params);
    void read(const SimState& state, const Parameters& params) const;

    // Used for debugging
    void write(const std::vector<double>& data, const std::string& dset_name);
    ~OutputFile();
};
#endif
