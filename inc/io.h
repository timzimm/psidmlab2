#ifndef __IO__
#define __IO__
#include "hdf5.h"
#include "psidm.h"

class OutputFile {
   private:
    hid_t file;
    hid_t complex_type;

   public:
    OutputFile(const Parameters& params);
    void write(const SimState& state, const Parameters& params);
    void read(const SimState& state, const Parameters& params) const;
    ~OutputFile();
};

#endif
