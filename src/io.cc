#include "io.h"
#include "common.h"

// TODO: Add simulation parameters to constructor
// TODO: Add write overload for analysis function
// TODO: MPI
// TODO: read function
OutputFile::OutputFile(const Parameters& params) {
    file = H5Fcreate(params.out_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                     H5P_DEFAULT);
    auto state_group = H5Gcreate(file, "SimulationState", H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    auto obs_group =
        H5Gcreate(file, "Observables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#ifndef NDEBUG
    auto debug_group =
        H5Gcreate(file, "Debug", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    // HDF5 doesn't know about complex numbers. Let's define them
    complex_type = H5Tcreate(H5T_COMPOUND, 2 * sizeof(double));
    H5Tinsert(complex_type, "r", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);

    H5Gclose(obs_group);
    H5Gclose(state_group);

#ifndef NDEBUG
    H5Gclose(debug_group);
#endif
}

// Creates new group n (current step number), appends two datasets to the
// n_group and inserts the simulation state.
void OutputFile::write(const SimState& state, const Parameters& params) {
    // Step 1: Create new group n to hold potentials and wavefunctions
    auto state_group = H5Gopen(file, "SimulationState", H5P_DEFAULT);
    auto n_group = H5Gcreate(state_group, std::to_string(state.n).c_str(),
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Step 2: Add step information to the n_group
    hsize_t dims = 2;
    auto dataspace_attr = H5Screate_simple(1, &dims, NULL);
    double attr_buffer[] = {state.tau, state.a};
    auto attr_time = H5Acreate2(n_group, "tau a(tau)", H5T_NATIVE_DOUBLE,
                                dataspace_attr, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &attr_buffer);

    // Step 3: Create new datasets in n_group
    // Dataspace is the same for psis and Vs. The only difference is the
    // datatype
    const int rank = 2;
    hsize_t dim_psiVs[] = {params.M, params.N};
    auto dataspace_psiVs = H5Screate_simple(rank, dim_psiVs, NULL);

    auto ds_Vs = H5Dcreate(n_group, "Vs", H5T_NATIVE_DOUBLE, dataspace_psiVs,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto ds_psis = H5Dcreate(n_group, "psis", complex_type, dataspace_psiVs,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write to file
    H5Dwrite(ds_Vs, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             state.Vs.data());
    H5Dwrite(ds_psis, complex_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             state.psis.data());

    // Close all file handles (dataspaces, sets, attributes)
    H5Dclose(ds_Vs);
    H5Sclose(dataspace_psiVs);
    H5Sclose(dataspace_attr);
    H5Dclose(ds_psis);
    H5Aclose(attr_time);
    H5Gclose(n_group);
    H5Gclose(state_group);
}

void OutputFile::write(const blaze::DynamicVector<double>& data,
                       const std::string& dset_name) {
    const int rank = 1;
    hsize_t dim_data = static_cast<hsize_t>(data.size());
    auto dataspace_data = H5Screate_simple(rank, &dim_data, NULL);
    auto debug_group = H5Gopen(file, "Debug", H5P_DEFAULT);
    auto ds_data =
        H5Dcreate(debug_group, dset_name.c_str(), H5T_NATIVE_DOUBLE,
                  dataspace_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             data.data());

    H5Dclose(ds_data);
    H5Sclose(dataspace_data);
    H5Gclose(debug_group);
}

OutputFile::~OutputFile() { H5Fclose(file); }
