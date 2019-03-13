#include "io.h"
#include "psidm.h"

// TODO: Add simulation parameters to constructor
// TODO: Add write overload for analysis function
// TODO: MPI
// TODO: read function
OutputFile::OutputFile(const Parameters& params) {
    file = H5Fcreate(params.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                     H5P_DEFAULT);
    auto state_group = H5Gcreate(file, "SimulationState", H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    auto obs_group =
        H5Gcreate(file, "Observables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // HDF5 doesn't know about complex numbers. Let's define them
    complex_type = H5Tcreate(H5T_COMPOUND, 2 * sizeof(double));
    H5Tinsert(complex_type, "r", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);

    H5Gclose(obs_group);
    H5Gclose(state_group);
}

// Creates new group n (current step number), appends two datasets to the
// n_group and inserts the simulation state.
void OutputFile::write(const SimState& state, const Parameters& params) {
    // Step 1: Create new group n to hold potentials and wavefunctions
    auto n_name = std::string("/SimulationState/") + std::to_string(state.n);
    auto n_group =
        H5Gcreate(file, n_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Step 2: Add step information to the n_group
    hsize_t dims = 2;
    auto dataspace_attr = H5Screate_simple(1, &dims, NULL);
    double attr_buffer[] = {state.t, state.a};
    auto attr_time = H5Acreate2(n_group, "tau a(tau)", H5T_NATIVE_DOUBLE,
                                dataspace_attr, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &attr_buffer);

    // Step 3: Create new datasets in n_group
    // Dataspace is the same for psis and Vs. The only difference is the
    // datatype
    const int rank = 2;
    hsize_t dim_psiVs[] = {params.M, params.N};
    auto dataspace_psiVs = H5Screate_simple(rank, dim_psiVs, NULL);

    auto ds_Vs_name = n_name + "/Vs";
    auto ds_Vs =
        H5Dcreate(file, ds_Vs_name.c_str(), H5T_NATIVE_DOUBLE, dataspace_psiVs,
                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto ds_psis_name = n_name + "/psis";
    auto ds_psis =
        H5Dcreate(file, ds_psis_name.c_str(), complex_type, dataspace_psiVs,
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
}

OutputFile::~OutputFile() { H5Fclose(file); }
