#ifndef __IO__
#define __IO__
#include <string>
#include <vector>
#include "blaze/math/DynamicVector.h"
#include "hdf5.h"

class HDF5File {
   private:
    hid_t file;          // Handle to HDF5 file
    hid_t complex_type;  // compound datatype for complex data

    // Strip off padding
    template <typename T>
    void matrix_dataspace_dim(const blaze::DynamicMatrix<T, blaze::rowMajor>& M,
                              hsize_t* dims) {
        dims[0] = M.rows();
        dims[1] = M.spacing();
    }

    template <typename T>
    void matrix_dataspace_dim(
        const blaze::DynamicMatrix<T, blaze::columnMajor>& M, hsize_t* dims) {
        dims[0] = M.columns();
        dims[1] = M.spacing();
    }

    template <typename T>
    void H5Dwrite_wrapper(hid_t dataset_id, hid_t type, hid_t mem_space_id,
                          hid_t file_space_id, hid_t xfer_plist_id,
                          const blaze::DynamicMatrix<T, blaze::rowMajor>& buf) {
        H5Dwrite(dataset_id, type, mem_space_id, file_space_id, xfer_plist_id,
                 buf.data());
    }

    // Column-Major matrix needs to be "transposed" onn write
    template <typename T>
    void H5Dwrite_wrapper(
        hid_t dataset_id, hid_t type, hid_t mem_space_id, hid_t file_space_id,
        hid_t xfer_plist_id,
        const blaze::DynamicMatrix<T, blaze::columnMajor>& buf) {
        // Write column by column ( otherwise the raw pointer is invalid)
        hsize_t offset_file[] = {0, 0};
        hsize_t count_file[] = {buf.rows(), 1};

        hsize_t offset_mem[] = {0, 0};
        hsize_t count_mem[] = {1, buf.rows()};
        for (int i = 0; i < buf.columns(); ++i) {
            // Select column in file
            offset_file[1] = i;
            H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset_file,
                                NULL, count_file, NULL);

            // Select row in memory
            offset_mem[0] = i;
            H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, offset_mem, NULL,
                                count_mem, NULL);

            H5Dwrite(dataset_id, type, mem_space_id, file_space_id,
                     xfer_plist_id, buf.data());
        }
    }

    template <typename T>
    hid_t H5Dcreate_wrapper(hid_t loc_id, const char* name, hid_t space_id,
                            hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id) {
        return H5Dcreate(loc_id, name, H5T_NATIVE_DOUBLE, space_id, lcpl_id,
                         dcpl_id, dapl_id);
    }

    hid_t create_groups_along_path(const std::string& path) {
        size_t end = 0;
        herr_t status;
        std::string subpath;
        while ((end = path.find("/", end + 1)) != std::string::npos) {
            subpath = path.substr(0, end);
            auto exists = H5Lexists(file, subpath.c_str(), H5P_DEFAULT);
            if (exists == 0) {
                auto group = H5Gcreate(file, subpath.c_str(), H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);
                H5Gclose(group);
            }
        }
        return H5Gopen(file, subpath.c_str(), H5P_DEFAULT);
    }

   public:
    HDF5File(const std::string& filename,
             const std::vector<std::string>& init_groups)
        : file(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT)),
          complex_type(H5Tcreate(H5T_COMPOUND, 2 * sizeof(double))) {
        // HDF5 doesn't know about complex numbers. Let's define them
        H5Tinsert(complex_type, "r", 0, H5T_NATIVE_DOUBLE);
        H5Tinsert(complex_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);

        for (const auto& group_name : init_groups) {
            auto group = H5Gcreate(file, group_name.c_str(), H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose(group);
        }
    }

    template <typename T>
    void write(const std::string& ds_path, const blaze::DynamicVector<T> data) {
        auto parent_group = create_groups_along_path(ds_path);

        // Create new dataspace and dataset for data vector
        const int rank = 1;
        hsize_t dim = data.size();

        // Discards padding automatically
        auto dataspace = H5Screate_simple(rank, &dim, NULL);
        auto dataset =
            H5Dcreate_wrapper<T>(parent_group, ds_path.c_str(), dataspace,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Write data into dataset
        // Query dataset type
        auto type = H5Dget_type(dataset);
        H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close all handles
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    template <typename T, bool SO>
    void write(const std::string& ds_path,
               const blaze::DynamicMatrix<T, SO> data) {
        auto parent_group = create_groups_along_path(ds_path);

        // Create new dataspaces in both the file and the memory area
        const int rank = 2;
        hsize_t dim_file[] = {data.rows(), data.columns()};
        auto dataspace_file = H5Screate_simple(rank, dim_file, NULL);
        auto dataset =
            H5Dcreate_wrapper<T>(parent_group, ds_path.c_str(), dataspace_file,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t dim[2];
        // Mem and file dataspace are not necessarily the same due to padding
        matrix_dataspace_dim(data, dim);
        auto dataspace_matrix = H5Screate_simple(rank, dim, NULL);
        hsize_t offset[] = {0, 0};
        // Discard padding
        H5Sselect_hyperslab(dataspace_matrix, H5S_SELECT_SET, offset, NULL,
                            dim_file, NULL);

        // Write data into dataset
        // Query dataset type
        auto type = H5Dget_type(dataset);
        H5Dwrite_wrapper(dataset, type, dataspace_matrix, dataspace_file,
                         H5P_DEFAULT, data);

        // Close all handles
        H5Sclose(dataspace_matrix);
        H5Sclose(dataspace_file);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    ~HDF5File() { H5Fclose(file); };
};

template <>
hid_t HDF5File::H5Dcreate_wrapper<std::complex<double>>(
    hid_t loc_id, const char* name, hid_t space_id, hid_t lcpl_id,
    hid_t dcpl_id, hid_t dapl_id) {
    return H5Dcreate(loc_id, name, complex_type, space_id, lcpl_id, dcpl_id,
                     dapl_id);
}

#endif
