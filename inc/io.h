#ifndef __IO__
#define __IO__
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include "blaze/math/DynamicVector.h"
#include "hdf5.h"

// Type Traits
template <typename T>
struct is_complex : public std::false_type {};
template <typename T>
struct is_complex<std::complex<T>> : public std::true_type {};

template <typename T>
struct is_string : public std::false_type {};
template <>
struct is_string<std::string> : public std::true_type {};

// Definition of the File Interface. The project uses HDF5 as file format. HDF5
// is the quasi-standard in scientific computing. It handles multidimensional,
// heterogeneous data in a self-explaining way, very similar to a file structure
// in your OS. Moreover, it is fast and can handle async, parallel, compressed
// I/O if required (currently not).

class HDF5File {
   private:
    hid_t file;          // Handle to HDF5 file
    hid_t complex_type;  // compound datatype for complex data
    hid_t str_type;      // compound datatype for string data

    // Computes correct memory dataspace due to potential padding of blaze
    // matrices.
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

    // HDF5 always assumes row major data. This is of course wrong for blaze's
    // columnMajor matrices. Hence, a distinction has to be made indepedent of
    // the type.
    // Row-Major matrix is trivial. Just forward to the C-API.
    template <typename T>
    void H5Dwrite_wrapper(hid_t dataset_id, hid_t type, hid_t mem_space_id,
                          hid_t file_space_id, hid_t xfer_plist_id,
                          const blaze::DynamicMatrix<T, blaze::rowMajor>& buf) {
        H5Dwrite(dataset_id, type, mem_space_id, file_space_id, xfer_plist_id,
                 buf.data());
    }

    // Column-Major matrix needs to be "transposed" on write
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

    // mkdir -p for HDF5
    // TODO Regex
    hid_t create_groups_along_path(const std::string& path) {
        size_t end = 0;
        herr_t status;
        std::string subpath = "/";
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
             const std::vector<std::string>& init_groups = {})
        : file(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT)),
          complex_type(H5Tcreate(H5T_COMPOUND, 2 * sizeof(double))),
          str_type(H5Tcopy(H5T_C_S1)) {
        // HDF5 doesn't know about complex numbers. Let's define them
        H5Tinsert(complex_type, "r", 0, H5T_NATIVE_DOUBLE);
        H5Tinsert(complex_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);

        // Same goes for variable strings
        /* H5Tset_size(str_type, H5T_VARIABLE); */

        for (const auto& group_name : init_groups) {
            auto group = H5Gcreate(file, group_name.c_str(), H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose(group);
        }
    }

    // Write generic vector to file at path ds_path. In this context generic
    // means:
    //
    // T = std::complex<FT> , FT
    // FT = double, float, int... (anything with sizeof(FT) < sizeof(double))
    // TF = blaze::columnVector, blaze::rowVector
    //
    // If groups along ds_path are missing, we generate them in the process

    template <typename T, bool TF>
    void write(const std::string& ds_path,
               const blaze::DynamicVector<T, TF>& data) {
        // Split path into parent group and dataset name
        auto split_pos = ds_path.find_last_of("/");
        auto parent_name = ds_path.substr(0, split_pos + 1);
        auto set_name = ds_path.substr(split_pos + 1);

        auto parent_group = create_groups_along_path(parent_name);

        // Create new dataspace and dataset for data vector
        const int rank = 1;
        hsize_t dim = data.size();

        // Discards padding automatically
        auto dataspace = H5Screate_simple(rank, &dim, NULL);

        hid_t datatype;
        if constexpr (std::is_floating_point_v<T>)
            datatype = H5T_NATIVE_DOUBLE;
        else if constexpr (is_complex<T>())
            datatype = complex_type;
        auto dataset =
            H5Dcreate(parent_group, ds_path.c_str(), datatype, dataspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close all handles
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    // Write generic matrix to file. In this context generic means:
    //
    // T = std::complex<FT> , FT
    // FT = double, float, int... (anything with sizeof(FT) < sizeof(double))
    // SO = blaze::rowMajor, blaze::columnMajor

    template <typename T, bool SO>
    void write(const std::string& ds_path,
               const blaze::DynamicMatrix<T, SO>& data) {
        // Split path into parent group and dataset name
        auto split_pos = ds_path.find_last_of("/");
        auto parent_name = ds_path.substr(0, split_pos + 1);
        auto set_name = ds_path.substr(split_pos + 1);

        auto parent_group = create_groups_along_path(ds_path);

        // Create new dataspaces in both the file and the memory area
        const int rank = 2;
        hsize_t dim_file[] = {data.rows(), data.columns()};
        auto dataspace_file = H5Screate_simple(rank, dim_file, NULL);

        hid_t datatype;
        if constexpr (std::is_floating_point_v<T>)
            datatype = H5T_NATIVE_DOUBLE;
        else if constexpr (is_complex<T>())
            datatype = complex_type;
        auto dataset =
            H5Dcreate(parent_group, ds_path.c_str(), datatype, dataspace_file,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t dim[2];
        // Mem and file dataspace are not necessarily the same due to padding
        matrix_dataspace_dim(data, dim);
        auto dataspace_matrix = H5Screate_simple(rank, dim, NULL);
        hsize_t offset[] = {0, 0};
        // Discard padding
        H5Sselect_hyperslab(dataspace_matrix, H5S_SELECT_SET, offset, NULL,
                            dim_file, NULL);

        H5Dwrite_wrapper(dataset, datatype, dataspace_matrix, dataspace_file,
                         H5P_DEFAULT, data);

        // Close all handles
        H5Sclose(dataspace_matrix);
        H5Sclose(dataspace_file);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    // Adds scalar attribute of type T to object (group, dataset) at path. If
    // the object does not exist yet, we generate a group and attach the
    // attribute to it.
    template <typename T>
    void add_scalar_attribute(const std::string& path, const std::string& name,
                              const T& value, hid_t loc = -1) {
        // We need linear order of members within class types
        static_assert(std::is_standard_layout<T>::value);

        bool ext_managed = true;
        if (loc == -1) {
            auto exists = H5Lexists(file, path.c_str(), H5P_DEFAULT);
            // Open object whatever it is
            if (exists) loc = H5Oopen(file, path.c_str(), H5P_DEFAULT);
            // Assume caller wants a group
            else {
                loc = create_groups_along_path(path + "/");
                ext_managed = false;
            }
        }

        hsize_t dim = 1;
        hid_t attr_space = H5Screate_simple(1, &dim, nullptr);

        hid_t attr_type = H5T_NATIVE_INT;

        auto c_ptr = reinterpret_cast<const void*>(&value);
        if constexpr (std::is_floating_point_v<T>)
            attr_type = H5T_NATIVE_DOUBLE;
        else if constexpr (is_complex<T>())
            attr_type = complex_type;
        else if constexpr (is_string<T>()) {
            attr_type = str_type;
            H5Tset_size(str_type, value.size());
            c_ptr = reinterpret_cast<const void*>(value.c_str());
        }

        hid_t attr = H5Acreate(loc, name.c_str(), attr_type, attr_space,
                               H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, attr_type, c_ptr);

        H5Aclose(attr);

        if (!ext_managed) H5Gclose(loc);
    }

    template <typename T>
    void add_scalar_attribute(const std::string& path,
                              const std::map<std::string, T>& attr_map) {
        hid_t loc;
        auto exists = H5Lexists(file, path.c_str(), H5P_DEFAULT);
        // Open object whatever it is
        if (exists) loc = H5Oopen(file, path.c_str(), H5P_DEFAULT);
        // Assume caller wants a group
        else
            loc = create_groups_along_path(path + "/");
        hsize_t dim = 1;
        auto attr_space = H5Screate_simple(1, &dim, nullptr);

        for (const auto& pair : attr_map) {
            const std::string& key = pair.first;
            const T& val = pair.second;
            add_scalar_attribute(path, key, val, loc);
        }

        H5Gclose(loc);
    }

    ~HDF5File() { H5Fclose(file); };
};

#endif
