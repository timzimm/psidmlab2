#ifndef __HDF5_FILE__
#define __HDF5_FILE__

#include "logging.h"

#include <blaze/math/DynamicVector.h>
#include <blaze/util/typetraits/IsBoolean.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsInteger.h>
#include <hdf5.h>
#include <string>
#include <type_traits>
#include <vector>

// Type Traits
template <typename T>
struct is_string : public std::false_type {};
template <>
struct is_string<std::string> : public std::true_type {};

class HDF5File {
   private:
    hid_t file;  // Handle to HDF5 file

    // Computes correct memory dataspace due to potential padding of blaze
    // matrices.
    template <typename T>
    void matrix_dataspace_dim(const blaze::DynamicMatrix<T, blaze::rowMajor>& M,
                              hsize_t* dims) const {
        dims[0] = M.rows();
        dims[1] = M.spacing();
    }

    template <typename T>
    void matrix_dataspace_dim(
        const blaze::DynamicMatrix<T, blaze::columnMajor>& M,
        hsize_t* dims) const {
        dims[0] = M.columns();
        dims[1] = M.spacing();
    }

    // mkdir -p. Caller is responsible to close group!
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
    // File access types. We distinguish:
    // Access::Read     =   Read access to existing file.
    //                      Fails if file is damaged or does not exist.
    // Access::Write    =   Read/Write access to existing file.
    //                      Fails if file is damaged or does not exist.
    // Access::NewFile  =   Create new file.
    //                      Truncates content if file already exists.
    enum class Access { Read, Write, NewFile };

    HDF5File(const std::string& filename, const Access& strategy) {
        // Store old error handler
        herr_t (*old_func)(void*);
        void* old_client_data;
        H5Eget_auto1(&old_func, &old_client_data);
        // Turn off error stack printing
        H5Eset_auto1(nullptr, nullptr);

        // Create new file
        if (strategy == Access::NewFile) {
            file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                             H5P_DEFAULT);
        }
        // Open existing file
        if (strategy == Access::Write) {
            file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (file < 0) {
                std::cerr << ERRORTAG("Cannot open file for writing. Damaged?")
                          << std::endl;
                exit(1);
            }
        }
        // Read existing file
        if (strategy == Access::Read) {
            file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file < 0) {
                std::cerr << ERRORTAG("Cannot open file for reading. Damaged?")
                          << std::endl;
                exit(1);
            }
        }

        H5Eset_auto1(old_func, old_client_data);
    }

    /************************* WRITING FUNCTIONS **************************/

    // Write generic vector to file at path ds_path. In this context generic
    // means:
    //
    // T = double, int
    // TF = blaze::columnVector, blaze::rowVector
    //
    // If groups along ds_path are missing, we generate them in the process
    template <typename T, bool TF>
    void write(const std::string& ds_path,
               const blaze::DynamicVector<T, TF>& data) {
        // We only deal with Ts as specified above
        static_assert(blaze::IsDouble_v<T> ||
                          (blaze::IsInteger_v<T> && std::is_signed_v<T>),
                      "Datatype not supported");

        // Split path into parent group and dataset name
        auto split_pos = ds_path.find_last_of("/");
        auto parent_name = ds_path.substr(0, split_pos + 1);
        auto set_name = ds_path.substr(split_pos + 1);

        auto parent_group = create_groups_along_path(parent_name);

        const int rank = 1;
        hsize_t dim = data.size();
        auto dataspace = H5Screate_simple(rank, &dim, nullptr);

        // Select correct memory and file types
        hid_t filetype, memtype;
        if constexpr (blaze::IsDouble_v<T>) {
            filetype = H5T_IEEE_F64BE;  // 8-byte float, big endian
            memtype = H5T_NATIVE_DOUBLE;
        } else {
            filetype = H5T_STD_I32BE;  // 4-byte integer, big endian
            memtype = H5T_NATIVE_INT;
        }
        auto dataset =
            H5Dcreate(parent_group, ds_path.c_str(), filetype, dataspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Using H5S_all for mem_space_id makes HDF5 use the file dataspace
        // (file_space_id) for the memory data space. Since file_space_id is
        // H5S_ALL as well, the entire dataspace in the file is used. This
        // excludes potential padding automatically.
        H5Dwrite(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close all handles
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    // Write generic matrix to file. In this context generic means:
    //
    // T = double, int
    // SO = blaze::rowMajor, blaze::columnMajor
    //
    // If groups along ds_path are missing, we generate them in the process

    template <typename T, bool SO>
    void write(const std::string& ds_path,
               const blaze::DynamicMatrix<T, SO>& data) {
        // We only deal with Ts as specified above
        static_assert(blaze::IsDouble_v<T> ||
                          (blaze::IsInteger_v<T> && std::is_signed_v<T>),
                      "Datatype not supported");

        // Split path into parent group and dataset name
        auto split_pos = ds_path.find_last_of("/");
        auto parent_name = ds_path.substr(0, split_pos + 1);
        auto set_name = ds_path.substr(split_pos + 1);

        auto parent_group = create_groups_along_path(ds_path);

        const int rank = 2;

        // The file space excludes padding
        hsize_t dim_file[] = {
            (SO == blaze::columnMajor) ? data.columns() : data.rows(),
            (SO == blaze::columnMajor) ? data.rows() : data.columns()};

        // Mem and file dataspace are not necessarily the same due to
        // padding.
        hsize_t dim_mem[2];
        matrix_dataspace_dim(data, dim_mem);

        auto dataspace_file = H5Screate_simple(rank, dim_file, nullptr);
        auto dataspace_matrix = H5Screate_simple(rank, dim_mem, nullptr);

        // Discard padding in memory
        const hsize_t offset[] = {0, 0};
        H5Sselect_hyperslab(dataspace_matrix, H5S_SELECT_SET, offset, nullptr,
                            dim_file, nullptr);

        // Select correct memory and file types
        hid_t filetype, memtype;
        if constexpr (blaze::IsDouble_v<T>) {
            filetype = H5T_IEEE_F64BE;
            memtype = H5T_NATIVE_DOUBLE;
        } else {
            filetype = H5T_STD_I32BE;
            memtype = H5T_NATIVE_INT;
        }

        auto dataset =
            H5Dcreate(parent_group, ds_path.c_str(), filetype, dataspace_file,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset, memtype, dataspace_matrix, dataspace_file,
                 H5P_DEFAULT, data.data());

        // If columnMajor data is written to file, the user needs to be
        // aware that an additional transpose is required
        if (SO == blaze::columnMajor)
            write_scalar_attribute(ds_path, "transpose", true, dataset);
        else
            write_scalar_attribute(ds_path, "transpose", false, dataset);

        // Close all handles
        H5Sclose(dataspace_matrix);
        H5Sclose(dataspace_file);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    // (Over)writes scalar attribute of type T to object (group, dataset) at
    // path.
    //
    // T = double, (unsigned) int, bool, std::string
    //
    // If the object does not exist yet, we generate a group and attach the
    // attribute to it. If the attribute exists, it gets overwritten.
    template <typename T>
    void write_scalar_attribute(const std::string& path,
                                const std::string& name, const T& value,
                                hid_t loc = -1) {
        // We only deal with Ts as specified above
        static_assert(blaze::IsDouble_v<T> || blaze::IsInteger_v<T> ||
                          blaze::IsBoolean_v<T> || is_string<T>(),
                      "Datatype not supported");

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

        // Select correct memory and file types
        hid_t filetype, memtype;
        if constexpr (blaze::IsDouble_v<T>) {
            filetype = H5T_IEEE_F64BE;
            memtype = H5T_NATIVE_DOUBLE;
        } else if constexpr (blaze::IsInteger_v<T> && std::is_signed_v<T>) {
            filetype = H5T_STD_I32BE;
            memtype = H5T_NATIVE_INT;
        } else if constexpr (blaze::IsInteger_v<T> && std::is_unsigned_v<T>) {
            filetype = H5T_STD_U32BE;
            memtype = H5T_NATIVE_UINT;
        } else if constexpr (blaze::IsBoolean_v<T>) {
            filetype = H5T_NATIVE_HBOOL;
            memtype = filetype;
        } else if constexpr (is_string<T>()) {
            hid_t str_type = H5Tcopy(H5T_C_S1);
            H5Tset_size(str_type, value.size());
            filetype = str_type;
            memtype = filetype;
            c_ptr = reinterpret_cast<const void*>(value.c_str());
        }

        // Overwrite if attribute exists
        auto exists = H5Aexists(loc, name.c_str());
        if (exists > 0) {
            H5Adelete(loc, name.c_str());
        }

        hid_t attr = H5Acreate(loc, name.c_str(), filetype, attr_space,
                               H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, memtype, c_ptr);

        if constexpr (is_string<T>()) {
            H5Tclose(memtype);
        }
        H5Aclose(attr);
        if (!ext_managed) H5Oclose(loc);
    }

    /************************* READING FUNCTIONS **************************/

    // Reads vector at path into double precision blaze column-vector
    blaze::DynamicVector<double> read_vector(const std::string& path) const {
        hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);

        // Vector => rank=1
        int rank = 1;
        hsize_t dim_file[rank];

        // File side determines size of memory vector
        hid_t file_space = H5Dget_space(ds);
        H5Sget_simple_extent_dims(file_space, dim_file, nullptr);

        // Memory side
        blaze::DynamicVector<double> data(dim_file[0]);

        H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                data.data());

        H5Sclose(file_space);
        H5Dclose(ds);

        return data;
    }

    // Reads row-major matrix at path into double precision row-major blaze
    // matrix
    blaze::DynamicMatrix<double> read_matrix(const std::string& path) const {
        hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);

        // Matrix => rank=2
        int rank = 2;
        hsize_t dim_file[rank];
        hsize_t dim_mem[rank];

        // File side
        hid_t file_space = H5Dget_space(ds);
        H5Sget_simple_extent_dims(file_space, dim_file, nullptr);

        // Memory side
        blaze::DynamicMatrix<double> data(dim_file[0], dim_file[1]);
        matrix_dataspace_dim(data, dim_mem);
        hid_t mem_space = H5Screate_simple(rank, dim_mem, nullptr);
        hsize_t mem_offset[] = {0, 0};
        H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr,
                            dim_file, nullptr);

        H5Dread(ds, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT,
                data.data());

        H5Sclose(file_space);
        H5Sclose(mem_space);
        H5Dclose(ds);

        return data;
    }

    // Read scalar attribute of type T from object (group, dataset) at path.
    //
    // T = double, (unsigned) int,
    template <typename T>
    T read_scalar_attribute(const std::string& path,
                            const std::string& name) const {
        // We only deal with Ts as specified above
        static_assert(blaze::IsDouble_v<T> || blaze::IsInteger_v<T>,
                      "Datatype not supported");
        hid_t loc = H5Oopen(file, path.c_str(), H5P_DEFAULT);
        hid_t attr = H5Aopen(loc, name.c_str(), H5P_DEFAULT);

        hid_t space = H5Aget_space(attr);

        // Select correct memory file type
        hid_t memtype;
        if constexpr (blaze::IsDouble_v<T>) {
            memtype = H5T_NATIVE_DOUBLE;
        } else if constexpr (blaze::IsInteger_v<T> && std::is_signed_v<T>) {
            memtype = H5T_NATIVE_INT;
        } else if constexpr (blaze::IsInteger_v<T> && std::is_unsigned_v<T>) {
            memtype = H5T_NATIVE_UINT;
        }

        T attribute_val;
        H5Aread(attr, memtype, &attribute_val);

        H5Sclose(space);
        H5Aclose(attr);
        H5Oclose(loc);

        return attribute_val;
    }
    /************************** Miscellaneous **************************/

    // Return vector of absolute paths to all object names linked to root
    // Behaves like UNIX ls /some/path
    std::vector<std::string> ls(const std::string& root) const {
        std::vector<std::string> content;
        // Does root even exist?
        if (!H5Lexists(file, root.c_str(), H5P_DEFAULT)) {
            return content;
        }
        hid_t root_obj = H5Oopen(file, root.c_str(), H5P_DEFAULT);
        auto push_back_names = [](hid_t, const char* name, const H5L_info_t*,
                                  void* data) {
            auto ptr = reinterpret_cast<std::vector<std::string>*>(data);
            ptr->push_back(name);
            return 0;
        };
        H5Literate(root_obj, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
                   push_back_names, &content);

        std::string newroot = (root == "/") ? "" : root;
        for (auto& s : content) {
            s = newroot + "/" + s;
        }

        H5Oclose(root_obj);

        return content;
    }

    ~HDF5File() { H5Fclose(file); };
};

#endif
