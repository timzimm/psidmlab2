#ifndef __IO__
#define __IO__

#include "H5Apublic.h"
#include "H5Fpublic.h"
#include "H5Ppublic.h"
#include "logging.h"

#include <blaze/math/DynamicVector.h>
#include <blaze/util/typetraits/IsBoolean.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsComplexDouble.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsInteger.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <hdf5.h>
#include <cstring>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

template <typename Head, typename... Tail>
struct Last {
    using T = typename Last<Tail...>::T;
};

template <typename Head>
struct Last<Head> {
    using T = Head;
};

// Reads columns of stream s into the supplied vectors.
// Passed in vectors already have to be large enough to hold the columns
// of s. If Container holds complex elmenents, even numbered (0-based) columns
// are interpreted as real part and odd numbered columns as imaginary parts.
template <typename... Ts>
void fill_from_file(std::istream& s, Ts&... vecs) {
    using ElementType = blaze::RemoveReference_t<decltype(
        std::declval<typename Last<Ts...>::T>()[0])>;

    const int N = std::max({std::size(vecs)...});
    if constexpr (blaze::IsComplex_v<ElementType>) {
        auto real = [](std::complex<double>& z) -> double& {
            return reinterpret_cast<double(&)[2]>(z)[0];
        };

        auto imag = [](std::complex<double>& z) -> double& {
            return reinterpret_cast<double(&)[2]>(z)[1];
        };
        for (int i = 0; i < N; ++i) {
            double drop;
            size_t pos = s.tellg();
            (void(s >> real(vecs[i]) >> drop), ...);
            s.seekg(pos);
            (void(s >> drop >> imag(vecs[i])), ...);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            (s >> ... >> vecs[i]);
        }
    }
}

// Type Traits
template <typename T>
struct is_string : public std::false_type {};
template <>
struct is_string<std::string> : public std::true_type {};

class HDF5File {
   private:
    hid_t file;          // Handle to HDF5 file
    hid_t complex_type;  // compound datatype for complex data

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
    enum class Access { Read, Write };
    HDF5File(const std::string& filename, const Access& strategy) {
        // Store old error handler
        herr_t (*old_func)(void*);
        void* old_client_data;
        H5Eget_auto1(&old_func, &old_client_data);
        // Turn off error stack printing
        /* H5Eset_auto1(NULL, NULL); */

        if (strategy == Access::Write) {
            file = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT,
                             H5P_DEFAULT);
            if (file < 0) {
                // File exists. Open it in rw mode
                file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                if (file < 0) {
                    std::cerr
                        << ERRORTAG("Cannot open file for writing. Damaged?")
                        << std::endl;
                    exit(1);
                }
            }
        }
        if (strategy == Access::Read) {
            // Read implies the file already exists. Open it i r mode.
            file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file < 0) {
                std::cerr << ERRORTAG("Cannot open file for reading. Damaged?")
                          << std::endl;
                exit(1);
            }
        }

        H5Eset_auto1(old_func, old_client_data);

        // HDF5 doesn't know about complex numbers. Let's define them
        complex_type = H5Tcreate(H5T_COMPOUND, 2 * sizeof(double));
        H5Tinsert(complex_type, "r", 0, H5T_NATIVE_DOUBLE);
        H5Tinsert(complex_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);
    }

    // Write generic vector to file at path ds_path. In this context generic
    // means:
    //
    // T = std::complex<double> , double, int
    // TF = blaze::columnVector, blaze::rowVector
    //
    // If groups along ds_path are missing, we generate them in the process

    template <typename T, bool TF>
    void write(const std::string& ds_path,
               const blaze::DynamicVector<T, TF>& data) {
        // We only deal with Ts as specified above
        static_assert(blaze::IsComplexDouble_v<T> || blaze::IsDouble_v<T> ||
                          (blaze::IsInteger_v<T> && std::is_signed_v<T>),
                      "Datatype not supported");

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

        // Select correct memory and file types
        hid_t filetype, memtype;
        if constexpr (blaze::IsComplexDouble_v<T>) {
            filetype = complex_type;
            memtype = complex_type;
        } else if constexpr (blaze::IsDouble_v<T>) {
            filetype = H5T_IEEE_F64BE;
            memtype = H5T_NATIVE_DOUBLE;
        } else {
            filetype = H5T_STD_I32BE;
            memtype = H5T_NATIVE_INT;
        }
        auto dataset =
            H5Dcreate(parent_group, ds_path.c_str(), filetype, dataspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close all handles
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Gclose(parent_group);
    }

    // Write generic matrix to file. In this context generic means:
    //
    // T = std::complex<double> , double, int
    // SO = blaze::rowMajor, blaze::columnMajor
    //
    // If groups along ds_path are missing, we generate them in the process

    template <typename T, bool SO>
    void write(const std::string& ds_path,
               const blaze::DynamicMatrix<T, SO>& data) {
        // We only deal with Ts as specified above
        static_assert(blaze::IsComplexDouble_v<T> || blaze::IsDouble_v<T> ||
                          (blaze::IsInteger_v<T> && std::is_signed_v<T>),
                      "Datatype not supported");

        // Split path into parent group and dataset name
        auto split_pos = ds_path.find_last_of("/");
        auto parent_name = ds_path.substr(0, split_pos + 1);
        auto set_name = ds_path.substr(split_pos + 1);

        auto parent_group = create_groups_along_path(ds_path);

        // Create new dataspaces in both the file and the memory area
        const int rank = 2;

        hsize_t dim_file[] = {
            (SO == blaze::columnMajor) ? data.columns() : data.rows(),
            (SO == blaze::columnMajor) ? data.rows() : data.columns()};

        // Mem and file dataspace are not necessarily the same due to
        // padding.
        hsize_t dim[2];
        matrix_dataspace_dim(data, dim);

        auto dataspace_file = H5Screate_simple(rank, dim_file, NULL);
        auto dataspace_matrix = H5Screate_simple(rank, dim, NULL);

        // Discard padding in memory
        const hsize_t offset[] = {0, 0};
        H5Sselect_hyperslab(dataspace_matrix, H5S_SELECT_SET, offset, NULL,
                            dim_file, NULL);

        // Select correct memory and file types
        hid_t filetype, memtype;
        if constexpr (blaze::IsComplexDouble_v<T>) {
            filetype = complex_type;
            memtype = complex_type;
        } else if constexpr (blaze::IsDouble_v<T>) {
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
    // T = std::complex<double>, double, (unsigned) int, bool, std::string
    //
    // If the object does not exist yet, we generate a group and attach the
    // attribute to it.
    template <typename T>
    void write_scalar_attribute(const std::string& path,
                                const std::string& name, const T& value,
                                hid_t loc = -1) {
        // We only deal with Ts as specified above
        static_assert(blaze::IsComplexDouble_v<T> || blaze::IsDouble_v<T> ||
                          blaze::IsInteger_v<T> || blaze::IsBoolean_v<T> ||
                          is_string<T>(),
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
        if constexpr (blaze::IsComplexDouble_v<T>) {
            filetype = complex_type;
            memtype = complex_type;
        } else if constexpr (blaze::IsDouble_v<T>) {
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
    // Return vector of all object names linked to root
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
        H5Literate(root_obj, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                   push_back_names, &content);

        std::string newroot = (root == "/") ? "" : root;
        for (auto& s : content) {
            s = newroot + "/" + s;
        }

        H5Oclose(root_obj);

        return content;
    }
    // Reads vector at path into double precision blaze column-vector
    blaze::DynamicVector<double> read_vector(const std::string& path) const {
        hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);

        // Vector => rank=1
        int rank = 1;
        hsize_t dim_file[rank];
        hsize_t dim_mem[rank];

        // File side
        hid_t file_space = H5Dget_space(ds);
        H5Sget_simple_extent_dims(file_space, dim_file, nullptr);

        // Memory side
        blaze::DynamicVector<double> data(dim_file[0]);
        dim_mem[0] = dim_file[0];
        hid_t mem_space = H5Screate_simple(rank, dim_mem, nullptr);

        H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                data.data());

        H5Sclose(file_space);
        H5Sclose(mem_space);
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

    ~HDF5File() {
        H5Fclose(file);
        H5Tclose(complex_type);
    };
};

#endif
