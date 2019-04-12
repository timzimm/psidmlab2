#ifndef __PARAMETERS__
#define __PARAMETERS__
#include <boost/lexical_cast.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include "logging.h"

using namespace boost::multi_index;

class Parameters {
    struct Parameter_t {
        std::string key;
        std::string value;
    };

    using ParameterMultiMap_t = multi_index_container<
        Parameter_t,
        indexed_by<sequenced<>, hashed_unique<member<Parameter_t, std::string,
                                                     &Parameter_t::key>>>>;
    ParameterMultiMap_t vault;

   public:
    Parameters(const std::string& path_to_ini) {
        using namespace boost::property_tree;

        // Tree representation of the input file
        ptree tree;

        // Parse the ini file
        ini_parser::read_ini(path_to_ini, tree);

        // Parameters is a "flat" representation of ptree.
        // To get to this flattened hierarchy, we need to iterate over
        // the individual ini sections and insert the paramters into our
        // vault.
        for (const auto& section : tree) {
            for (const auto& parameter : section.second) {
                insert(parameter.first, parameter.second.data());
            }
        }
    }

    // Write access to container via "unordered_map" interface
    void insert(const std::string& key, const std::string& value) {
        auto& string_index = vault.template get<1>();
        string_index.insert({key, value});
    }

    // Read access to container via "unordered_map" interface
    template <typename T>
    void get(const std::string& key, T& var) const {
        using boost::bad_lexical_cast;
        using boost::lexical_cast;
        auto& string_index = vault.template get<1>();
        auto iter = string_index.find(key);
        // Bail out if key not present
        if (iter == string_index.end()) {
            std::cout << ERRORTAG("Parameter " << key << " not found")
                      << std::endl;
            exit(1);
        }
        // Bail out if cast fails
        try {
            var = lexical_cast<T>(iter->value);
        } catch (bad_lexical_cast) {
            std::cout << ERRORTAG("Cast not possible") << std::endl;
            exit(1);
        }
    }
    template <typename T>
    T get(const std::string& key) const {
        T temp;
        get(key, temp);
        return temp;
    }

    // Print to stream
    friend std::ostream& operator<<(std::ostream& s, const Parameters& p) {
        // Printing to stream should be in the order we read in the INI file.
        // Otherwise it can be quiet confusing. Therefore, we employ vaults
        // list-like interface for which only the insertion order matters.
        auto& list_index = p.vault.template get<0>();
        for (const auto& p_struct : list_index) {
            s << std::left << std::setw(24) << p_struct.key << "|" << std::right
              << std::setw(24) << p_struct.value << std::endl;
        };

        return s;
    }
};

#endif
