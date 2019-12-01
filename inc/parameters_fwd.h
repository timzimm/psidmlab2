#ifndef __PARAM_FORWARD__
#define __PARAM_FORWARD__
#ifndef INCLUDE_NLOHMANN_JSON_FWD_HPP_
#define INCLUDE_NLOHMANN_JSON_FWD_HPP_

#include <cstdint>  // int64_t, uint64_t
#include <map>
#include <memory>  // allocator
#include <string>
#include <vector>

namespace nlohmann {
template <typename T = void, typename SFINAE = void>
struct adl_serializer;

template <template <typename U, typename V, typename... Args> class ObjectType =
              std::map,
          template <typename U, typename... Args> class ArrayType = std::vector,
          class StringType = std::string, class BooleanType = bool,
          class NumberIntegerType = std::int64_t,
          class NumberUnsignedType = std::uint64_t,
          class NumberFloatType = double,
          template <typename U> class AllocatorType = std::allocator,
          template <typename T, typename SFINAE = void> class JSONSerializer =
              adl_serializer>
class basic_json;

template <typename BasicJsonType>
class json_pointer;

using json = basic_json<>;
}  // namespace nlohmann

#endif  // INCLUDE_NLOHMANN_JSON_FWD_HPP_

using Parameters = nlohmann::json;
#endif

