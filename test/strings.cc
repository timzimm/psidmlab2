#include <gtest/gtest.h>
#include <string>

namespace {
TEST(Strings, Equality) {
    std::string s = "Hello World";
    EXPECT_EQ(s, s);
}
TEST(Strings, NotEquality) {
    std::string s = "Hello World";
    EXPECT_NE(s, s + s);
}
}  // namespace
