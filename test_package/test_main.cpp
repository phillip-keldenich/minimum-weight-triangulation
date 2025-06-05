#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <CGAL_MWT/version.h>
#include <doctest/doctest.h>
#include <regex>

TEST_CASE("[cgal-mwt] Test version") {
    std::string version = mwt::get_version();
    std::smatch match;
    std::regex version_regex("[0-9]+\\.[0-9]+\\.[0-9]+");
    CHECK(std::regex_match(version, match, version_regex));
}
