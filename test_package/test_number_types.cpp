#include <CGAL/Exact_rational.h>
#include <CGAL/FPU.h>
#include <CGAL/Interval_nt.h>
#include <CGAL_MWT/CGAL_Rational_aux.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/rational.hpp>
#include <doctest/doctest.h>
#include <gmpxx.h>
#include <random>

TEST_CASE("test CGAL interval construction from large integer") {
    CGAL::Protect_FPU_rounding round_protect;
    CGAL::Interval_nt_advanced interval1(1);
    CHECK(interval1.inf() == 1.0);
    CHECK(interval1.sup() == 1.0);

    std::int64_t large_int = INT64_C(9007199254740993);
    CGAL::Interval_nt_advanced interval2 = large_int;
    CHECK(interval2.inf() <= 9007199254740992.0);
    CHECK(interval2.sup() >= 9007199254740994.0);

    std::uint64_t large_uint = UINT64_C(9223372036854827931);
    CGAL::Interval_nt_advanced interval3(large_uint);
    CHECK(interval3.inf() <= 9223372036854827008.0);
    CHECK(interval3.sup() >= 9223372036854829056.0);
}

TEST_CASE_TEMPLATE("test CGAL::Exact_rational auxiliary functions (different possible implementations)", NT,
                   CGAL::Exact_rational, boost::multiprecision::cpp_rational, mpq_class,
                   CGAL::Quotient<boost::multiprecision::cpp_int>) {
    CHECK(mwt::rational_is_integer(NT(0)));
    CHECK(mwt::rational_is_zero(NT(0)));
    CHECK(!mwt::rational_is_zero(NT(2.5)));
    CHECK(!mwt::rational_is_integer(NT(2.5)));
    CHECK(mwt::rational_is_integer(NT(2)));
    CHECK(mwt::rational_to_platform_int(NT(2)));
    CHECK(mwt::rational_to_platform_int(NT(2)).value() == 2);
    NT huge_val = NT(std::numeric_limits<long>::max());
    huge_val *= huge_val;
    huge_val *= huge_val;
    huge_val *= huge_val;
    huge_val += huge_val;
    huge_val /= NT(std::numeric_limits<long>::max());
    CHECK(!mwt::rational_to_platform_int(huge_val));
    CHECK(!mwt::rational_is_zero(huge_val));
    CHECK(mwt::rational_is_integer(huge_val));
}

TEST_CASE("test conversions from int64 to CGAL type with square roots") {
    std::int64_t large_val1 = INT64_C(4611686018434622827);
    std::int64_t large_val2 = INT64_C(-4611686018434622827);
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_ft1 = mwt::rational_to_sqrt_nt(large_val1);
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_ft2 = mwt::rational_to_sqrt_nt(large_val2);
    CORE::BigInt large_val_bi("-4611686018434622827");
    CHECK(-large_val_bi == cgal_ft1);
    CHECK(large_val_bi == cgal_ft2);

    std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<std::int64_t> dist(std::numeric_limits<std::int64_t>::min(),
                                                     std::numeric_limits<std::int64_t>::max());
    for(std::size_t i = 0; i < 10'000; ++i) {
        std::int64_t val = dist(gen);
        std::ostringstream oss;
        oss << val;
        CORE::BigInt bi(oss.str());
        CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_num = mwt::rational_to_sqrt_nt(val);
        CHECK(bi == cgal_num);
    }
}

TEST_CASE_TEMPLATE("test conversions to CGAL type with square roots", NT, CGAL::Exact_rational,
                   boost::multiprecision::cpp_rational, mpq_class, CGAL::Quotient<boost::multiprecision::cpp_int>) {
    NT zero(0);
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_zero = mwt::rational_to_sqrt_nt(zero);
    CHECK(cgal_zero == 0);

    NT one(1);
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_one = mwt::rational_to_sqrt_nt(one);
    CHECK(cgal_one == 1);

    NT numden(11, 17);
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_numden = mwt::rational_to_sqrt_nt(numden);
    CHECK(cgal_numden * 17 == 11);

    NT largenum("29485739458039487510394019834570139485713094857130495831745134513");
    NT largeden("348573279458423429837429854279845995378645395834756983475693874598437");
    NT ld = largenum / largeden;
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_ld = mwt::rational_to_sqrt_nt(ld);
    CORE::BigInt ld_bn("29485739458039487510394019834570139485713094857130495831745134513");
    CORE::BigInt ld_bd("348573279458423429837429854279845995378645395834756983475693874598437");
    CORE::BigRat ld_br(ld_bn, ld_bd);
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT cgal_ld_br(ld_br);
    CHECK(cgal_ld == cgal_ld_br);
}
