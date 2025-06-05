#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>
#include <CGAL_MWT/Rational_or_int.h>
#include <doctest/doctest.h>

TEST_CASE("[RationalOrInt] test overflow detecting addition") {
    std::int64_t large = std::numeric_limits<std::int64_t>::max() - 5;
    std::int64_t small = 5;
    std::int64_t mlarge = std::numeric_limits<std::int64_t>::min() + 5;
    std::int64_t msmall = -5;

    auto does_not_overflow = [](std::int64_t lhs, auto &&rhs, std::optional<std::int64_t> expected = {}) -> bool {
        std::int64_t result;
        bool overflow = mwt::overflow_check_detail::add_to_i64(lhs, rhs, result);
        CHECK(!overflow);
        if(!overflow) {
            if(expected) {
                CHECK(result == *expected);
            } else {
                CHECK(result == lhs + std::int64_t(rhs));
            }
        }
        return !overflow;
    };

    auto does_overflow = [](std::int64_t lhs, auto &&rhs) -> bool {
        std::int64_t result;
        bool overflow = mwt::overflow_check_detail::add_to_i64(lhs, rhs, result);
        CHECK(overflow);
        return overflow;
    };

    CHECK(does_not_overflow(large, small));
    CHECK(does_not_overflow(small, large));
    CHECK(does_not_overflow(mlarge, msmall));
    CHECK(does_not_overflow(msmall, mlarge));
    CHECK(does_overflow(large, std::numeric_limits<std::int64_t>::max()));
    CHECK(does_overflow(std::numeric_limits<std::int64_t>::max(), large));
    CHECK(does_overflow(mlarge, std::numeric_limits<std::int64_t>::min()));
    CHECK(does_overflow(std::numeric_limits<std::int64_t>::min(), mlarge));
    CHECK(does_overflow(mlarge, -6));
    CHECK(does_overflow(-6, mlarge));
    CHECK(does_overflow(large, 6));
    CHECK(does_overflow(6, large));
    CHECK(does_overflow(-1, std::numeric_limits<std::int64_t>::min()));
    CHECK(does_not_overflow(std::numeric_limits<std::int64_t>::min(), std::numeric_limits<std::uint64_t>::max(),
                            std::numeric_limits<std::int64_t>::max()));
    CHECK(does_overflow(5, std::uint64_t(std::numeric_limits<std::int64_t>::max() - 4)));
    CHECK(does_not_overflow(5, std::uint64_t(std::numeric_limits<std::int64_t>::max() - 5)));
}

TEST_CASE_TEMPLATE("[RationalOrInt] construction", SourceType, char, unsigned char, short, unsigned short, int,
                   unsigned int, long, unsigned long, long long, unsigned long long, CGAL::Exact_rational, double,
                   float, long double) {
    mwt::RationalOrInt<CGAL::Exact_rational> vneg(-std::numeric_limits<std::int64_t>::max());
    CHECK(vneg.is_platform_int());
    CHECK(vneg.get_platform_int() == -std::numeric_limits<std::int64_t>::max());
    CGAL::Exact_predicates_exact_constructions_kernel::FT vneg_ft(-9223372036854775808.0);
    vneg_ft += 1.0;
    CHECK(vneg.cast_to<CGAL::Exact_predicates_exact_constructions_kernel::FT>() == vneg_ft);

    SourceType value{5};
    mwt::RationalOrInt<CGAL::Exact_rational> rational_value(value);
    CHECK(!rational_value.is_zero());
    CHECK(rational_value.is_platform_int());
    CHECK(rational_value.get_platform_int() == 5);
    CHECK(rational_value.unchecked_platform_int() == 5);
    CHECK(rational_value.as_rational() == 5);

    mwt::RationalOrInt<CGAL::Exact_rational> zero_value;
    CHECK(zero_value.is_zero());
    CHECK(zero_value.is_platform_int());
    CHECK(zero_value.get_platform_int() == 0);

    mwt::RationalOrInt<CGAL::Exact_rational> copied(zero_value);
    copied = rational_value;
    CHECK(!copied.is_zero());
    CHECK(copied.is_platform_int());
    CHECK(copied.get_platform_int() == 5);

    mwt::RationalOrInt<CGAL::Exact_rational> moved(std::move(copied));
    CHECK(moved.is_platform_int());
    CHECK(moved.get_platform_int() == 5);

    mwt::RationalOrInt<CGAL::Exact_rational> move_to;
    move_to = std::move(moved);
    CHECK(move_to.is_platform_int());
    CHECK(move_to.get_platform_int() == 5);

    CHECK(!rational_value.is_zero());
    CHECK(rational_value.is_platform_int());
    CHECK(rational_value.get_platform_int() == 5);
    CHECK(rational_value.unchecked_platform_int() == 5);
    CHECK(rational_value.as_rational() == 5);

    if constexpr(std::is_integral_v<SourceType>) {
        mwt::RationalOrInt<CGAL::Exact_rational> pair_constructed{SourceType{2}, SourceType{3}};
        CHECK(!pair_constructed.is_zero());
        CHECK(!pair_constructed.is_platform_int());
        CHECK(pair_constructed.get_rational() == CGAL::Exact_rational(2, 3));
    }
}

TEST_CASE_TEMPLATE("[RationalOrInt] implicit conversion", SourceType, char, unsigned char, short, unsigned short, int,
                   unsigned int, long, unsigned long, long long, unsigned long long, CGAL::Exact_rational, double,
                   float, long double) {
    SourceType value{10};
    auto implicit_convert = [](SourceType value) -> mwt::RationalOrInt<CGAL::Exact_rational> { return value; };
    auto rational = implicit_convert(value);
    CHECK(!rational.is_zero());
    CHECK(rational.is_platform_int());
    CHECK(rational.get_platform_int() == 10);
    CHECK(rational.unchecked_platform_int() == 10);
}

TEST_CASE("[RationalOrInt] long double handling") {
    if constexpr (sizeof(long double) != sizeof(double)) {
        long double x = 5.0000000000000000004336808689942018L;
        std::cout << sizeof(x) << std::endl;
        std::cout << x << std::endl;
        mwt::RationalOrInt<CGAL::Exact_rational> rational_value(x);
        CHECK(double(x) == 5);
        CHECK(!rational_value.is_zero());
        CHECK(!rational_value.is_platform_int());
    } else {
        volatile long double x = 5.0000000000000000004336808689942018L;
        volatile long double y = 5.0L;
        CHECK(x == y);
        mwt::RationalOrInt<CGAL::Exact_rational> rational_value(x);
        CHECK(double(x) == 5);
        CHECK(!rational_value.is_zero());
        CHECK(rational_value.is_platform_int());
    }
}

TEST_CASE("[RationalOrInt] addition operator+=") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    auto x = val1;
    x += val2;
    CHECK(x.get_rational() == CGAL::Exact_rational(28, 5));
    x += val3;
    CHECK(x.as_rational() == CGAL::Exact_rational(30, 5));
    CHECK(x.get_platform_int() == 6);

    auto y = val1;
    y += val1;
    CHECK(y.get_platform_int() == 10);

    auto z = val1;
    z += 5;
    CHECK(z.is_platform_int());
    CHECK(z.get_platform_int() == 10);
    z += 7.3;
    CHECK(!z.is_platform_int());
    CHECK(z.as_rational() == CGAL::Exact_rational(10) + 7.3);
}

TEST_CASE("[RationalOrInt] addition operator+") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    CHECK((val2 + val3).is_platform_int());
    CHECK((val2 + val3).get_platform_int() == 1);
    CHECK((val1 + val2 + val3).get_platform_int() == 6);
    CHECK((val1 + val1).get_platform_int() == 10);
    CHECK(val1.get_platform_int() == 5);

    char c = 1;
    unsigned char uc = 2;
    short s = 3;
    unsigned short us = 4;
    int i = 5;
    unsigned int ui = 6;
    long l = 7;
    unsigned long ul = 8;
    long long ll = 9;
    unsigned long long ull = 10;
    double d = 11.1;
    float f = 12.0f;
    long double ld = 13.0L;

    CHECK((val1 + c).get_platform_int() == 6);
    CHECK((val1 + uc).get_platform_int() == 7);
    CHECK((val1 + s).get_platform_int() == 8);
    CHECK((val1 + us).get_platform_int() == 9);
    CHECK((val1 + i).get_platform_int() == 10);
    CHECK((val1 + ui).get_platform_int() == 11);
    CHECK((val1 + l).get_platform_int() == 12);
    CHECK((val1 + ul).get_platform_int() == 13);
    CHECK((val1 + ll).get_platform_int() == 14);
    CHECK((val1 + ull).get_platform_int() == 15);
    CHECK((val1 + d).get_rational() == CGAL::Exact_rational(5) + 11.1);
    CHECK((val1 + f).get_platform_int() == 17);
    CHECK((val1 + ld).get_platform_int() == 18);

    CHECK((c + val1).get_platform_int() == 6);
    CHECK((uc + val1).get_platform_int() == 7);
    CHECK((s + val1).get_platform_int() == 8);
    CHECK((us + val1).get_platform_int() == 9);
    CHECK((i + val1).get_platform_int() == 10);
    CHECK((ui + val1).get_platform_int() == 11);
    CHECK((l + val1).get_platform_int() == 12);
    CHECK((ul + val1).get_platform_int() == 13);
    CHECK((ll + val1).get_platform_int() == 14);
    CHECK((ull + val1).get_platform_int() == 15);
    CHECK((d + val1).get_rational() == CGAL::Exact_rational(5) + 11.1);
    CHECK((f + val1).get_platform_int() == 17);
    CHECK((ld + val1).get_platform_int() == 18);
}

TEST_CASE("[RationalOrInt] subtraction operator-=") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    auto x = val1;
    x -= val2;
    CHECK(x.get_rational() == CGAL::Exact_rational(22, 5));
    x -= val3;
    CHECK(x.as_rational() == 4);
    CHECK(x.get_platform_int() == 4);

    auto y = val1;
    y -= val1;
    CHECK(y.get_platform_int() == 0);
    CHECK(y.is_zero());

    auto z = val1;
    z -= 5;
    CHECK(z.is_platform_int());
    CHECK(z.get_platform_int() == 0);
    CHECK(z.is_zero());
    z -= 7.3;
    CHECK(!z.is_platform_int());
    CHECK(z.as_rational() == -7.3);
}

TEST_CASE("[RationalOrInt] subtraction operator-") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    CHECK(!(val2 - val3).is_platform_int());
    CHECK((val2 - val3).get_rational() == CGAL::Exact_rational(1, 5));
    CHECK((val1 - val2 - val3).get_platform_int() == 4);
    CHECK((val1 - val1).get_platform_int() == 0);
    CHECK(val1.get_platform_int() == 5);

    char c = 1;
    unsigned char uc = 2;
    short s = 3;
    unsigned short us = 4;
    int i = 5;
    unsigned int ui = 6;
    long l = 7;
    unsigned long ul = 8;
    long long ll = 9;
    unsigned long long ull = 10;
    double d = 11.1;
    float f = 12.0f;
    long double ld = 13.0L;

    CHECK((val1 - c).get_platform_int() == 4);
    CHECK((val1 - uc).get_platform_int() == 3);
    CHECK((val1 - s).get_platform_int() == 2);
    CHECK((val1 - us).get_platform_int() == 1);
    CHECK((val1 - i).get_platform_int() == 0);
    CHECK((val1 - ui).get_platform_int() == -1);
    CHECK((val1 - l).get_platform_int() == -2);
    CHECK((val1 - ul).get_platform_int() == -3);
    CHECK((val1 - ll).get_platform_int() == -4);
    CHECK((val1 - ull).get_platform_int() == -5);
    CHECK((val1 - d).get_rational() == CGAL::Exact_rational(5) - 11.1);
    CHECK((val1 - f).get_platform_int() == -7);
    CHECK((val1 - ld).get_platform_int() == -8);

    CHECK((c - val1).get_platform_int() == -4);
    CHECK((uc - val1).get_platform_int() == -3);
    CHECK((s - val1).get_platform_int() == -2);
    CHECK((us - val1).get_platform_int() == -1);
    CHECK((i - val1).get_platform_int() == 0);
    CHECK((ui - val1).get_platform_int() == 1);
    CHECK((l - val1).get_platform_int() == 2);
    CHECK((ul - val1).get_platform_int() == 3);
    CHECK((ll - val1).get_platform_int() == 4);
    CHECK((ull - val1).get_platform_int() == 5);
    CHECK((d - val1).get_rational() == 11.1 - CGAL::Exact_rational(5));
    CHECK((f - val1).get_platform_int() == 7);
    CHECK((ld - val1).get_platform_int() == 8);
}

TEST_CASE("[RationalOrInt] multiplication operator*=") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    auto x = val1;
    x *= val2;
    CHECK(x.get_platform_int() == 3);

    auto y = val1;
    y *= val1;
    y *= val3;
    CHECK(y.get_platform_int() == 10);

    auto z = val1;
    z *= 5;
    z *= val2;
    CHECK(z.as_rational() == 15);
    CHECK(z.get_platform_int() == 15);
    z *= 7.3;
    CHECK(!z.is_platform_int());
    CHECK(z.get_rational() == CGAL::Exact_rational(15) * 7.3);
}

TEST_CASE("[RationalOrInt] multiplication operator*") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    CHECK((val1 * val2).is_platform_int());
    CHECK(!(val2 * val3).is_platform_int());
    CHECK((val2 * val3).get_rational() == CGAL::Exact_rational(6, 25));
    CHECK((val1 * val2 * val3 * val1).get_platform_int() == 6);
    CHECK((val1 * val1).get_platform_int() == 25);
    CHECK(val1.get_platform_int() == 5);

    char c = 1;
    unsigned char uc = 2;
    short s = 3;
    unsigned short us = 4;
    int i = 5;
    unsigned int ui = 6;
    long l = 7;
    unsigned long ul = 8;
    long long ll = 9;
    unsigned long long ull = 10;
    double d = 11.1;
    float f = 12.0f;
    long double ld = 13.0L;

    CHECK(val1 * c == 5);
    CHECK(val1 * uc == 10);
    CHECK(val1 * s == 15);
    CHECK(val1 * us == 20);
    CHECK(val1 * i == 25);
    CHECK(val1 * ui == 30);
    CHECK(val1 * l == 35);
    CHECK(val1 * ul == 40);
    CHECK(val1 * ll == 45);
    CHECK(val1 * ull == 50);
    CHECK(val1 * d == 5 * CGAL::Exact_rational(11.1));
    CHECK(val1 * f == 60);
    CHECK(val1 * ld == 65);

    CHECK(c * val1 == 5);
    CHECK(uc * val1 == 10);
    CHECK(s * val1 == 15);
    CHECK(us * val1 == 20);
    CHECK(i * val1 == 25);
    CHECK(ui * val1 == 30);
    CHECK(l * val1 == 35);
    CHECK(ul * val1 == 40);
    CHECK(ll * val1 == 45);
    CHECK(ull * val1 == 50);
    CHECK(d * val1 == 5 * CGAL::Exact_rational(11.1));
    CHECK(f * val1 == 60);
    CHECK(ld * val1 == 65);
}

TEST_CASE("[RationalOrInt] division") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(168);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 7);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(21);

    auto x = val1;
    x /= val2;
    CHECK(x.get_platform_int() == 392);
    x /= val3;
    CHECK(x.get_rational() == CGAL::Exact_rational(392, 21));
    CHECK(!x.is_platform_int());
    CHECK(!!x);
    CHECK(val1 == 168);
    CHECK(val2 / val2 == 1);
    CHECK((val2 / val2).is_platform_int());
    auto z = val1;
    z /= val3;
    CHECK(z.is_platform_int());
    CHECK(z.get_platform_int() == 8);
}

TEST_CASE("[RationalOrInt] comparison operators") {
    mwt::RationalOrInt<CGAL::Exact_rational> val1(5);
    mwt::RationalOrInt<CGAL::Exact_rational> val2(3, 5);
    mwt::RationalOrInt<CGAL::Exact_rational> val3(2, 5);

    CHECK(val1 == val1);
    CHECK(!(val1 != val1));
    CHECK(val1 != val2);
    CHECK(val1 != val3);
    CHECK(val2 != val3);
    CHECK(val2 < val1);
    CHECK(val2 <= val1);
    CHECK(val1 <= val1);
    CHECK(val2 > val3);
    CHECK(val2 >= val3);
    CHECK(val2 >= val2);

    CHECK(val1 == 5);
    CHECK(val1 != 3);
    CHECK(val1 != 2);
    CHECK(val2 != 2);
    CHECK(3 < val1);
    CHECK(3 <= val1);
    CHECK(val1 <= 6);
    CHECK(val2 > -2);
    CHECK(val2 >= -2);
    CHECK(val3 >= -3);
    CHECK(-val1 < std::uint64_t(3));

    char c = 5;
    unsigned char uc = 2;
    short s = 3;
    unsigned short us = 4;
    int i = 5;
    unsigned int ui = 6;
    long l = 7;
    unsigned long ul = 8;
    long long ll = 9;
    unsigned long long ull = 10;
    double d = 11.1;
    float f = 12.0f;
    long double ld = 13.0L;

    CHECK(val1 == c);
    CHECK(val1 != uc);
    CHECK(val1 != s);
    CHECK(val2 != us);
    CHECK(val1 <= i);
    CHECK(val1 >= i);
    CHECK(val1 < l);
    CHECK(val2 > -s);
    CHECK(val2 <= s);
    CHECK(val2 >= -i);
    CHECK(val3 < i);
    CHECK(val1 < ul);
    CHECK(val1 != ll);
    CHECK(val2 != ull);
    CHECK(val1 < d);
    CHECK(val1 + val1 + val1 > f);
    CHECK(val1 + val1 + val1 >= ld);

    CHECK(c == val1);
    CHECK(uc != val1);
    CHECK(s != val1);
    CHECK(us != val2);
    CHECK(i >= val1);
    CHECK(i <= val1);
    CHECK(l > val1);
    CHECK(-s < val2);
    CHECK(s >= val2);
    CHECK(-i <= val2);
    CHECK(i > val3);
    CHECK(ul > val1);
    CHECK(ll != val1);
    CHECK(ull != val2);
    CHECK(d > val1);
    CHECK(f < val1 + val1 + val1);
    CHECK(ld <= val1 + val1 + val1);
}
