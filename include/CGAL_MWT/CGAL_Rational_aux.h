#ifndef CGAL_MWT_CGAL_RATIONAL_AUX_H_INCLUDED_
#define CGAL_MWT_CGAL_RATIONAL_AUX_H_INCLUDED_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/number_utils.h>
#include <cstdint>
#include <sstream>

/**
 * Unfortunately, CGAL::Exact_rational and alike types do not have a lot of
 * well-defined special operations (such as checks for integrality),
 * because they are completely different types depending on how CGAL was compiled.
 * This header provides tools to detect the actual underlying type through
 * duck typing, i.e., without relying on the actual implementation headers to be present,
 * and to provide access to some of their internals.
 *
 * Furthermore, the CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT
 * type does not interoperate well with almost any other exact number type to
 * the point that creating one from another is not possible without resorting to
 * extreme measures such as string conversion which also carry the risk of precision loss.
 * In this header, we define (hopefully sensibly efficient) functions to make the types
 * interoperable.
 */

namespace mwt {

/**
 * Helper for static_assert (avoid hard errors,
 * but provide a reliably false value with a pseudo-dependency
 * on some type parameters T).
 */
template<class... T> constexpr bool always_false = false && (false && ... && std::is_void_v<T>);

namespace is_boost_number_detail {

/**
 * Check for boost::multiprecision::number-like type.
 * Basically, test for the presence of ::backend_type type,
 * .backend() member function, and operator+.
 */
template<typename NumType, typename = typename NumType::backend_type,
         typename = decltype(std::declval<NumType &>().backend()),
         typename = decltype(std::declval<NumType &>() + std::declval<const NumType &>())>
std::true_type is_boost_number_impl(const NumType *);
std::false_type is_boost_number_impl(...);

}

namespace is_gmp_backend_detail {

struct UniqueType {};

/**
 * Provide a fake mpq_init with recognizable return type
 * to check for the presence of the real, preferred mpq_init.
 */
template<typename X = UniqueType> static X mpq_init(...) {
    static_assert(always_false<X>, "This function template should never be instantiated");
    return {};
}

/**
 * Provide a fake mpz_init with recognizable return type
 * to check for the presence of the real, preferred mpq_init.
 */
template<typename X = UniqueType> static X mpz_init(...) {
    static_assert(always_false<X>, "This function template should never be instantiated");
    return {};
}

#ifndef mpq_numref
/**
 * Provide a fake mpq_numref with recognizable return type
 * to allow using the real mpq_numref without hard error from the
 * compiler even if the real gmp header is missing/not included by CGAL/Exact_rational.h.
 * Only defined if mpq_denref is not some macro.
 */
template<typename X = UniqueType> static X mpq_numref(...) {
    static_assert(always_false<X>, "This function template should never be instantiated");
    return {};
}
#endif

#ifndef mpq_denref
/**
 * Provide a fake mpq_denref with recognizable return type
 * to allow using the real mpq_numref without hard error from the
 * compiler even if the real gmp header is missing/not included by CGAL/Exact_rational.h.
 * Only defined if mpq_denref is not some macro.
 */
template<typename X = UniqueType> static X mpq_denref(...) {
    static_assert(always_false<X>, "This function template should never be instantiated");
    return {};
}
#endif

#ifndef mpz_cmp_si
/**
 * Provide a fake mpz_cmp_si
 * to allow using the real mpz_cmp_si without hard error from the
 * compiler even if the real gmp header is missing/not included by CGAL/Exact_rational.h.
 * Only defined if mpz_cmp_si is not some macro
 * (in which case the header is assumed to be present).
 */
template<typename A1> static A1 mpz_cmp_si(A1, ...) {
    static_assert(always_false<A1>, "This function template should never be instantiated");
    return {};
}
#endif

#ifndef mpz_fits_slong_p
/**
 * Provide a fake mpz_fits_slong_p
 * to allow using the real mpz_fits_slong_p without hard error from the
 * compiler even if the real gmp header is missing/not included by CGAL/Exact_rational.h.
 * Only defined if mpz_fits_slong_p is not some macro
 * (in which case the header is assumed to be present).
 */
template<typename X = UniqueType> static X mpz_fits_slong_p(...) {
    static_assert(always_false<X>, "This function template should never be instantiated");
    return {};
}
}
#endif

#ifndef mpz_get_si

/**
 * Provide a fake mpz_get_si
 * to allow using the real mpz_get_si without hard error from the
 * compiler even if the real gmp header is missing/not included by CGAL/Exact_rational.h.
 * Only defined if mpz_get_si is not some macro
 * (in which case the header is assumed to be present).
 */
template<typename X = UniqueType> static X mpz_get_si(...) {
    static_assert(always_false<X>, "This function template should never be instantiated");
    return {};
}

#endif

template<typename NumType> bool raw_mpq_rational_is_integer(const NumType &data) {
    auto denref = mpq_denref(data);
    static_assert(!std::is_convertible_v<decltype(denref), UniqueType>, "actual mpq_denref not found");
    return mpz_cmp_si(denref, long(1)) == 0;
}

template<typename NumType> std::optional<std::int64_t> raw_mpz_to_platform_int(const NumType &data) {
    if(!mpz_fits_slong_p(data))
        return std::nullopt;
    return std::int64_t(mpz_get_si(data));
}

template<typename NumType>
CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT raw_mpz_to_sqrt_nt(const NumType &data) {
    return CORE::Expr(CORE::BigInt{data});
}

template<typename NumType> std::optional<std::int64_t> raw_mpq_to_platform_int(const NumType &data) {
    if(!raw_mpq_rational_is_integer(data))
        return std::nullopt;
    auto numref = mpq_numref(data);
    static_assert(!std::is_convertible_v<decltype(numref), UniqueType>, "actual mpq_numref not found");
    return raw_mpz_to_platform_int(numref);
}

template<typename NumType>
CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT raw_mpq_to_sqrt_nt(const NumType &data) {
    return CORE::Expr(CORE::BigRat{data});
}

template<typename NumType, typename D = decltype(std::declval<NumType &>().data()),
         std::enable_if_t<!std::is_convertible_v<decltype(mpq_init(std::declval<D>())), UniqueType>, int> = 0>
std::true_type is_gmp_backend_impl(const NumType *);
std::false_type is_gmp_backend_impl(...);

template<typename NumType, typename D = decltype(std::declval<NumType &>().data()),
         std::enable_if_t<!std::is_convertible_v<decltype(mpz_init(std::declval<D>())), UniqueType>, int> = 0>
std::true_type is_gmpz_backend_impl(const NumType *);
std::false_type is_gmpz_backend_impl(...);

}

namespace is_mpq_class_v_detail {

/**
 * Check for gmpxx::mpq_class-like type.
 * Basically, test for the presence of the .get_mpq_t() member function
 * and operator+.
 */
template<typename NumType, typename = decltype(std::declval<NumType &>().get_mpq_t()),
         typename = decltype(std::declval<NumType &>() + std::declval<const NumType &>())>
std::true_type is_mpq_class_impl(const NumType *);
std::false_type is_mpq_class_impl(...);

template<typename NumType, typename = decltype(std::declval<NumType &>().get_mpz_t()),
         typename = decltype(std::declval<NumType &>() + std::declval<const NumType &>())>
std::true_type is_mpz_class_impl(const NumType *);
std::false_type is_mpz_class_impl(...);

}

namespace is_cpp_int_detail {

template<typename NumType,
         typename B = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<NumType &>().backend())>>,
         typename = typename NumType::backend_type, typename = typename B::checked_type>
std::true_type is_cpp_int_impl(const NumType *);
std::false_type is_cpp_int_impl(...);

}

namespace is_boost_rational_detail {

/**
 * Avoid hard errors from the compiler if no numerator(x)
 * function matching NumType exists, by providing a useless
 * and un-callable declaration of a numerator(x) function template.
 */
template<typename FakeNumeratorType, std::enable_if_t<std::is_void_v<FakeNumeratorType>, int> = 0>
void numerator(const FakeNumeratorType &);

/**
 * Avoid hard errors from the compiler if no denominator(x)
 * function matching NumType exists, by providing a useless
 * and un-callable declaration of a denominator(x) function template.
 */
template<typename FakeNumeratorType, std::enable_if_t<std::is_void_v<FakeNumeratorType>, int> = 0>
void denominator(const FakeNumeratorType &);

template<typename NumType,
         typename B = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<NumType &>().backend())>>,
         typename = typename NumType::backend_type,
         std::enable_if_t<!std::is_void_v<decltype(numerator(std::declval<NumType &>()))>, int> = 0,
         std::enable_if_t<!std::is_void_v<decltype(denominator(std::declval<NumType &>()))>, int> = 0>
std::true_type is_boost_rational_impl(const NumType *);
std::false_type is_boost_rational_impl(...);

}

namespace is_cgal_quotient_detail {

/**
 * Check for CGAL::Quotient-like type.
 * Basically, test for the presence of the .numerator() member function,
 * .denominator() member function, operator+, and ::NT type.
 */
template<typename NumType, typename = decltype(std::declval<NumType &>().numerator()),
         typename = decltype(std::declval<NumType &>().denominator()),
         typename = decltype(std::declval<NumType &>() + std::declval<const NumType &>()),
         typename = typename NumType::NT>
std::true_type is_cgal_quotient_impl(const NumType *);
std::false_type is_cgal_quotient_impl(...);

}

namespace has_exact_method_detail {

template<typename NumType, typename = decltype(std::declval<NumType &>().exact())>
std::true_type has_exact_method_impl(const NumType *);
std::false_type has_exact_method_impl(...);

}

/**
 * Should be true only for boost::multiprecision::number-like types,
 * which can have very different backends themselves; this also includes
 * cpp_int and cpp_rational, which are a special case of boost::multiprecision::number.
 */
template<typename NumType>
constexpr static bool is_boost_number_v =
    decltype(is_boost_number_detail::is_boost_number_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for boost::multiprecision::cpp_int-like types.
 */
template<typename NumType>
constexpr static bool is_boost_cpp_int_v =
    decltype(is_cpp_int_detail::is_cpp_int_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for boost::multiprecision::cpp_rational-like types,
 */
template<typename NumType>
constexpr static bool is_boost_rational_v =
    decltype(is_boost_rational_detail::is_boost_rational_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for boost::multiprecision::gmp_backend-like types,
 * which expose the underlying gmp data structure via a .data() method.
 */
template<typename NumType>
constexpr static bool is_boost_gmp_backend_v =
    decltype(is_gmp_backend_detail::is_gmp_backend_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for boost::multiprecision::gmp_int backend-like types,
 * which expose the underlying gmp data structure via a .data() method.
 */
template<typename NumType>
constexpr static bool is_boost_gmp_int_backend_v =
    decltype(is_gmp_backend_detail::is_gmpz_backend_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for gmpxx::mpq_class-like types.
 */
template<typename NumType>
constexpr static bool is_mpq_class_v =
    decltype(is_mpq_class_v_detail::is_mpq_class_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for gmpxx::mpz_class-like types.
 */
template<typename NumType>
constexpr static bool is_mpz_class_v =
    decltype(is_mpq_class_v_detail::is_mpz_class_impl(std::declval<const NumType *>()))::value;

/**
 * Should be true only for CGAL::Quotient-like types.
 */
template<typename NumType>
constexpr static bool is_cgal_quotient_v =
    decltype(is_cgal_quotient_detail::is_cgal_quotient_impl(std::declval<const NumType *>()))::value;

template<typename NumType>
constexpr static bool has_exact_method_v =
    decltype(has_exact_method_detail::has_exact_method_impl(std::declval<const NumType *>()))::value;

/**
 * Check if the value is actually zero;
 * may be more efficient than a == 0 comparison.
 */
template<typename NumType> inline bool rational_is_zero(const NumType &num) {
    if constexpr(is_boost_number_v<NumType>) {
        return num.is_zero();
    } else if constexpr(is_mpq_class_v<NumType>) {
        return num.get_num() == 0;
    } else {
        return num == 0;
    }
}

/**
 * Check if the given rational number is actually integer.
 */
template<typename NumType> inline bool rational_is_integer(const NumType &num) {
    if constexpr(is_boost_rational_v<NumType>) {
        return denominator(num) == 1;
    } else if constexpr(is_boost_cpp_int_v<NumType>) {
        static_assert(always_false<NumType>,
                      "rational_is_integer called on cpp_int, which is not a rational number type");
        return true;
    } else if constexpr(is_boost_number_v<NumType>) {
        return rational_is_integer(num.backend());
    } else if constexpr(is_boost_gmp_backend_v<NumType>) {
        auto data = num.data();
        return is_gmp_backend_detail::raw_mpq_rational_is_integer(data);
    } else if constexpr(is_mpq_class_v<NumType>) {
        return num.get_den() == 1;
    } else if constexpr(is_cgal_quotient_v<NumType>) {
        if(num.denominator() == 1 || num.numerator() == 0)
            return true;
        auto x = CGAL::gcd(num.numerator(), num.denominator());
        return x == num.denominator();
    } else {
        static_assert(always_false<NumType>, "rational_is_integer not implemented for this type");
    }
}

/**
 * As a fallback method, convert from the given rational or integer
 * to a CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT
 * via a string.
 */
template<typename NumType>
inline CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT
rational_to_sqrt_nt_via_string(const NumType &num) {
    std::ostringstream oss;
    oss << num;
    std::string result = oss.str();
    auto remove_char = [](char c) { return c == ' ' || c == '(' || c == ')'; };
    result.erase(std::remove_if(result.begin(), result.end(), remove_char), result.end());
    auto idx = result.find('/');
    if(idx != std::string::npos) {
        result[idx] = '\0';
        return CORE::BigRat(CORE::BigInt(result.c_str()), CORE::BigInt(result.c_str() + idx + 1));
    } else {
        return CORE::BigInt(result);
    }
}

/**
 * Depending on whether std::int64_t is long or long long,
 * we also need an extra method to convert std::int64_t into
 * a CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT.
 */
template<typename NumType, typename D = std::decay_t<NumType>,
         std::enable_if_t<std::is_same_v<D, std::int64_t>, int> = 0>
CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT rational_to_sqrt_nt(NumType x) {
    using Result = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT;
    if constexpr(std::is_convertible_v<D, Result>) {
        return x;
    } else {
        if(x >= std::numeric_limits<long>::min() && x <= std::numeric_limits<long>::max()) {
            return Result(long(x));
        }
        if(x < 0) {
            std::uint64_t val(x);
            val = -val;
            std::uint32_t v1(val >> 32);
            std::uint32_t v2(val);
            return Result(v1) * -4294967296.0 - v2;
        } else {
            std::uint64_t val(x);
            std::uint32_t v1(val >> 32);
            std::uint32_t v2(val);
            return Result(v1) * 4294967296.0 + v2;
        }
    }
}

template<typename NumType, std::enable_if_t<!std::is_integral_v<std::decay_t<NumType>>, int> = 0>
inline CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT rational_to_sqrt_nt(const NumType &num) {
    if constexpr(is_boost_cpp_int_v<NumType> || is_boost_rational_v<NumType>) {
        return rational_to_sqrt_nt_via_string(num);
    } else if constexpr(is_boost_number_v<NumType>) {
        return rational_to_sqrt_nt(num.backend());
    } else if constexpr(is_boost_gmp_backend_v<NumType>) {
        return is_gmp_backend_detail::raw_mpq_to_sqrt_nt(num.data());
    } else if constexpr(is_mpq_class_v<NumType>) {
        return CORE::BigRat(CORE::BigInt(num.get_num_mpz_t()), CORE::BigInt(num.get_den_mpz_t()));
    } else if constexpr(is_cgal_quotient_v<NumType>) {
        return rational_to_sqrt_nt_via_string(num);
    } else {
        static_assert(always_false<NumType>, "rational_to_sqrt_nt not implemented for this type");
    }
}

/**
 * Attempt to turn a big integer number type into a platform
 * integer value; if the number is too large or too small,
 * return std::nullopt.
 */
template<typename NumType> inline std::optional<std::int64_t> bigint_to_platform_int(const NumType &num) {
    if constexpr(is_mpz_class_v<NumType>) {
        if(num.fits_slong_p()) {
            return std::int64_t(num.get_si());
        }
        return std::nullopt;
    } else if constexpr(is_boost_cpp_int_v<NumType>) {
        if(num > std::numeric_limits<std::int64_t>::max() || num < std::numeric_limits<std::int64_t>::min()) {
            return std::nullopt;
        }
        return static_cast<std::int64_t>(num);
    } else if constexpr(is_boost_number_v<NumType>) {
        return bigint_to_platform_int(num.backend());
    } else if constexpr(is_boost_gmp_int_backend_v<NumType>) {
        return is_gmp_backend_detail::raw_mpz_to_platform_int(num.data());
    } else {
        static_assert(always_false<NumType>, "bigint_to_platform_int not implemented for this type");
    }
}

/**
 * Try to turn the rational number into a platform integer type.
 * If the number is not an integer or the conversion would overflow,
 * return std::nullopt.
 */
template<typename NumType> inline std::optional<std::int64_t> rational_to_platform_int(const NumType &num) {
    if constexpr(is_boost_rational_v<NumType>) {
        if(num.is_zero())
            return std::int64_t(0);
        if(denominator(num) != 1)
            return std::nullopt;
        return bigint_to_platform_int(numerator(num));
    } else if constexpr(is_boost_cpp_int_v<NumType>) {
        return bigint_to_platform_int(num);
    } else if constexpr(is_boost_number_v<NumType>) {
        return rational_to_platform_int(num.backend());
    } else if constexpr(is_boost_gmp_backend_v<NumType>) {
        return is_gmp_backend_detail::raw_mpq_to_platform_int(num.data());
    } else if constexpr(is_mpq_class_v<NumType>) {
        if(num.get_den() != 1)
            return std::nullopt;
        return bigint_to_platform_int(num.get_num());
    } else if constexpr(is_cgal_quotient_v<NumType>) {
        if(num.numerator() == 0)
            return std::int64_t(0);
        if(num.denominator() == 1) {
            return bigint_to_platform_int(num.numerator());
        }
        auto x = CGAL::gcd(num.numerator(), num.denominator());
        if(x != num.denominator())
            return std::nullopt;
        x = num.numerator() / x;
        return bigint_to_platform_int(x);
    } else {
        static_assert(always_false<NumType>, "rational_to_platform_int not implemented for this type");
    }
}

/**
 * Check if CGAL considers the given type to be
 * a RealEmbeddable type or a possibly cv-qualified reference to it.
 */
template<typename CGALNumberType> static constexpr bool is_real_embeddable() {
    using D = std::decay_t<CGALNumberType>;
    using Traits = CGAL::Real_embeddable_traits<D>;
    return std::is_same_v<typename Traits::Is_real_embeddable, CGAL::Tag_true>;
}

static_assert(is_real_embeddable<CGAL::Exact_rational>());
static_assert(is_real_embeddable<std::uint64_t>());
static_assert(is_real_embeddable<std::int64_t>());
static_assert(is_real_embeddable<double>());
static_assert(is_real_embeddable<CGAL::Epeck::FT>());
static_assert(is_real_embeddable<CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT>());

/**
 * Provide to_interval for number types for which CGAL provides it.
 */
template<typename CGALNumberType, std::enable_if_t<is_real_embeddable<CGALNumberType>(), int> = 0>
CGAL::Interval_nt<false> to_interval(const CGALNumberType &num) {
    return CGAL::to_interval(num);
}

/**
 * Provide to_interval for number types for which CGAL does not
 * provide it, but which have a .to_interval() method.
 * The main candidate for this is our RationalOrInt type.
 */
template<typename CGALNumberType, std::enable_if_t<!is_real_embeddable<CGALNumberType>(), int> = 0,
         typename D = decltype(std::declval<const CGALNumberType &>().to_interval()),
         std::enable_if_t<std::is_convertible_v<D, CGAL::Interval_nt_advanced>, int> = 0>
CGAL::Interval_nt<false> to_interval(const CGALNumberType &num) {
    return num.to_interval();
}

/**
 * exact_to_interval is an auxiliary function to convert a number
 * to an interval, but invokes .exact() on lazy number types before.
 */

/**
 * Provide exact_to_interval for types with .exact() and mwt::to_interval.
 */
template<typename NumberType, std::enable_if_t<has_exact_method_v<NumberType>, int> = 0,
         typename = decltype(mwt::to_interval(std::declval<const NumberType &>()))>
CGAL::Interval_nt<false> exact_to_interval(const NumberType &num) {
    num.exact();
    return mwt::to_interval(num);
}

/**
 * Provide exact_to_interval for types without .exact() but with mwt::to_interval.
 */
template<typename NumberType, std::enable_if_t<!has_exact_method_v<NumberType>, int> = 0,
         typename = decltype(mwt::to_interval(std::declval<const NumberType &>()))>
CGAL::Interval_nt<false> exact_to_interval(const NumberType &num) {
    return mwt::to_interval(num);
}
}

#endif
