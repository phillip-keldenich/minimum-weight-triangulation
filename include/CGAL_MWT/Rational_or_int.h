#ifndef CGAL_MWT_RATIONAL_OR_INT_H_INCLUDED_
#define CGAL_MWT_RATIONAL_OR_INT_H_INCLUDED_

#include "CGAL_Rational_aux.h"
#include <CGAL/CORE_BigRat.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <variant>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace mwt {

namespace overflow_check_detail {

/**
 * Overflow-checking addition of two signed 64-bit integers.
 * Called if there is no type-generic builtin available.
 */
static inline bool add_s64s64_overflow(std::int64_t lhs, std::int64_t rhs, std::int64_t &result) {
#ifdef _MSC_VER
    return _add_overflow_i64(0, lhs, rhs, &result) != 0;
#else
    if(lhs >= 0) {
        if(std::numeric_limits<std::int64_t>::max() - lhs < rhs) {
            return true;
        }
    } else {
        if(rhs < std::numeric_limits<std::int64_t>::min() - lhs) {
            return true;
        }
    }
    result = lhs + rhs;
    return false;
#endif
}

/**
 * Overflow-checking subtraction of two signed 64-bit integers.
 * Called if there is no type-generic builtin available.
 */
static inline bool sub_s64s64_overflow(std::int64_t lhs, std::int64_t rhs, std::int64_t &result) {
#ifdef _MSC_VER
    return _sub_overflow_i64(0, lhs, rhs, &result) != 0;
#else
    if(rhs == std::numeric_limits<std::int64_t>::min()) {
        if(lhs >= 0) {
            return true;
        } else {
            result = lhs - rhs;
            return false;
        }
    }
    return add_s64s64_overflow(lhs, -rhs, result);
#endif
}

/**
 * Overflow-checking multiplication of two signed 64-bit integers.
 * Called if there is no type-generic builtin available.
 */
template<typename T = void>
static inline bool mul_s64s64_overflow(std::int64_t lhs, std::int64_t rhs, std::int64_t &result) {
#ifdef _MSC_VER
    return _mul_overflow_i64(lhs, rhs, &result) != 0;
#else
    static_assert(std::is_void_v<T> && !std::is_void_v<T>, "mul_s64s64_overflow not implemented without intrinsics");
#endif
}

/**
 * Overflow-checking addition of some number type to a 64-bit signed integer.
 * Uses a type-generic builtin if available.
 */
template<typename IntType> bool add_to_i64(std::int64_t value, IntType add, std::int64_t &result) {
#if (defined(__GNUC__) && __GNUC__ >= 5) || defined(__clang__)
    // if possible, use builtin to make this really efficient (add + branch)
    return __builtin_add_overflow(value, add, &result);
#else
    if constexpr(std::is_unsigned_v<IntType>) {
        // IntType is an unsigned integer type
        constexpr std::uint64_t umax(std::numeric_limits<std::int64_t>::max());
        if constexpr(umax >= std::numeric_limits<IntType>::max()) {
            // can reduce this to signed 64-bit addition in any case
            return add_s64s64_overflow(value, std::int64_t(add), result);
        } else {
            if(add > umax) {
                if(value >= 0 || -std::uint64_t(value) < add - umax) {
                    return true;
                } else {
                    result = std::int64_t(IntType(value) + add);
                    return false;
                }
            } else {
                return add_s64s64_overflow(value, std::int64_t(add), result);
            }
        }
    } else {
        // IntType is a signed integer type
        if constexpr(std::numeric_limits<std::int64_t>::max() >= std::numeric_limits<IntType>::max()) {
            return add_s64s64_overflow(value, std::int64_t(add), result);
        } else {
            if(value >= 0) {
                if(std::numeric_limits<std::int64_t>::max() - value < add) {
                    return true;
                }
            } else {
                if(add < std::numeric_limits<std::int64_t>::min() - value) {
                    return true;
                }
            }
            result = value + add;
            return false;
        }
    }
#endif
}

template<typename IntType> bool sub_from_i64(std::int64_t value, IntType sub, std::int64_t &result) {
#if (defined(__GNUC__) && __GNUC__ >= 5) || defined(__clang__)
    // if possible, use builtin to make this really efficient (add + branch)
    return __builtin_sub_overflow(value, sub, &result);
#else
    if constexpr(std::is_unsigned_v<IntType>) {
        // IntType is an unsigned integer type
        constexpr std::uint64_t umax(std::numeric_limits<std::int64_t>::max());
        if constexpr(umax >= std::numeric_limits<IntType>::max()) {
            // can reduce this to signed 64-bit addition in any case
            return sub_s64s64_overflow(value, std::int64_t(sub), result);
        } else {
            if(sub > umax) {
                if(value < 0 || sub - 1 > umax + std::uint64_t(value)) {
                    return true;
                } else {
                    std::int64_t v = std::int64_t(sub - umax - 1);
                    v = value - v;
                    v -= std::numeric_limits<std::int64_t>::max();
                    result = v - 1;
                    return false;
                }
            } else {
                return sub_s64s64_overflow(value, std::int64_t(sub), result);
            }
        }
    } else {
        // IntType is a signed integer type
        if constexpr(std::numeric_limits<std::int64_t>::max() >= std::numeric_limits<IntType>::max()) {
            return sub_s64s64_overflow(value, std::int64_t(sub), result);
        } else {
            static_assert(std::numeric_limits<std::int64_t>::max() >= std::numeric_limits<IntType>::max(),
                          "sub_from_i64 not fully implemented without intrinsics");
        }
    }
#endif
}

template<typename IntType> bool mul_to_i64(std::int64_t value, IntType mul, std::int64_t &result) {
#if (defined(__GNUC__) && __GNUC__ >= 5) || defined(__clang__)
    return __builtin_mul_overflow(value, mul, &result);
#else
    if constexpr(std::is_unsigned_v<IntType>) {
        // IntType is an unsigned integer type
        constexpr std::uint64_t umax(std::numeric_limits<std::int64_t>::max());
        if constexpr(umax >= std::numeric_limits<IntType>::max()) {
            // can reduce this to signed 64-bit multiplication in any case
            return mul_s64s64_overflow(value, std::int64_t(mul), result);
        } else {
            if(mul > umax) {
                if(value == 0) {
                    result = 0;
                    return false;
                }
                if(value == -1 && mul == umax + 1) {
                    result = std::numeric_limits<std::int64_t>::min();
                    return false;
                }
                return true;
            } else {
                return mul_s64s64_overflow(value, std::int64_t(mul), result);
            }
        }
    } else {
        // IntType is a signed integer type
        if constexpr(std::numeric_limits<std::int64_t>::max() >= std::numeric_limits<IntType>::max()) {
            return mul_s64s64_overflow(value, std::int64_t(mul), result);
        } else {
            static_assert(std::numeric_limits<std::int64_t>::max() >= std::numeric_limits<IntType>::max(),
                          "mul_to_i64 not fully implemented without intrinsics");
        }
    }
#endif
}

template<typename IntType> bool div_i64_is_exact(std::int64_t value, IntType div, std::int64_t &result) {
    if(div == 0) {
        throw std::logic_error("division by zero");
    }
    if constexpr(std::is_signed_v<IntType>) {
        if(value == std::numeric_limits<std::int64_t>::min() && div == -1) {
            return false;
        }
    }
    if(value % div != 0) {
        return false;
    }
    result = value / div;
    return true;
}

template<typename IntType> bool conversion_to_i64(IntType value, std::int64_t &result) {
    if constexpr(std::is_signed_v<IntType>) {
        if constexpr(std::numeric_limits<std::int64_t>::min() > std::numeric_limits<IntType>::min() ||
                     std::numeric_limits<std::int64_t>::max() < std::numeric_limits<IntType>::max()) {
            if(value > std::numeric_limits<std::int64_t>::max() || value < std::numeric_limits<std::int64_t>::min()) {
                return false;
            }
        }
    } else {
        constexpr std::uint64_t umax(std::numeric_limits<std::int64_t>::max());
        if constexpr(umax < std::numeric_limits<IntType>::max()) {
            if(value > umax) {
                return false;
            }
        }
    }
    result = std::int64_t(value);
    return true;
}

} // namespace overflow_check_detail

template<typename TargetIntType, typename SourceIntType,
         std::enable_if_t<!std::is_integral_v<TargetIntType> || !std::is_integral_v<SourceIntType>, int> = 0>
constexpr static inline bool nooverflow_convertible() {
    return false;
}

template<typename TargetIntType, typename SourceIntType, std::enable_if_t<std::is_integral_v<TargetIntType>, int> = 0,
         std::enable_if_t<std::is_integral_v<SourceIntType>, int> = 0>
constexpr static inline bool nooverflow_convertible() {
    constexpr bool target_signed = std::is_signed_v<TargetIntType>;
    constexpr bool source_signed = std::is_signed_v<SourceIntType>;
    constexpr bool same_sign = (target_signed == source_signed);
    if constexpr(same_sign) {
        return (std::numeric_limits<TargetIntType>::max() >= std::numeric_limits<SourceIntType>::max() &&
                std::numeric_limits<TargetIntType>::min() <= std::numeric_limits<SourceIntType>::min());
    } else if(source_signed) {
        return false;
    } else {
        using UnsignedTarget = std::make_unsigned_t<TargetIntType>;
        constexpr UnsignedTarget target_max(std::numeric_limits<TargetIntType>::max());
        return target_max >= std::numeric_limits<SourceIntType>::max();
    }
}

template<class... Ts> struct multi_visitor : Ts... {
    using Ts::operator()...;
};
template<class... Ts> multi_visitor(Ts...) -> multi_visitor<Ts...>;

/**
 * Arbitrary precision rational with
 * heap-storage-free special-case handling
 * for values that fit into platform integers.
 */
template<typename ArbitraryPrecisionRationalType> class RationalOrInt {
  public:
    using IntType = std::int64_t;
    using RationalType = ArbitraryPrecisionRationalType;
    using Raw = std::variant<IntType, RationalType>;

    /**
     * RationalOrInt is implicitly noexcept-constructible from
     * any integer type that cannot possibly overflow a std::int64_t.
     */
    template<typename OtherIntType, std::enable_if_t<nooverflow_convertible<IntType, OtherIntType>(), int> = 0>
    /*implicit*/ RationalOrInt(OtherIntType value) noexcept : m_value(IntType(value)) {}

    /**
     * RationalOrInt is implicitly constructible from any builtin integer type.
     */
    template<typename OtherIntType, std::enable_if_t<!nooverflow_convertible<IntType, OtherIntType>(), int> = 0,
             std::enable_if_t<std::is_integral_v<OtherIntType>, int> = 0>
    /*implicit*/ RationalOrInt(OtherIntType value) {
        std::int64_t result;
        if(overflow_check_detail::conversion_to_i64(value, result)) {
            m_value = IntType(result);
        } else {
            m_value = RationalType(value);
        }
    }

    /**
     * RationalOrInt is implicitly constructible from any type that is (implicitly) convertible to RationalType.
     */
    template<typename OtherNumType, typename D = std::decay_t<OtherNumType>,
             std::enable_if_t<!nooverflow_convertible<IntType, D>(), int> = 0,
             std::enable_if_t<!std::is_integral_v<D>, int> = 0,
             std::enable_if_t<std::is_constructible_v<RationalType, D>, int> = 0,
             std::enable_if_t<std::is_convertible_v<D, RationalType>, int> = 0>
    /*implicit*/ RationalOrInt(OtherNumType &&value)
        : m_value(std::in_place_type<RationalType>, std::forward<OtherNumType>(value)) {
        p_to_int_if_possible();
    }

    /**
     * RationalOrInt is explicitly constructible from any type that is explicitly convertible to RationalType.
     */
    template<typename OtherNumType, typename D = std::decay_t<OtherNumType>,
             std::enable_if_t<!nooverflow_convertible<IntType, D>(), int> = 0,
             std::enable_if_t<!std::is_integral_v<D>, int> = 0,
             std::enable_if_t<std::is_constructible_v<RationalType, D>, int> = 0,
             std::enable_if_t<!std::is_convertible_v<D, RationalType>, int> = 0>
    explicit RationalOrInt(OtherNumType &&value)
        : m_value(std::in_place_type<RationalType>, std::forward<OtherNumType>(value)) {
        p_to_int_if_possible();
    }

    /**
     * RationalOrInt can be constructed from two integers.
     */
    template<typename OtherIntType1, typename OtherIntType2,
             std::enable_if_t<std::is_integral_v<OtherIntType1>, int> = 0,
             std::enable_if_t<std::is_integral_v<OtherIntType2>, int> = 0>
    RationalOrInt(OtherIntType1 num, OtherIntType2 den) {
        if(den == 0) {
            throw std::logic_error("division by zero");
        } else if(num == 0) {
            m_value = IntType(0);
        } else if(den == -1 && num == std::numeric_limits<OtherIntType1>::min()) {
            m_value.template emplace<RationalType>(num, den);
            p_to_int_if_possible();
        } else if(den == 1) {
            p_convert_or_rational(num);
        } else if(den == -1) {
            p_convert_or_rational(-num);
        } else if(num % den == 0) {
            p_convert_or_rational(num / den);
        } else {
            m_value.template emplace<RationalType>(num, den);
        }
    }

    RationalOrInt() noexcept : m_value(IntType(0)) {}
    RationalOrInt(const RationalOrInt &) = default;
    RationalOrInt(RationalOrInt &&) noexcept = default;
    RationalOrInt &operator=(const RationalOrInt &) = default;
    RationalOrInt &operator=(RationalOrInt &&) noexcept = default;
    ~RationalOrInt() = default;

    /**
     * Check if *this currently holds
     * a platform integer value.
     */
    bool is_platform_int() const noexcept { return std::holds_alternative<IntType>(m_value); }

    /**
     * Get the value as platform int (throws an
     * exception if it doesn't hold a platform int).
     */
    IntType get_platform_int() const { return std::get<IntType>(m_value); }

    /**
     * Get the value as rational (throws an
     * exception if it doesn't hold a rational).
     */
    const RationalType &get_rational() const { return std::get<RationalType>(m_value); }

    /**
     * Get the value as rational.
     * A copy is created in any case;
     * returns a valid rational value even
     * if the variant holds a platform int.
     */
    RationalType as_rational() const {
        return visit(multi_visitor{[](const IntType &i) -> RationalType { return RationalType(i); },
                                   [](const RationalType &r) -> RationalType { return r; }});
    }

    CGAL::Sign sign() const {
        return visit(multi_visitor{
            [](IntType i) -> CGAL::Sign { return (i == 0) ? CGAL::ZERO : (i > 0 ? CGAL::POSITIVE : CGAL::NEGATIVE); },
            [](const RationalType &r) -> CGAL::Sign { return CGAL::sign(r); }});
    }

    /**
     * Get the value of *this as platform int,
     * without checking whether it holds a platform int
     * (i.e., with UB if it doesn't).
     */
    IntType unchecked_platform_int() const noexcept { return *std::get_if<IntType>(&m_value); }

    /**
     * Get the value of *this as rational,
     * without checking whether it holds a rational
     * (i.e., with UB if it doesn't).
     */
    const RationalType &unchecked_rational() const noexcept { return *std::get_if<RationalType>(&m_value); }

    /**
     * Visit the value with a callable that can handle
     * both IntType and RationalType values.
     */
    template<typename MultiCallable> auto visit(MultiCallable &&callable) const {
#ifndef NDEBUG
        return std::visit(std::forward<MultiCallable>(callable), m_value);
#else
        if(std::holds_alternative<IntType>(m_value)) {
            return std::forward<MultiCallable>(callable)(*std::get_if<IntType>(&m_value));
        } else {
            return std::forward<MultiCallable>(callable)(*std::get_if<RationalType>(&m_value));
        }
#endif
    }

    /**
     * Visit the value with a callable that can handle
     * both IntType and RationalType values.
     */
    template<typename MultiCallable> auto visit(MultiCallable &&callable) {
#ifndef NDEBUG
        return std::visit(std::forward<MultiCallable>(callable), m_value);
#else
        if(std::holds_alternative<IntType>(m_value)) {
            return std::forward<MultiCallable>(callable)(*std::get_if<IntType>(&m_value));
        } else {
            return std::forward<MultiCallable>(callable)(*std::get_if<RationalType>(&m_value));
        }
#endif
    }

    /**
     * Return an interval guaranteed to contain the value.
     */
    CGAL::Interval_nt_advanced to_interval() const {
        return visit(multi_visitor{[](const RationalType &v) { return mwt::to_interval(v); },
                                   [](IntType v) { return CGAL::Interval_nt_advanced(v); }});
    }

    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT to_sqrt_nt() const {
        return visit([](const auto &v) { return mwt::rational_to_sqrt_nt(v); });
    }

    bool is_integer() const noexcept {
        return visit(multi_visitor{[](IntType) -> bool { return true; },
                                   [](const RationalType &r) -> bool { return mwt::rational_is_integer(r); }});
    }

    /**
     * Cast the value to the given type;
     * assumes the CastTarget is constructible
     * from double (preserving the actual value)
     * and from RationalType.
     * The reason for not accepting int-like
     * constructors on the CastTarget is that
     * it may end up invoking a double constructor
     * or conversion operator with an integer that
     * is not representable as double, resulting
     * in silent loss of precision.
     */
    template<typename CastTarget> CastTarget cast_to() const {
        return visit(multi_visitor{[](IntType x) -> CastTarget {
                                       if(x <= INT64_C(9007199254740992) && x >= -INT64_C(9007199254740992)) {
                                           return CastTarget(double(x));
                                       } else {
                                           std::uint64_t uval(x);
                                           if(x < 0) {
                                               uval = ~uval + 1;
                                               return CastTarget(-double(std::uint32_t(uval >> 32))) * 4294967296.0 -
                                                      double(std::uint32_t(uval));
                                           } else {
                                               return CastTarget(double(std::uint32_t(uval >> 32))) * 4294967296.0 +
                                                      double(std::uint32_t(uval));
                                           }
                                       }
                                   },
                                   [](const RationalType &r) -> CastTarget { return CastTarget(r); }});
    }

    /**
     * Check if the encoded number is zero.
     */
    bool is_zero() const noexcept {
        return visit(multi_visitor{[](const IntType &i) { return i == 0; },
                                   [](const RationalType &r) { return rational_is_zero(r); }});
    }

    /**
     * Get the raw std::variant containing an int or rational.
     */
    const Raw &raw() const noexcept { return m_value; }

    /**
     * Negate the stored value.
     */
    void negate() {
        visit(multi_visitor{[&](IntType &i) {
                                if(i != std::numeric_limits<IntType>::min()) {
                                    i = -i;
                                } else {
                                    RationalType &r = m_value.template emplace<RationalType>(i);
                                    r = -r;
                                }
                            },
                            [&](RationalType &r) { r = -r; }});
    }

    /**
     * Addition operator.
     */
    RationalOrInt &operator+=(const RationalOrInt &other) {
        std::visit(multi_visitor{[this](IntType &lhs, IntType rhs) {
                                     if(overflow_check_detail::add_to_i64(lhs, rhs, lhs)) {
                                         RationalType value(lhs);
                                         value += rhs;
                                         m_value.template emplace<RationalType>(std::move(value));
                                     }
                                 },
                                 [this](IntType &lhs, const RationalType &rhs) {
                                     RationalType value(lhs);
                                     value += rhs;
                                     m_value.template emplace<RationalType>(std::move(value));
                                 },
                                 [this](RationalType &lhs, IntType rhs) { lhs += rhs; },
                                 [this](RationalType &lhs, const RationalType &rhs) {
                                     lhs += rhs;
                                     auto pi = rational_to_platform_int(lhs);
                                     if(pi) {
                                         m_value.template emplace<IntType>(*pi);
                                     }
                                 }},
                   m_value, other.m_value);
        return *this;
    }

    template<typename OtherIntType, std::enable_if_t<std::is_integral_v<std::decay_t<OtherIntType>>, int> = 0>
    RationalOrInt &operator+=(OtherIntType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                if(overflow_check_detail::add_to_i64(platform_int, other, platform_int)) {
                                    RationalType value(platform_int);
                                    value += other;
                                    m_value.template emplace<RationalType>(std::move(value));
                                }
                            },
                            [&](RationalType &rational) { rational += other; }});
        return *this;
    }

    template<typename OtherRatType, std::enable_if_t<!std::is_integral_v<std::decay_t<OtherRatType>>, int> = 0,
             std::enable_if_t<!std::is_same_v<std::decay_t<OtherRatType>, RationalOrInt>, int> = 0,
             std::enable_if_t<std::is_convertible_v<OtherRatType, RationalType>, int> = 0>
    RationalOrInt &operator+=(OtherRatType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                RationalType value(platform_int);
                                value += other;
                                m_value.template emplace<RationalType>(std::move(value));
                            },
                            [&](RationalType &rational) { rational += other; }});
        p_to_int_if_possible();
        return *this;
    }

    /**
     * Subtraction operator.
     */
    RationalOrInt &operator-=(const RationalOrInt &other) {
        std::visit(multi_visitor{[this](IntType &lhs, IntType rhs) {
                                     if(overflow_check_detail::sub_from_i64(lhs, rhs, lhs)) {
                                         RationalType value(lhs);
                                         value -= rhs;
                                         m_value.template emplace<RationalType>(std::move(value));
                                     }
                                 },
                                 [this](IntType &lhs, const RationalType &rhs) {
                                     RationalType value(lhs);
                                     value -= rhs;
                                     m_value.template emplace<RationalType>(std::move(value));
                                 },
                                 [this](RationalType &lhs, IntType rhs) { lhs -= rhs; },
                                 [this](RationalType &lhs, const RationalType &rhs) {
                                     lhs -= rhs;
                                     auto pi = rational_to_platform_int(lhs);
                                     if(pi) {
                                         m_value.template emplace<IntType>(*pi);
                                     }
                                 }},
                   m_value, other.m_value);
        return *this;
    }

    template<typename OtherIntType, std::enable_if_t<std::is_integral_v<std::decay_t<OtherIntType>>, int> = 0>
    RationalOrInt &operator-=(OtherIntType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                if(overflow_check_detail::sub_from_i64(platform_int, other, platform_int)) {
                                    RationalType value(platform_int);
                                    value -= other;
                                    m_value.template emplace<RationalType>(std::move(value));
                                }
                            },
                            [&](RationalType &rational) { rational -= other; }});
        return *this;
    }

    template<typename OtherRatType, std::enable_if_t<!std::is_integral_v<std::decay_t<OtherRatType>>, int> = 0,
             std::enable_if_t<!std::is_same_v<std::decay_t<OtherRatType>, RationalOrInt>, int> = 0,
             std::enable_if_t<std::is_convertible_v<OtherRatType, RationalType>, int> = 0>
    RationalOrInt &operator-=(OtherRatType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                RationalType value(platform_int);
                                value -= other;
                                m_value.template emplace<RationalType>(std::move(value));
                            },
                            [&](RationalType &rational) { rational -= other; }});
        p_to_int_if_possible();
        return *this;
    }

    /**
     * Multiplication operator.
     */
    RationalOrInt &operator*=(const RationalOrInt &other) {
        std::visit(multi_visitor{[this](IntType &lhs, IntType rhs) {
                                     if(overflow_check_detail::mul_to_i64(lhs, rhs, lhs)) {
                                         IntType x = lhs;
                                         RationalType value(lhs);
                                         value *= rhs;
                                         m_value.template emplace<RationalType>(std::move(value));
                                     }
                                 },
                                 [this](IntType &lhs, const RationalType &rhs) {
                                     RationalType value(lhs);
                                     value *= rhs;
                                     m_value.template emplace<RationalType>(std::move(value));
                                     p_to_int_if_possible();
                                 },
                                 [this](RationalType &lhs, IntType rhs) {
                                     lhs *= rhs;
                                     p_to_int_if_possible();
                                 },
                                 [this](RationalType &lhs, const RationalType &rhs) {
                                     lhs *= rhs;
                                     p_to_int_if_possible();
                                 }},
                   m_value, other.m_value);
        return *this;
    }

    template<typename OtherIntType, std::enable_if_t<std::is_integral_v<std::decay_t<OtherIntType>>, int> = 0>
    RationalOrInt &operator*=(OtherIntType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                if(overflow_check_detail::mul_to_i64(platform_int, other, platform_int)) {
                                    RationalType value(platform_int);
                                    value *= other;
                                    m_value.template emplace<RationalType>(std::move(value));
                                }
                            },
                            [&](RationalType &rational) {
                                rational *= other;
                                p_to_int_if_possible();
                            }});
        return *this;
    }

    template<typename OtherRatType, std::enable_if_t<!std::is_integral_v<std::decay_t<OtherRatType>>, int> = 0,
             std::enable_if_t<!std::is_same_v<std::decay_t<OtherRatType>, RationalOrInt>, int> = 0,
             std::enable_if_t<std::is_convertible_v<OtherRatType, RationalType>, int> = 0>
    RationalOrInt &operator*=(OtherRatType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                RationalType value(platform_int);
                                value *= other;
                                m_value.template emplace<RationalType>(std::move(value));
                            },
                            [&](RationalType &rational) { rational *= other; }});
        p_to_int_if_possible();
        return *this;
    }

    /**
     * Division operator.
     */
    RationalOrInt &operator/=(const RationalOrInt &other) {
        std::visit(multi_visitor{[this](IntType &lhs, IntType rhs) {
                                     if(!overflow_check_detail::div_i64_is_exact(lhs, rhs, lhs)) {
                                         IntType lhscopy(lhs);
                                         m_value.template emplace<RationalType>(lhscopy, rhs);
                                     }
                                 },
                                 [this](IntType lhs, const RationalType &rhs) {
                                     RationalType &value = m_value.template emplace<RationalType>(lhs);
                                     value /= rhs;
                                     p_to_int_if_possible();
                                 },
                                 [this](RationalType &lhs, IntType rhs) {
                                     lhs /= rhs;
                                     p_to_int_if_possible();
                                 },
                                 [this](RationalType &lhs, const RationalType &rhs) {
                                     lhs /= rhs;
                                     p_to_int_if_possible();
                                 }},
                   m_value, other.m_value);
        return *this;
    }

    template<typename OtherIntType, std::enable_if_t<std::is_integral_v<std::decay_t<OtherIntType>>, int> = 0>
    RationalOrInt &operator/=(OtherIntType &&other) {
        visit(multi_visitor{[&](IntType &platform_int) {
                                if(!overflow_check_detail::div_i64_is_exact(platform_int, other, platform_int)) {
                                    IntType v = platform_int;
                                    RationalType &value = m_value.template emplace<RationalType>(v);
                                    value /= other;
                                }
                            },
                            [&](RationalType &rational) {
                                rational /= other;
                                p_to_int_if_possible();
                            }});
        return *this;
    }

    template<typename OtherRatType, std::enable_if_t<!std::is_integral_v<std::decay_t<OtherRatType>>, int> = 0,
             std::enable_if_t<!std::is_same_v<std::decay_t<OtherRatType>, RationalOrInt>, int> = 0,
             std::enable_if_t<std::is_convertible_v<OtherRatType, RationalType>, int> = 0>
    RationalOrInt &operator/=(OtherRatType &&other) {
        visit(multi_visitor{[&](IntType platform_int) {
                                RationalType &value = m_value.template emplace<RationalType>(platform_int);
                                value /= other;
                            },
                            [&](RationalType &rational) { rational /= other; }});
        p_to_int_if_possible();
        return *this;
    }

    explicit operator bool() const noexcept { return !is_zero(); }

    bool operator!() const noexcept { return is_zero(); }

    template<typename OtherNumType,
             typename D = decltype(std::declval<const RationalType &>() == std::declval<const OtherNumType &>()),
             std::enable_if_t<std::is_convertible_v<D, bool>, int> = 0>
    bool operator==(const OtherNumType &other) const {
        using RawOther = std::decay_t<OtherNumType>;
        auto compare_rational = [&](const RationalType &r) -> bool { return r == other; };
        auto compare_int = [&](IntType i) -> bool {
            if constexpr(std::is_integral_v<RawOther>) {
                if constexpr(std::is_signed_v<RawOther>) {
                    return i == other;
                } else {
                    return i >= 0 && i == other;
                }
            } else {
                RationalType tmp(i);
                return tmp == other;
            }
        };
        return visit(multi_visitor{compare_int, compare_rational});
    }

    bool operator==(const RationalOrInt &other) const {
        return std::visit([](const auto &v1, const auto &v2) { return v1 == v2; }, m_value, other.m_value);
    }

    template<typename OtherNumType,
             typename = decltype(std::declval<const RationalOrInt &>() == std::declval<const OtherNumType &>())>
    bool operator!=(const OtherNumType &other) const {
        return !(*this == other);
    }

    template<typename OtherNumType,
             typename D = decltype(std::declval<const RationalType &>() < std::declval<const OtherNumType &>()),
             std::enable_if_t<std::is_convertible_v<D, bool>, int> = 0>
    bool operator<(const OtherNumType &other) const {
        using RawOther = std::decay_t<OtherNumType>;
        auto compare_rational = [&](const RationalType &r) -> bool { return r < other; };
        auto compare_int = [&](IntType i) -> bool {
            if constexpr(std::is_integral_v<RawOther>) {
                if constexpr(std::is_signed_v<RawOther>) {
                    return i < other;
                } else {
                    return i < 0 || i < other;
                }
            } else {
                RationalType tmp(i);
                return tmp < other;
            }
        };
        return visit(multi_visitor{compare_int, compare_rational});
    }

    bool operator<(const RationalOrInt &other) const {
        return std::visit([](const auto &v1, const auto &v2) { return v1 < v2; }, m_value, other.m_value);
    }

    template<typename OtherNumType,
             typename D = decltype(std::declval<const RationalType &>() < std::declval<const OtherNumType &>()),
             std::enable_if_t<std::is_convertible_v<D, bool>, int> = 0>
    bool operator>(const OtherNumType &other) const {
        using RawOther = std::decay_t<OtherNumType>;
        auto compare_rational = [&](const RationalType &r) -> bool { return r > other; };
        auto compare_int = [&](IntType i) -> bool {
            if constexpr(std::is_integral_v<RawOther>) {
                if constexpr(std::is_signed_v<RawOther>) {
                    return i > other;
                } else {
                    return i >= 0 && i > other;
                }
            } else {
                RationalType tmp(i);
                return tmp > other;
            }
        };
        return visit(multi_visitor{compare_int, compare_rational});
    }

    bool operator>(const RationalOrInt &other) const {
        return std::visit([](const auto &v1, const auto &v2) { return v1 > v2; }, m_value, other.m_value);
    }

    template<typename OtherNumType,
             typename = decltype(std::declval<const RationalOrInt &>() < std::declval<const OtherNumType &>())>
    bool operator>=(const OtherNumType &other) const {
        return !(*this < other);
    }

    template<typename OtherNumType,
             typename = decltype(std::declval<const RationalOrInt &>() < std::declval<const OtherNumType &>())>
    bool operator<=(const OtherNumType &other) const {
        return !(*this > other);
    }

    RationalOrInt operator+() const { return *this; }

    RationalOrInt operator-() const {
        return visit(multi_visitor{[](const IntType &i) {
                                       if(i == std::numeric_limits<IntType>::min()) {
                                           RationalType value(i);
                                           value = -value;
                                           return RationalOrInt(std::move(value));
                                       }
                                       return RationalOrInt(-i);
                                   },
                                   [](const RationalType &r) { return RationalOrInt(-r); }});
    }

  private:
    template<typename OtherIntType> void p_convert_or_rational(OtherIntType i) {
        std::int64_t result;
        if(overflow_check_detail::conversion_to_i64(i, result)) {
            m_value = IntType(result);
        } else {
            m_value = RationalType(i);
        }
    }

    void p_to_int_if_possible() {
        auto *r = std::get_if<RationalType>(&m_value);
        auto pi = rational_to_platform_int(*r);
        if(pi) {
            m_value.template emplace<IntType>(*pi);
        }
    }

    Raw m_value;
};

namespace rational_or_int_interoperable_detail {

template<typename OtherType, typename RationalType,
         typename D = decltype(std::declval<RationalOrInt<RationalType> &>() += std::declval<const OtherType &>()),
         std::enable_if_t<std::is_same_v<D, RationalOrInt<RationalType> &>, int> = 0>
std::true_type is_rational_or_int_interoperable_impl(OtherType *, RationalOrInt<RationalType> *);

template<typename OtherType> std::false_type is_rational_or_int_interoperable_impl(OtherType *, ...);

} // namespace rational_or_int_interoperable_detail

template<typename RationalType, typename OtherType>
constexpr bool is_rational_or_int_interoperable_v =
    decltype(rational_or_int_interoperable_detail::is_rational_or_int_interoperable_impl(
        (OtherType *)nullptr, (RationalOrInt<RationalType> *)nullptr))::value;

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator+(const RationalOrInt<RationalType> &lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(lhs);
    result += rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator+(RationalOrInt<RationalType> &&lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(std::move(lhs));
    result += rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator+(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    RationalOrInt<RationalType> result(rhs);
    result += lhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator+(const OtherType &lhs, RationalOrInt<RationalType> &&rhs) {
    RationalOrInt<RationalType> result(std::move(rhs));
    result += lhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator-(const RationalOrInt<RationalType> &lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(lhs);
    result -= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator-(RationalOrInt<RationalType> &&lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(std::move(lhs));
    result -= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator-(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    RationalOrInt<RationalType> result(rhs);
    result.negate();
    result += lhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator-(const OtherType &lhs, RationalOrInt<RationalType> &&rhs) {
    RationalOrInt<RationalType> result(std::move(rhs));
    result.negate();
    result += lhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator*(const RationalOrInt<RationalType> &lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(lhs);
    result *= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator*(RationalOrInt<RationalType> &&lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(std::move(lhs));
    result *= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator*(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    RationalOrInt<RationalType> result(rhs);
    result *= lhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator*(const OtherType &lhs, RationalOrInt<RationalType> &&rhs) {
    RationalOrInt<RationalType> result(std::move(rhs));
    result *= lhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator/(const RationalOrInt<RationalType> &lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(lhs);
    result /= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0>
inline RationalOrInt<RationalType> operator/(RationalOrInt<RationalType> &&lhs, const OtherType &rhs) {
    RationalOrInt<RationalType> result(std::move(lhs));
    result /= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline RationalOrInt<RationalType> operator/(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    RationalOrInt<RationalType> result(lhs);
    result /= rhs;
    return result;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline bool operator==(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    return rhs == lhs;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline bool operator!=(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    return rhs != lhs;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline bool operator<(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    return rhs > lhs;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline bool operator>(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    return rhs < lhs;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline bool operator<=(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    return rhs >= lhs;
}

template<typename RationalType, typename OtherType,
         std::enable_if_t<is_rational_or_int_interoperable_v<RationalType, OtherType>, int> = 0,
         std::enable_if_t<!std::is_same_v<RationalOrInt<RationalType>, std::decay_t<OtherType>>, int> = 0>
inline bool operator>=(const OtherType &lhs, const RationalOrInt<RationalType> &rhs) {
    return rhs <= lhs;
}

template<typename PrintableRational,
         typename = decltype(std::declval<std::ostream &>() << std::declval<const PrintableRational &>())>
inline std::ostream &operator<<(std::ostream &out, const RationalOrInt<PrintableRational> &value) {
    value.visit([&](const auto &a) { out << a; });
    return out;
}

} // namespace mwt

#endif
