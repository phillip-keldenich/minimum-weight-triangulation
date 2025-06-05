#ifndef CGAL_MWT_COMPARE_WEIGHTS_H_INCLUDED_
#define CGAL_MWT_COMPARE_WEIGHTS_H_INCLUDED_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL/Lazy_exact_nt.h>
#include <boost/multiprecision/number.hpp>
#include <type_traits>
#include <utility>

namespace mwt {

using ExactWithSqrtFT = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT;

namespace detail {

template<typename T> using Normalize = std::remove_cv_t<std::remove_reference_t<T>>;
template<typename T>
constexpr bool IsIntOrFloat = std::is_floating_point_v<Normalize<T>> || std::is_integral_v<Normalize<T>>;

/**
 * Test if CGAL::exact can be applied to T.
 */
template<typename T> struct HasCGALExactT {
    template<typename U> static std::true_type test(Normalize<decltype(CGAL::exact(std::declval<const U &>()))> *);
    template<typename U> static std::false_type test(...);
    using type = decltype(test<T>(nullptr));
    static constexpr bool value = type::value;
};
template<typename T> constexpr bool HasCGALExact = HasCGALExactT<T>::value;

/**
 * Test if T supports calling .mpq() on it.
 */
template<typename T> struct HasMPQ {
    template<typename U> static std::true_type test(Normalize<decltype(std::declval<const U &>().mpq())> *);
    template<typename U> static std::false_type test(...);
    using type = decltype(test<T>(nullptr));
    static constexpr bool value = type::value;
};

template<typename T> struct IsBoostNumber : std::false_type {};
template<typename Backend, boost::multiprecision::expression_template_option ET>
struct IsBoostNumber<boost::multiprecision::number<Backend, ET>> : std::true_type {};

template<typename ExactNumber> ExactWithSqrtFT exact_number_to_exact_with_sqrt(const ExactNumber &n) {
    if constexpr(IsBoostNumber<ExactNumber>::value) {
        // handle CGAL using boost multiprecision numbers
        const auto &backend = n.backend();
        return CORE::BigRat(backend.data());
    } else if constexpr(HasMPQ<ExactNumber>::value) {
        // should cover CGAL directly using GMP/MPIR as backend
        return CORE::BigRat(n.mpq());
    } else {
        static_assert(IsBoostNumber<ExactNumber>::value || HasMPQ<ExactNumber>::value,
                      "Unsupported exact number type for to_exact_with_sqrt!");
        throw std::logic_error("Unsupported exact number type for to_exact_with_sqrt!");
    }
}

} // namespace detail

/**
 * Conversion routine (not apparently possible in CGAL) for
 * CGAL number types to CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT.
 */
template<typename Number> ExactWithSqrtFT to_exact_with_sqrt(const Number &n) {
    using namespace mwt::detail;
    using N = Normalize<Number>;
    if constexpr(IsIntOrFloat<N>) {
        return ExactWithSqrtFT(n);
    } else if constexpr(std::is_same_v<N, ExactWithSqrtFT>) {
        return n;
    } else if constexpr(HasCGALExact<N>) {
        auto exact = CGAL::exact(n);
        return exact_number_to_exact_with_sqrt(exact);
    } else {
        return exact_number_to_exact_with_sqrt(n);
    }
}

/**
 * Class to compare sums of edge weights exactly.
 * Essentially, a few additional filters on top of
 * CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt,
 * mostly aiming at reducing the time taken on exactly
 * equal sums of edge weights that can happen due to
 * symmetry; we aim to filter out any squared edge weights that
 * appear in both sums.
 */
template<typename InputKernel> class CompareWeights {
  public:
    using InputPoint = CGAL::Point_2<InputKernel>;
    using InputSegment = CGAL::Segment_2<InputKernel>;
    using InputCoord = detail::Normalize<decltype(std::declval<InputPoint>().x())>;
    using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;
    using ExactFT = ExactKernel::FT;
    using SqrtKernel = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
    using SqrtFT = SqrtKernel::FT;

    /**
     * Clear both LHS and RHS.
     */
    void reset() {
        m_lhs.clear();
        m_rhs.clear();
    }

    /**
     * Add a distance to the left-hand side.
     */
    void add_lhs(const InputPoint &p1, const InputPoint &p2) { p_add<true>(p1, p2); }

    /**
     * Add a distance to the left-hand side.
     */
    void add_lhs(const InputSegment &s) { add_lhs(s.source(), s.target()); }

    /**
     * Add a distance to the right-hand side.
     */
    void add_rhs(const InputPoint &p1, const InputPoint &p2) { p_add<false>(p1, p2); }

    /**
     * Add a distance to the right-hand side.
     */
    void add_rhs(const InputSegment &s) { add_rhs(s.source(), s.target()); }

    /**
     * Get the sign (POSITIVE, NEGATIVE, ZERO) of LHS - RHS.
     */
    CGAL::Sign sign() {
        p_filter_common();
        if(p_have_definitive_sign()) {
            return p_get_definitive_sign();
        }
        return p_definitive_answer<Sign>();
    }

    /**
     * After computing the sign, get an interval that contains the result.
     */
    CGAL::Interval_nt_advanced compute_interval() {
        CGAL::Protect_FPU_rounding protect;
        if(m_rhs.empty()) {
            return p_sum_up_weights(m_lhs);
        }
        if(m_lhs.empty()) {
            return -p_sum_up_weights(m_rhs);
        }
        return p_definitive_answer<ComputeInterval>();
    }

    /**
     * Check if LHS - RHS is positive.
     */
    bool is_positive() {
        p_filter_common();
        if(p_have_definitive_sign()) {
            return p_get_definitive_sign() == CGAL::POSITIVE;
        }
        return p_definitive_answer<Positive>();
    }

    /**
     * Check if LHS - RHS is negative.
     */
    bool is_negative() {
        p_filter_common();
        if(p_have_definitive_sign()) {
            return p_get_definitive_sign() == CGAL::NEGATIVE;
        }
        return p_definitive_answer<Negative>();
    }

  private:
    CGAL::Interval_nt_advanced p_sum_up_weights(const std::vector<ExactFT> &squared_weights) {
        SqrtFT comp(0);
        std::for_each(squared_weights.begin(), squared_weights.end(),
                      [&](const auto &w) { comp += CGAL::sqrt(to_exact_with_sqrt(w)); });
        double d1 = -std::numeric_limits<double>::infinity(), d2 = std::numeric_limits<double>::infinity();
        comp.doubleInterval(d1, d2);
        return CGAL::Interval_nt_advanced(d1, d2);
    }

    bool p_have_definitive_sign() { return m_lhs.empty() || m_rhs.empty(); }

    CGAL::Sign p_get_definitive_sign() {
        if(m_lhs.empty()) {
            return m_rhs.empty() ? CGAL::ZERO : CGAL::NEGATIVE;
        }
        return CGAL::POSITIVE;
    }

    void p_filter_common() {
        std::sort(m_lhs.begin(), m_lhs.end());
        std::sort(m_rhs.begin(), m_rhs.end());
        auto lhs_cur = m_lhs.begin();
        auto rhs_cur = m_rhs.begin();
        auto lhs_end = m_lhs.end();
        auto rhs_end = m_rhs.end();

        // find the first equal element, if there is one.
        // needed to avoid self-move-assignment, which may or
        // may not be safe depending on the implementation
        // of the move-assignment operator of the number type
        // (and it is also faster than moving everything before a match)
        bool found_equal = false;
        while(lhs_cur != lhs_end && rhs_cur != rhs_end) {
            auto cmp_result = CGAL::compare(*lhs_cur, *rhs_cur);
            if(cmp_result == CGAL::SMALLER) {
                ++lhs_cur;
            } else if(cmp_result == CGAL::LARGER) {
                ++rhs_cur;
            } else {
                found_equal = true;
                break;
            }
        }
        if(!found_equal)
            return;

        auto lhs_out = lhs_cur++;
        auto rhs_out = rhs_cur++;
        while(lhs_cur != lhs_end && rhs_cur != rhs_end) {
            auto cmp_result = CGAL::compare(*lhs_cur, *rhs_cur);
            if(cmp_result == CGAL::SMALLER) {
                *lhs_out++ = std::move(*lhs_cur++);
            } else if(cmp_result == CGAL::LARGER) {
                *rhs_out++ = std::move(*rhs_cur++);
            } else {
                ++lhs_cur;
                ++rhs_cur;
            }
        }
        lhs_out = std::move(lhs_cur, lhs_end, lhs_out);
        rhs_out = std::move(rhs_cur, rhs_end, rhs_out);
        m_lhs.erase(lhs_out, lhs_end);
        m_rhs.erase(rhs_out, rhs_end);
    }

    template<typename ResultOp> typename ResultOp::result_type p_definitive_answer() {
        SqrtFT comp(0);
        std::for_each(m_lhs.begin(), m_lhs.end(),
                      [&](const auto &lhs) { comp += CGAL::sqrt(to_exact_with_sqrt(lhs)); });
        std::for_each(m_rhs.begin(), m_rhs.end(),
                      [&](const auto &rhs) { comp -= CGAL::sqrt(to_exact_with_sqrt(rhs)); });
        return ResultOp()(comp);
    }

    struct Sign {
        using result_type = CGAL::Sign;

        result_type operator()(const SqrtFT &comp) const { return CGAL::sign(comp); }
    };

    struct Positive {
        using result_type = bool;

        result_type operator()(const SqrtFT &comp) const { return CGAL::is_positive(comp); }
    };

    struct Negative {
        using result_type = bool;

        result_type operator()(const SqrtFT &comp) const { return CGAL::is_negative(comp); }
    };

    struct ComputeInterval {
        using result_type = CGAL::Interval_nt_advanced;

        result_type operator()(const SqrtFT &comp) const {
            double d1 = -std::numeric_limits<double>::infinity(), d2 = std::numeric_limits<double>::infinity();
            comp.doubleInterval(d1, d2);
            return CGAL::Interval_nt_advanced(d1, d2);
        }
    };

    template<bool LHS> void p_add(const InputPoint &p1, const InputPoint &p2) {
        if(p1 == p2)
            return;
        ExactFT x(p1.x());
        ExactFT y(p1.y());
        x -= ExactFT(p2.x());
        y -= ExactFT(p2.y());
        x *= x;
        y *= y;
        x += y;
        if constexpr(LHS) {
            m_lhs.push_back(x);
        } else {
            m_rhs.push_back(x);
        }
    }

    std::vector<ExactFT> m_lhs; // the left-hand side squared distances
    std::vector<ExactFT> m_rhs; // the right-hand side squared distances
};

} // namespace mwt

#endif
