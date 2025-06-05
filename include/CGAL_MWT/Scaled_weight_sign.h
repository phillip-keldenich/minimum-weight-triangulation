#ifndef CGAL_MWT_SCALED_WEIGHT_SIGN_H_INCLUDED_
#define CGAL_MWT_SCALED_WEIGHT_SIGN_H_INCLUDED_

#include "Compare_weights.h"
#include "Rational_or_int.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/FPU.h>
#include <vector>

namespace mwt {

using ExactWithSqrtFT = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT;

template<typename Halfedge_, typename InnerNT, typename ExactFT_> class Scaled_weight_sign {
  public:
    using RationalOrInt = mwt::RationalOrInt<InnerNT>;
    using ExactFT = ExactFT_;
    using Halfedge = Halfedge_;
    using Interval = CGAL::Interval_nt_advanced;

    Scaled_weight_sign() = default;

    template<typename CoefficientType>
    void add_halfedge_with_coefficient(CoefficientType &&coeff, const Halfedge *halfedge) {
        if(halfedge->source() == halfedge->target())
            return;
        m_edge_entries.emplace_back(std::forward<CoefficientType>(coeff), halfedge);
    }

    void clear() noexcept {
        m_edge_entries.clear();
        m_squared_weights.clear();
    }

    CGAL::Sign compute_sign() {
        p_merge_edges();
        if(m_edge_entries.empty())
            return CGAL::ZERO;
        p_compute_squared_weights();
        p_merge_squared_weights();
        if(m_squared_weights.empty())
            return CGAL::ZERO;
        p_merge_coefficients_into_weights();
        if(m_negative_squared_weights.empty())
            return CGAL::POSITIVE;
        if(m_positive_squared_weights.empty())
            return CGAL::NEGATIVE;
        p_positive_negative_merge();
        auto interval_sign = p_compute_interval_sign();
        if(interval_sign)
            return *interval_sign;
        return p_compute_definitive_answer();
    }

  private:
    template<typename Iterator> static ExactWithSqrtFT p_range_to_log_height_sqrt_nt(Iterator begin, Iterator end) {
        ExactWithSqrtFT result(0);
        auto count = std::distance(begin, end);
        if(count < 8) {
            for(auto it = begin; it != end; ++it) {
                result += CGAL::sqrt(to_exact_with_sqrt(*it));
            }
        } else {
            auto mid = begin + count / 2;
            ExactWithSqrtFT result = p_range_to_log_height_sqrt_nt(begin, mid);
            result += p_range_to_log_height_sqrt_nt(mid, end);
        }
        return result;
    }

    CGAL::Sign p_compute_definitive_answer() const {
        ExactWithSqrtFT diff =
            p_range_to_log_height_sqrt_nt(m_positive_squared_weights.begin(), m_positive_squared_weights.end());
        diff -= p_range_to_log_height_sqrt_nt(m_negative_squared_weights.begin(), m_negative_squared_weights.end());
        return CGAL::sign(diff);
    }

    std::optional<CGAL::Sign> p_compute_interval_sign() {
        CGAL::Protect_FPU_rounding protect_round;
        Interval diff(0);
        std::size_t max_index = (std::min)(m_negative_squared_weights.size(), m_positive_squared_weights.size());
        for(std::size_t i = 0; i < max_index; ++i) {
            Interval spos = mwt::exact_to_interval(m_positive_squared_weights[i]);
            Interval sneg = mwt::exact_to_interval(m_negative_squared_weights[i]);
            spos = CGAL::sqrt(spos);
            sneg = CGAL::sqrt(sneg);
            spos -= sneg;
            diff += spos;
        }
        for(std::size_t i = max_index; i < m_negative_squared_weights.size(); ++i) {
            Interval sneg = mwt::exact_to_interval(m_negative_squared_weights[i]);
            sneg = CGAL::sqrt(sneg);
            diff -= sneg;
        }
        for(std::size_t i = max_index; i < m_positive_squared_weights.size(); ++i) {
            Interval spos = mwt::exact_to_interval(m_positive_squared_weights[i]);
            spos = CGAL::sqrt(spos);
            diff += spos;
        }
        if(diff.sup() < 0)
            return CGAL::NEGATIVE;
        if(diff.inf() > 0)
            return CGAL::POSITIVE;
        if(diff.is_point() && diff.inf() == 0)
            return CGAL::ZERO;
        return std::nullopt;
    }

    void p_compute_squared_weights() {
        m_squared_weights.reserve(m_edge_entries.size());
        for(const auto &eentry : m_edge_entries) {
            m_squared_weights.emplace_back(eentry.coefficient, eentry.halfedge);
        }
    }

    struct EdgeEntry {
        RationalOrInt coefficient;
        const Halfedge *halfedge;

        EdgeEntry(const RationalOrInt &coeff, const Halfedge *he) : coefficient(coeff), halfedge(he->primary_edge()) {}

        EdgeEntry(RationalOrInt &&coeff, const Halfedge *he) : coefficient(coeff), halfedge(he->primary_edge()) {}

        static bool compare_halfedge_id(const EdgeEntry &a, const EdgeEntry &b) { return a.halfedge < b.halfedge; }
    };

    struct SquaredWeightEntry {
        SquaredWeightEntry(const RationalOrInt &coeff, const Halfedge *he)
            : coefficient(coeff), squared_weight(p_squared_weight(*he)) {}

        RationalOrInt coefficient;
        ExactFT squared_weight;

        static bool compare_by_squared_weight(const SquaredWeightEntry &a, const SquaredWeightEntry &b) {
            return a.squared_weight < b.squared_weight;
        }

      private:
        static ExactFT p_squared_weight(const Halfedge &he) {
            const auto &p1 = he.source();
            const auto &p2 = he.target();
            ExactFT dx(p1.x());
            ExactFT dy(p1.y());
            dx -= p2.x();
            dy -= p2.y();
            dx *= dx;
            dy *= dy;
            dx += dy;
            return dx;
        }
    };

    std::vector<EdgeEntry> m_edge_entries;
    std::vector<SquaredWeightEntry> m_squared_weights;
    std::vector<ExactFT> m_negative_squared_weights;
    std::vector<ExactFT> m_positive_squared_weights;

    template<typename IteratorType, typename ExtractPrevious, typename ComparePrevious, typename Combine>
    IteratorType p_merge_elements(IteratorType begin, IteratorType end, ExtractPrevious &&extract,
                                  ComparePrevious &&compare_previous, Combine &&combine) {
        if(begin == end)
            return end;
        auto prev_value = extract(*begin);
        auto it = begin;
        bool found = false;
        for(++it; it != end; ++it) {
            auto this_value = extract(*it);
            if(compare_previous(prev_value, this_value)) {
                found = true;
                --it;
                break;
            }
            prev_value = std::move(this_value);
        }
        if(!found)
            return end;
        auto out_it = it;
        for(++it; it != end; ++it) {
            auto this_value = extract(*it);
            if(compare_previous(prev_value, this_value)) {
                combine(*out_it, *it);
            } else {
                ++out_it;
                *out_it = std::move(*it);
                prev_value = extract(*out_it);
            }
        }
        return ++out_it;
    }

    /**
     * Merge copies of the same edge into one coefficient.
     * Remove zero coefficients.
     */
    void p_merge_edges() {
        std::sort(m_edge_entries.begin(), m_edge_entries.end(), EdgeEntry::compare_halfedge_id);
        auto new_end = p_merge_elements(
            m_edge_entries.begin(), m_edge_entries.end(), [](const EdgeEntry &e) { return e.halfedge; },
            [](const Halfedge *prev, const Halfedge *cur) -> bool { return prev == cur; },
            [](EdgeEntry &out, const EdgeEntry &in) { out.coefficient += in.coefficient; });
        m_edge_entries.erase(new_end, m_edge_entries.end());
        m_edge_entries.erase(std::remove_if(m_edge_entries.begin(), m_edge_entries.end(),
                                            [](const EdgeEntry &e) { return !e.coefficient; }),
                             m_edge_entries.end());
    }

    /**
     * Merge edges with exactly the same squared weight into one coefficient.
     */
    void p_merge_squared_weights() {
        std::sort(m_squared_weights.begin(), m_squared_weights.end(), SquaredWeightEntry::compare_by_squared_weight);
        auto new_end = p_merge_elements(
            m_squared_weights.begin(), m_squared_weights.end(),
            [](const SquaredWeightEntry &e) { return &e.squared_weight; },
            [](const ExactFT *prev, const ExactFT *cur) { return *prev == *cur; },
            [](SquaredWeightEntry &out, const SquaredWeightEntry &in) { out.coefficient += in.coefficient; });
        m_squared_weights.erase(new_end, m_squared_weights.end());
        m_squared_weights.erase(std::remove_if(m_squared_weights.begin(), m_squared_weights.end(),
                                               [](const SquaredWeightEntry &e) { return !e.coefficient; }),
                                m_squared_weights.end());
    }

    /**
     * Merge the coefficients into the squared weights.
     */
    void p_merge_coefficients_into_weights() {
        m_negative_squared_weights.clear();
        m_positive_squared_weights.clear();
        for(const auto &e : m_squared_weights) {
            auto rational = e.coefficient.as_rational();
            ExactFT *out;
            if(e.coefficient < 0) {
                m_negative_squared_weights.push_back(e.squared_weight);
                rational = -rational;
                out = &m_negative_squared_weights.back();
            } else {
                m_positive_squared_weights.push_back(e.squared_weight);
                out = &m_positive_squared_weights.back();
            }
            rational *= rational;
            *out *= rational;
        }
    }

    /**
     * Eliminate common terms that occur both positive and negative.
     */
    void p_positive_negative_merge() {
        std::sort(m_negative_squared_weights.begin(), m_negative_squared_weights.end());
        std::sort(m_positive_squared_weights.begin(), m_positive_squared_weights.end());
        auto neg_it = m_negative_squared_weights.begin(), neg_end = m_negative_squared_weights.end();
        auto pos_it = m_positive_squared_weights.begin(), pos_end = m_positive_squared_weights.end();

        auto advance_pos = [&]() { return ++pos_it == pos_end; };

        auto advance_neg = [&]() { return ++neg_it == neg_end; };

        auto advance_both = [&]() {
            *pos_it++ = 0;
            *neg_it++ = 0;
            return pos_it == pos_end || neg_it == neg_end;
        };

        bool any_match = false;
        for(;;) {
            auto comp_result = CGAL::compare(*neg_it, *pos_it);
            if(comp_result == CGAL::SMALLER) {
                if(advance_neg())
                    break;
            } else if(comp_result == CGAL::LARGER) {
                if(advance_pos())
                    break;
            } else {
                any_match = true;
                if(advance_both())
                    break;
            }
        }
        if(any_match) {
            auto is_zero = [](const auto &x) { return x == 0; };
            auto new_pos_end =
                std::remove_if(m_positive_squared_weights.begin(), m_positive_squared_weights.end(), is_zero);
            auto new_neg_end =
                std::remove_if(m_negative_squared_weights.begin(), m_negative_squared_weights.end(), is_zero);
            m_positive_squared_weights.erase(new_pos_end, m_positive_squared_weights.end());
            m_negative_squared_weights.erase(new_neg_end, m_negative_squared_weights.end());
        }
    }
};

} // namespace mwt

#endif
