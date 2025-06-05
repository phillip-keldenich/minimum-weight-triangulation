#ifndef CGAL_MWT_TRAITS_H_INCLUDED_
#define CGAL_MWT_TRAITS_H_INCLUDED_

#include <boost/optional.hpp>
#include <type_traits>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/predicates/kernel_ftC2.h>

namespace mwt {

/**
 * Predicate to compute what we call 'pseudo angles'.
 * We need to store and compare intervals of polar angles.
 * The actual value is not important as long as the relative
 * ordering is correct. This predicate is much
 * faster than atan2 or dot-product angle calculations.
 * The predicate gets compiled to branch-free assembly for floating point types
 * and is probably the fastest way to do it.
 * It maps: [0, PI] --> [0, 2] and (PI, 2 PI] --> [-2, 0].
 * Discontinuity is at PI were it jumps from 2 to -2, similar to atan2.
 */
template<typename FT> struct Compute_pseudo_angle {
    static constexpr double min_polar_angle = -2.0;
    static constexpr double max_polar_angle = 2.0;

    template<typename T = FT>
    FT operator()(typename std::enable_if<std::is_floating_point<T>::value, T>::type x,
                  typename std::enable_if<std::is_floating_point<T>::value, T>::type y) const {
        return std::copysign(FT{1} - x / (std::fabs(x) + std::fabs(y)), y);
    }

    template<typename T = FT>
    FT operator()(const typename std::enable_if<!std::is_floating_point<T>::value, T>::type &x,
                  const typename std::enable_if<!std::is_floating_point<T>::value, T>::type &y) const {
        FT r = CGAL::abs(FT{1} - x / (CGAL::abs(x) + CGAL::abs(y)));
        if constexpr(std::is_same_v<std::decay_t<decltype(y >= 0)>, bool>) {
            // non-interval type/type which returns certain bool on comparisons
            return y >= 0 ? r : -r;
        } else {
            // interval type/type which returns uncertain on comparisons
            auto c = (y >= 0);
            if(BOOST_LIKELY(c.is_certain())) {
                return c ? r : -r;
            }
            auto mr = -r;
            return FT((std::min)(r.inf(), mr.inf()), (std::max)(r.sup(), mr.sup()));
        }
    }
};

namespace detail {

/**
 * The Diamond test:
 *
 * Allows to exclude edges from consideration as MWT edges
 * based on local properties.
 * Proved by Drysdale, McElfresh and Snoeyink in
 *  "An improved diamond property for minimum-weight triangulation" (1998)
 *
 *  An edge e can only be part of the MWT if at least one of the adjacent isosceles triangles
 *  with base angle alpha = PI / 4.6 is empty, i.e. it contains no other points of the point set in question.
 *
 *  This means that, for a conservative exclusion criterion, we have to make the triangles
 *  as small as they can possibly be (from a numerical point of view,
 *  since the angle PI / 4.6 is not representable).
 */

/**
 * The largest double value that is at most as large as the angle
 * used for the diamond test (exactly given as decimal value).
 * The value is not used below since all the constants are computed at much higher precision.
 * If this value is changed because of a tighter analysis of the diamond property,
 * all constants below have to be updated, but this is the only change that should be necessary.
 * The real value of PI / 4.6 is about           0.68295492469343331270927030071293540960808030421197952...
 * according to WolframAlpha.                                     >=
 */
static constexpr double DIAMOND_PROPERTY_ANGLE = 0.6829549246934332185077209942392073571681976318359375;

/**
 * We actually use half the tan value of the angle DIAMOND_PROPERTY_ANGLE;
 * again, we should underapproximate by taking the largest double value that
 * is at most this value.
 * The real value is about             0.406780171881322475658325127067534835325828539597578886...
 * according to WolframAlpha.                            >=
 */
static constexpr double HALF_TAN_DPA = 0.406780171881322427651639372925274074077606201171875;

/**
 * Squared cosine of 2 * DIAMOND_PROPERTY_ANGLE (upper bound).
 * The real value is about                 0.0413943492472734910780972760171922531548529138942998916030902714...
 * We need to make this slightly larger:                      <=
 */
static constexpr double SQUARED_COS_2DPA = 0.041394349247273497238719386359662166796624660491943359375;

/**
 * Squared value of 2cos(DIAMOND_PROPERTY_ANGLE) (upper bounded).
 * We use this to compute the distance at which a dead sector becomes active.
 * The real value is about                 2.4069120261052675797560574441231568535556673337801938479...
 *                                                          <=
 */
static constexpr double SQUARED_2COS_DPA = 2.406912026105267887743366372887976467609405517578125;

/**
 * Squared value of cos(DIAMOND_PROPERTY_ANGLE) (upper bounded).
 * We use this to check whether we are in the case where the dead sector
 * is bounded only by the two inducing vertices.
 * The real value is about                0.6017280065263168949390143610307892133889168334450484619811...
 *                                                         <=
 */
static constexpr double SQUARED_COS_DPA = 0.60172800652631697193584159322199411690235137939453125;

/**
 * Sine of DIAMOND_PROPERTY_ANGLE (upper bounded).
 * We use this to compute the dead sector angles in the case
 * where it is not two adjacent points that delimit a
 * dead sector.
 * In that case, erring on the side of caution means
 * using a slightly smaller angle (i.e., slightly smaller sine).
 * The real value is about           0.63108794432605278936740013014331057420082924774022172267785922...
 *                                                     >=
 */
static constexpr double SIN_DPA = 0.63108794432605275215308893166366033256053924560546875;
static constexpr double SIN_DPA_UB = 0.6310879443260528631753913941793143749237060546875;

/**
 * Cosine of DIAMOND_PROPERTY_ANGLE (lower bounded).
 * We use this to compute the dead sector angles in the case
 * where it is not two adjacent points that delimit a dead sector.
 * In that case, erring on the side of caution means
 * using a slightly smaller angle (i.e., slightly larger cosine).
 * The real value is about           0.77571129070441980704110101096953689558772273032391012075560931...
 *                                                    <=
 */
static constexpr double COS_DPA = 0.7757112907044199090478286962024867534637451171875;
static constexpr double COS_DPA_LB = 0.77571129070441979802552623368683271110057830810546875;

/**
 * Work around a CGAL bug (missing <FT> in call to circumcenter_translateC2 in kernel_ftC2.h).
 */
template<typename FT>
static inline FT workaround_squared_radiusC2(const FT &px, const FT &py, const FT &qx, const FT &qy, const FT &rx,
                                             const FT &ry, FT &x, FT &y) {
    CGAL::circumcenter_translateC2<FT>(qx - px, qy - py, rx - px, ry - py, x, y);
    FT r2 = CGAL_NTS square(x) + CGAL_NTS square(y);
    x += px;
    y += py;
    return r2;
}

} // namespace detail

template<typename Kernel_> class Mwt_traits_2 : public Kernel_ {
  public:
    /**
     * The kernel we are extending by additional predicates.
     * Usually, we will use CGAL's EPICK kernel.
     */
    using Kernel = Kernel_;

    /**
     * A 2d point.
     */
    using Point_2 = typename Kernel::Point_2;

    /**
     * A 2d direction.
     */
    using Direction_2 = typename Kernel::Direction_2;

    /**
     * The imprecise filtering kernel we use for our
     * filtered predicates.
     */
    using FK = CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>;

    /**
     * The exact kernel we use for our filtered predicates,
     * should the imprecise kernel not suffice.
     */
    using EK = CGAL::Simple_cartesian<CGAL::Exact_rational>;

    /**
     * Conversion between the original cartesian kernel
     * (e.g., EPICK) and EK/FK.
     */
    using C2E = CGAL::Cartesian_converter<Kernel, EK>;
    using C2F = CGAL::Cartesian_converter<Kernel, FK>;
    using E2C = CGAL::Cartesian_converter<EK, Kernel>;
    using F2C = CGAL::Cartesian_converter<FK, Kernel>;

    /**
     * Left turn predicate.
     */
    using Left_turn_2 = typename Kernel::Left_turn_2;

    /**
     * Check if a given quadrilateral is convex.
     * a, c, d and c, a, b must be left turns.
     */
    struct Quadrilateral_convex {
        bool operator()(const Point_2 &a, const Point_2 &b, const Point_2 &c, const Point_2 &d) const noexcept {
            Left_turn_2 lt{};
            CGAL_precondition(lt(a, c, d) && lt(c, a, b));
            return lt(b, d, a) && lt(d, b, c);
        }
    };

  private:
    /**
     * Predicate to directly compare two angles.
     */
    struct Angle_less {
        Angle_less(const Point_2 &p) : origin{p} {}

        bool operator()(const Point_2 &p, const Point_2 &q) const {
            return Direction_2(p - origin) < Direction_2(q - origin);
        }

      private:
        const Point_2 origin;
    };

    /**
     * Actually compute a polar angle (using some given predicate template).
     * We use this with the pseudo angle, but could use actual angles instead.
     */
    template<typename K, template<typename T> class A> struct Compute_polar_angle {
        using FT = typename K::FT;
        using Point_2 = typename K::Point_2;
        using Angle = A<FT>;

        Compute_polar_angle() : origin{} {}

        Compute_polar_angle(const Point_2 &p) : origin{p} {}

        FT operator()(const Point_2 &p) const { return Angle{}(p.x() - origin.x(), p.y() - origin.y()); }

        FT operator()(const FT &x, const FT &y) const { return Angle{}(x - origin.x(), y - origin.y()); }

      private:
        const Point_2 origin;
    };

    /**
     * Less-than comparison for polar angles.
     */
    template<typename K> struct Less_polar_angle {
        using FT = typename K::FT;
        using Point_2 = typename K::Point_2;
        using result_type = bool;

        Less_polar_angle(const Point_2 &p) : origin{p} {}

        result_type operator()(const Point_2 &p, const Point_2 &q) const {
            return CGAL::compare_angle_with_x_axisC2(p.x() - origin.x(), p.y() - origin.y(), q.x() - origin.x(),
                                                     q.y() - origin.y()) == CGAL::Comparison_result::SMALLER;
        }

      private:
        const Point_2 origin;
    };

    /**
     * Predicate for constructing (simplified) dead sectors.
     * A (simplified) dead sector is a region of points
     * which must all fail the diamond test w.r.t. a
     * source point s (based on an existing pair of points l,r
     * that populate the triangles on either side of the diamond).
     * A simplified dead sector consists of an interval of polar angles
     * and a (squared) activation distance. All points t within
     * the polar angle interval w.r.t. s that are at least the activation
     * distance away from s must fail the diamond test.
     */
    template<typename K, template<typename T> class A> struct Construct_dead_sector {
        using FT = typename K::FT;
        using Point_2 = typename K::Point_2;
        using Sector_2 = std::pair<FT, FT>;
        using result_type = boost::optional<Sector_2>;
        using Angle = A<FT>;

        /**
         * Used to compute the activation distance (it is squared_distance_coefficient() * sqdist),
         * where sqdist is the distance to the farther point from s between l and r.
         */
        constexpr static double squared_distance_coefficient() { return detail::SQUARED_2COS_DPA; }

        /**
         * Compute the dead sector (as boost::optional, allowing us to return 'no sector').
         */
        result_type operator()(const Point_2 &s, const Point_2 &l, const Point_2 &r) const {
            result_type sector{}; // initially, 'no sector' (this is a boost::optional)
            const FT &px = s.x(), &py = s.y();
            const FT &rx = r.x(), &ry = r.y();
            const FT &lx = l.x(), &ly = l.y();
            const FT vlx = lx - px, vly = ly - py;
            const FT vrx = rx - px, vry = ry - py;

            // angle computations between vl and vr (dot product)
            const FT dot = vlx * vrx + vly * vry;
            if(dot <= FT{0}) {
                // the angle is definitely too great (>= PI/2)
                return sector;
            }
            CGAL_assertion(CGAL::left_turn(s, r, l));

            // cos^2( <(l,r) ) * (||l|| * ||r||)^2 = (l * r)^2 = dot^2
            const FT dot_squared = dot * dot;
            const FT nl = vlx * vlx + vly * vly; // ||l||
            const FT nr = vrx * vrx + vry * vry; // ||r||
            const FT denom = (nl * nr);          // ||l|| * ||r||

            FT lower{0}, upper{0};
            if(dot_squared <= detail::SQUARED_COS_2DPA * denom) {
                return sector; // angle >= 2 * alpha = PI/2.3: No dead sector
            } else if(dot_squared >= detail::SQUARED_COS_DPA * denom) {
                // angle <= alpha: The dead sector is bounded by the two vertices
                Angle angle{};
                lower = angle(vrx, vry);
                upper = angle(vlx, vly);
            } else { // alpha < angle < 2 * alpha
                // Rotate right vector to the left (ccw)
                const FT slx = detail::COS_DPA * vrx - detail::SIN_DPA * vry;
                const FT sly = detail::SIN_DPA * vrx + detail::COS_DPA * vry;

                // Rotate left vector to the right (cw)
                const FT srx = detail::COS_DPA * vlx + detail::SIN_DPA * vly;
                const FT sry = -detail::SIN_DPA * vlx + detail::COS_DPA * vly;

                // The dead sector is the overlapping between the rotated vectors.
                Angle angle{};
                lower = angle(srx, sry);
                upper = angle(slx, sly);
            }

            // Depending on the kernel, numerical errors might
            // generate degenerate sectors (lower >= upper).
            // Returning a smaller (or an empty) sector is always fine.
            if(lower < upper) {
                sector.emplace(lower, upper);
            }
            return sector;
        }
    };

    /**
     * Implementation of the actual diamond test as predicate.
     */
    template<typename K> struct Diamond_test {
        using FT = typename K::FT;
        using Point_2 = typename K::Point_2;
        using Segment_2 = typename K::Segment_2;
        using result_type = bool;

        /**
         * Check the diamond property for the edge. Returns true iff the edge e can be excluded,
         * i.e., if both triangles are occupied by the points l and r.
         * Assumes that l lies on the left side of e and r on the right side.
         */
        result_type operator()(const Segment_2 &e, const Point_2 &l, const Point_2 &r) const {
            return both_triangles_occupied(e.source(), e.target(), l, r);
        }

        /**
         * Check the diamond property for the edge Segment_2{source, target},
         * like operator()(Segment, Point, Point).
         */
        result_type operator()(const Point_2 &source, const Point_2 &target, const Point_2 &left,
                               const Point_2 &right) const {
            return both_triangles_occupied(source, target, left, right);
        }

        /**
         * Check the diamond property for one of the two isoceles triangles,
         * i.e., check if l occupies the left isoceles triangle of e.
         * Assumes that l lies on the left side of e.
         */
        result_type operator()(const Segment_2 &e, const Point_2 &l) const {
            return left_triangle_occupied(e.source(), e.target(), l);
        }

        /**
         * Check the diamond property for one of the two isoceles triangles,
         * like operator()(Segment_2{source,target}, left).
         */
        result_type operator()(const Point_2 &source, const Point_2 &target, const Point_2 &left) const {
            return left_triangle_occupied(source, target, left);
        }

      private:
        result_type left_triangle_occupied(const Point_2 &s, const Point_2 &t, const Point_2 &l) const {
            const FT &sx = s.x();
            const FT &sy = s.y();
            const FT &tx = t.x();
            const FT &ty = t.y();

            // Midpoint of edge st
            const FT mx = FT{0.5} * (sx + tx);
            const FT my = FT{0.5} * (sy + ty);

            // Orthogonal vector to st with equal length
            const FT ox = (ty - sy);
            const FT oy = (sx - tx);

            // Apex of left triangle
            const FT alx = mx - detail::HALF_TAN_DPA * ox;
            const FT aly = my - detail::HALF_TAN_DPA * oy;

            // Points must be in CCW order
            return in_triangle(sx, sy, tx, ty, alx, aly, l.x(), l.y());
        }

        result_type both_triangles_occupied(const Point_2 &s, const Point_2 &t, const Point_2 &l,
                                            const Point_2 &r) const {
            const FT &sx = s.x();
            const FT &sy = s.y();
            const FT &tx = t.x();
            const FT &ty = t.y();

            // Midpoint of edge st
            const FT mx = FT{0.5} * (sx + tx);
            const FT my = FT{0.5} * (sy + ty);

            // Orthogonal vector to st with equal length
            const FT ox = (ty - sy);
            const FT oy = (sx - tx);

            // Apex of right triangle
            const FT arx = mx + detail::HALF_TAN_DPA * ox;
            const FT ary = my + detail::HALF_TAN_DPA * oy;

            // Apex of left triangle
            const FT alx = mx - detail::HALF_TAN_DPA * ox;
            const FT aly = my - detail::HALF_TAN_DPA * oy;

            // Points must be in CCW order
            return in_triangle(sx, sy, tx, ty, alx, aly, l.x(), l.y()) &&
                   in_triangle(sx, sy, arx, ary, tx, ty, r.x(), r.y());
        }

        // Assumption: Triangle abc is given in CCW order
        result_type in_triangle(const FT &ax, const FT &ay, const FT &bx, const FT &by, const FT &cx, const FT &cy,
                                const FT &px, const FT &py) const {
            return CGAL::orientationC2(px, py, ax, ay, bx, by) == CGAL::Orientation::COUNTERCLOCKWISE &&
                   CGAL::orientationC2(px, py, bx, by, cx, cy) == CGAL::Orientation::COUNTERCLOCKWISE &&
                   CGAL::orientationC2(px, py, cx, cy, ax, ay) == CGAL::Orientation::COUNTERCLOCKWISE;
        }
    };

    /**
     * Predicate to compare distance between two pairs of points.
     */
    template<typename K> struct Less_distance {
        using FT = typename K::FT;
        using Point_2 = typename K::Point_2;
        using result_type = bool;

        result_type operator()(const Point_2 &p, const Point_2 &q, const Point_2 &s, const Point_2 &t) const {
            K kernel{};
            auto d1 = kernel.compute_squared_distance_2_object()(p, q);
            auto d2 = kernel.compute_squared_distance_2_object()(s, t);
            return d1 < d2;
        }
    };

  public:
    using Diamond_test_2 = CGAL::Filtered_predicate<Diamond_test<EK>, Diamond_test<FK>, C2E, C2F>;

    template<typename K> using Compute_angle = Compute_pseudo_angle<K>;

    using Compute_polar_angle_2 = Compute_polar_angle<Kernel, Compute_angle>;
    using Construct_dead_sector_2 = Construct_dead_sector<Kernel, Compute_angle>;

    using Angle_less_2 = Angle_less;
    using Less_polar_angle_2 = CGAL::Filtered_predicate<Less_polar_angle<EK>, Less_polar_angle<FK>, C2E, C2F>;
    using Less_distance_2 = CGAL::Filtered_predicate<Less_distance<EK>, Less_distance<FK>, C2E, C2F>;

    Angle_less_2 angle_less_2_object(const Point_2 &p) const { return Angle_less_2{p}; }

    Less_polar_angle_2 less_polar_angle_2_object(const Point_2 &p) const { return Less_polar_angle_2{p}; }

    Less_distance_2 less_distance_2_object() const { return Less_distance_2{}; }

    Diamond_test_2 diamond_test_2_object() const { return Diamond_test_2{}; }

    Construct_dead_sector_2 construct_dead_sector_2_object() const { return Construct_dead_sector_2{}; };

    Compute_polar_angle_2 compute_polar_angle_2_object(const Point_2 &p) const { return Compute_polar_angle_2{p}; }

    Quadrilateral_convex quadrilateral_convex_object() const { return Quadrilateral_convex{}; }

    static constexpr double min_polar_angle() { return Compute_angle<Kernel>::min_polar_angle; }
    static constexpr double max_polar_angle() { return Compute_angle<Kernel>::max_polar_angle; }
};

} // namespace mwt

#endif // MWT_TRAITS_H
