#ifndef CGAL_MWT_GENERATE_RANDOM_INTEGRAL_H_INCLUDED_
#define CGAL_MWT_GENERATE_RANDOM_INTEGRAL_H_INCLUDED_

#include <CGAL/number_utils.h>
#include <algorithm>
#include <cstdint>
#include <random>
#include <type_traits>
#include <vector>

namespace mwt {

template<typename Kernel, typename RNG, std::enable_if_t<!std::is_arithmetic_v<RNG>, int> = 0>
std::vector<typename Kernel::Point_2> generate_random_integral(const typename Kernel::Iso_rectangle_2 &region,
                                                               std::size_t num_points, RNG &gen) {
    auto to_i = [](const auto &x) { return std::int64_t(CGAL::to_double(x)); };
    std::uniform_int_distribution<std::int64_t> dist_x(to_i(region.xmin()), to_i(region.xmax()));
    std::uniform_int_distribution<std::int64_t> dist_y(to_i(region.ymin()), to_i(region.ymax()));
    std::vector<typename Kernel::Point_2> points;
    points.reserve(num_points);
    for(std::size_t i = 0; i < num_points; ++i) {
        points.emplace_back(dist_x(gen), dist_y(gen));
    }
    std::sort(points.begin(), points.end(), typename Kernel::Less_xy_2{});
    points.erase(std::unique(points.begin(), points.end()), points.end());
    return points;
}

template<typename Kernel>
std::vector<typename Kernel::Point_2> generate_random_integral(const typename Kernel::Iso_rectangle_2 &region,
                                                               std::size_t num_points, std::size_t seed) {
    std::mt19937_64 gen(seed);
    return generate_random_integral<Kernel>(region, num_points, gen);
}

} // namespace mwt

#endif
