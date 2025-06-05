#ifndef CGAL_MWT_TIME_UTIL_H_INCLUDED_
#define CGAL_MWT_TIME_UTIL_H_INCLUDED_

#include <chrono>

namespace mwt {

template<typename D> double to_secs(D duration) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}

template<typename Callable> double measure_time(Callable &&f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto end = std::chrono::high_resolution_clock::now();
    return to_secs(end - start);
}

} // namespace mwt

#endif
