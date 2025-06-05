#ifndef CGAL_MWT_LP_BACKENDS_H_INCLUDED_
#define CGAL_MWT_LP_BACKENDS_H_INCLUDED_

#include <mutex>

namespace mwt {

enum class LPStatus { OPTIMAL, INFEASIBLE, UNBOUNDED, TIMEOUT };

enum class BasicStatus { LOWER_BOUND, UPPER_BOUND, BASIC };

enum class ExtendedBasicStatus { LOWER_BOUND, UPPER_BOUND, BASIC, SUPERBASIC };

template<typename LPBackend> inline bool is_runtime_available() {
    static std::once_flag check_was_done;
    static bool is_available;
    std::call_once(check_was_done, [&]() { is_available = LPBackend::test_runtime_availability(true); });
    return is_available;
}

} // namespace mwt

#endif
