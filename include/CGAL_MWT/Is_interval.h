#ifndef CGAL_MWT_IS_EXACT_H_INCLUDED_
#define CGAL_MWT_IS_EXACT_H_INCLUDED_

#include <CGAL/Interval_nt.h>
#include <type_traits>

namespace mwt {

template<typename FT> struct is_interval_number_type : std::false_type {};

template<bool P> struct is_interval_number_type<CGAL::Interval_nt<P>> : std::true_type {};

template<typename FT> constexpr bool is_interval_number_type_v = is_interval_number_type<FT>::value;

} // namespace mwt

#endif // CGAL_MWT_IS_EXACT_H_INCLUDED_
