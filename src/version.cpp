#ifndef CGAL_MWT_VERSION
#error "VERSION NOT DEFINED ON COMMAND LINE!"
#endif

#include <CGAL_MWT/version.h>

namespace mwt {

std::string get_version() { return CGAL_MWT_VERSION; }

} // namespace mwt
