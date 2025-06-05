#ifndef CGAL_MWT_OUTPUT_H_INCLUDED_
#define CGAL_MWT_OUTPUT_H_INCLUDED_

#include "LMT_halfedge.h"
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <variant>

namespace mwt {

using InfoEntry = std::variant<std::string, std::int64_t, double, nlohmann::json>;
using InfoMap = std::map<std::string, InfoEntry>;

struct InfoEntryPrinter {
    explicit InfoEntryPrinter(std::ostream &output) : output(&output) {}

    void operator()(const std::string &str) const { *output << "\"" << str << "\""; }

    void operator()(double d) const {
        if(d == 0) {
            *output << "0.0";
        } else {
            *output << std::setprecision(19) << std::showpoint << d;
        }
    }

    void operator()(std::int64_t x) const { *output << x; }

    void operator()(const nlohmann::json &data) const { *output << data; }

    std::ostream *output;
};

inline std::ostream &operator<<(std::ostream &output, const InfoEntry &entry) {
    InfoEntryPrinter visitor(output);
    std::visit(visitor, entry);
    return output;
}

template<typename PointIterator> void output_only_points(std::ostream &output, PointIterator begin, PointIterator end) {
    output << std::setprecision(19);
    output << "{\n";
    output << "  \"points\": [\n";
    for(auto it = begin; it != end; ++it) {
        if(it != begin) {
            output << ",\n";
        }
        output << "    [" << it->x() << ", " << it->y() << "]";
    }
    output << "\n  ]\n";
    output << "}\n";
}

/**
 * Output the triangulation and additional information to the output stream
 * in JSON format; cannot use nlohmann::json directly because it does not
 * allow streaming the data but would require us to construct a JSON object
 * representing the triangulation in memory (which we don't have the spare memory for).
 */
template<typename Tree, typename Skeleton>
void output_triangulation(std::ostream &output, const Tree &tree, const Skeleton &skeleton,
                          bool include_possibility = false, const InfoMap &info_map = {}) {
    const auto *halfedges_begin = skeleton.halfedges_begin();
    const auto *halfedges_end = skeleton.halfedges_end();
    const auto points_begin = tree.points_begin();
    const auto points_end = tree.points_end();

    output << std::setprecision(19);
    output << "{\n";
    if(!info_map.empty()) {
        output << "  \"solution_info\": {";
        bool first = true;
        for(const auto &[k, v] : info_map) {
            if(!first)
                output << ", ";
            first = false;
            output << "\"" << k << "\": " << v;
        }
        output << "  },\n";
    }
    output << "  \"points\": [\n";
    for(auto it = points_begin; it != points_end; ++it) {
        if(it != points_begin) {
            output << ",\n";
        }
        output << "    [" << it->x() << ", " << it->y() << "]";
    }
    output << "\n  ],\n";
    output << "  \"edges\": [\n    ";
    bool first = true;
    for(auto it = halfedges_begin; it != halfedges_end; ++it) {
        if(it->status != mwt::LMTStatus::Certain && it->status != mwt::LMTStatus::CH) {
            if(!include_possibility && it->status == mwt::LMTStatus::Possible) {
                throw std::logic_error("There should not be any possible edges left!");
            }
            if(it->status == mwt::LMTStatus::Impossible) {
                continue;
            }
        }
        if(!it->is_primary()) {
            continue;
        }
        if(first) {
            first = false;
        } else {
            output << ",\n     ";
        }
        output << "[" << (it->source_handle() - &*points_begin) << ", " << (it->target_handle() - &*points_begin);
        if(include_possibility) {
            output << ", " << int(it->status);
        }
        output << "]";
    }
    output << "\n    ]\n";
    output << "}\n";
}

} // namespace mwt

#endif
