#ifndef CGAL_MWT_READ_INSTANCE_H_INCLUDED_
#define CGAL_MWT_READ_INSTANCE_H_INCLUDED_

#include <boost/exception/diagnostic_information.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <optional>
#include <regex>
#include <string>
#include <variant>

namespace mwt {

/**
 * Turn a JSON number (integer or floating point) into a CGAL number.
 */
template<typename Kernel> typename Kernel::FT json_number_to_ft(const nlohmann::json &json_number) {
    if(json_number.is_number_integer()) {
        std::int64_t i = json_number.get<std::int64_t>();
        return typename Kernel::FT(i);
    } else {
        double d = json_number.get<double>();
        return typename Kernel::FT(d);
    }
}

/**
 * Try to interpret a JSON array [x,y] as point.
 */
template<typename Kernel> typename Kernel::Point_2 interpret_json_array2(const nlohmann::json &json_array) {
    if(json_array.size() != 2) {
        throw std::runtime_error("Invalid JSON data: JSON array of length != 2 for point!");
    }
    if(!json_array[0].is_number() || !json_array[1].is_number()) {
        throw std::runtime_error("Invalid JSON data: JSON array contains non-numerical entries!");
    }
    return typename Kernel::Point_2(json_number_to_ft<Kernel>(json_array[0]), json_number_to_ft<Kernel>(json_array[1]));
}

/**
 * Try to interpret a JSON object with {"x": x, "y": y} as point.
 */
template<typename Kernel> typename Kernel::Point_2 interpret_json_object_xy(const nlohmann::json &json_object) {
    if(!json_object.count("x") || !json_object.count("y")) {
        throw std::runtime_error("Invalid JSON data: JSON object does not contain x and y!");
    }
    const auto &x = json_object["x"];
    const auto &y = json_object["y"];
    if(!x.is_number() || !y.is_number()) {
        throw std::runtime_error("Invalid JSON data: JSON object 'x' or 'y' contains non-numerical entries!");
    }
    return typename Kernel::Point_2(json_number_to_ft<Kernel>(x), json_number_to_ft<Kernel>(y));
}

template<typename PointType> std::vector<PointType> unique_points(std::vector<PointType> &&points) {
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    return std::vector<PointType>{std::move(points)};
}

template<typename Kernel> std::vector<typename Kernel::Point_2> interpret_json(const nlohmann::json &jsondata) {
    std::vector<typename Kernel::Point_2> points;
    if(jsondata.is_array()) {
        // assume an array of points
        if(jsondata.size() == 0) {
            return points;
        }
        // check how the points look
        const auto &first = jsondata[0];
        if(first.is_array()) {
            // assume each point is [x, y]
            for(const auto &point : jsondata) {
                points.emplace_back(interpret_json_array2<Kernel>(point));
            }
            return points;
        } else if(first.is_object()) {
            // assume each point has { "x": x, "y": y }
            for(const auto &point : jsondata) {
                points.emplace_back(interpret_json_object_xy<Kernel>(point));
            }
        } else {
            throw std::runtime_error("Invalid JSON data: JSON array contains neither arrays nor objects!");
        }
    } else if(jsondata.is_object()) {
        // assume an object containing a field with points
        const char *point_field_names[] = {"points", "Points", "POINTS"};
        for(const char *fname : point_field_names) {
            if(jsondata.count(fname)) {
                return interpret_json<Kernel>(jsondata[fname]);
            }
        }
        throw std::runtime_error("Invalid JSON data: JSON object does not appear to contain a field with points!");
    } else {
        throw std::runtime_error("Invalid JSON data: JSON data is neither an array nor an object!");
    }
    return points;
}

inline bool is_probably_gzip_compressed(std::istream &input) {
    input.seekg(0, std::ios::beg);
    char buf[2];
    input.read(buf, 2);
    bool result = (buf[0] == '\x1f' && buf[1] == '\x8b');
    input.seekg(0, std::ios::beg);
    return result;
}

inline bool is_probably_bzip2_compressed(std::istream &input) {
    input.seekg(0, std::ios::beg);
    char buf[3];
    input.read(buf, 3);
    bool result = (buf[0] == 'B' && buf[1] == 'Z' && buf[2] == 'h');
    input.seekg(0, std::ios::beg);
    return result;
}

template<typename Kernel>
std::variant<std::string, std::vector<typename Kernel::Point_2>> read_instance_json(std::istream &input) {
    try {
        nlohmann::json jsondata = nlohmann::json::parse(input);
        return unique_points(interpret_json<Kernel>(jsondata));
    } catch(const std::exception &json_err) {
        return std::string(json_err.what());
    }
}

template<typename Kernel>
bool parse_triple_line(const std::string &line, typename Kernel::FT &x, typename Kernel::FT &y) {
    /**
     * Ugly workaround for libc++ bug.
     * https://github.com/llvm/llvm-project/issues/38360
     */
    using FT = typename Kernel::FT;
    if constexpr(std::is_same_v<std::decay_t<FT>, double>) {
        std::istringstream line_parser(line);
        std::string idx, x_str, y_str;
        if(!(line_parser >> idx >> x_str >> y_str)) {
            return false;
        }
        char *end = nullptr;
        x = std::strtod(x_str.c_str(), &end);
        if(*end != '\0') {
            return false;
        }
        y = std::strtod(y_str.c_str(), &end);
        if(*end != '\0') {
            return false;
        }
        return true;
    } else {
        std::istringstream line_parser(line);
        typename Kernel::FT idx;
        return bool(line_parser >> idx >> x >> y);
    }
}

template<typename Kernel>
bool parse_pair_line(const std::string &line, typename Kernel::FT &x, typename Kernel::FT &y) {
    /**
     * Ugly workaround for libc++ bug.
     * https://github.com/llvm/llvm-project/issues/38360
     */
    using FT = typename Kernel::FT;
    if constexpr(std::is_same_v<std::decay_t<FT>, double>) {
        std::istringstream line_parser(line);
        std::string x_str, y_str;
        if(!(line_parser >> x_str >> y_str)) {
            return false;
        }
        char *end = nullptr;
        x = std::strtod(x_str.c_str(), &end);
        if(*end != '\0') {
            return false;
        }
        y = std::strtod(y_str.c_str(), &end);
        if(*end != '\0') {
            return false;
        }
        return true;
    } else {
        std::istringstream line_parser(line);
        return bool(line_parser >> x >> y);
    }
}

template<typename Kernel> std::vector<typename Kernel::Point_2> read_instance_text(std::istream &input) {
    std::vector<typename Kernel::Point_2> points;
    input.exceptions(std::ios::badbit);
    std::string line;
    std::optional<bool> is_triple_lines;
    while(std::getline(input, line)) {
        if(line.empty()) {
            continue;
        }
        std::size_t first_nonspace = line.find_first_not_of(" \t\r\n");
        if(first_nonspace == std::string::npos || line[first_nonspace] == '#') {
            continue;
        }
        if(line[first_nonspace] == 'E' || line[first_nonspace] == 'e') {
            if(line.find("EOF", first_nonspace) == first_nonspace ||
               line.find("eof", first_nonspace) == first_nonspace ||
               line.find("Eof", first_nonspace) == first_nonspace) {
                // support some EOF markers from TSPLIB files
                break;
            }
        }
        typename Kernel::FT x, y;
        if(!is_triple_lines) {
            is_triple_lines = parse_triple_line<Kernel>(line, x, y);
        }
        if(*is_triple_lines) {
            if(!parse_triple_line<Kernel>(line, x, y)) {
                throw std::runtime_error("Invalid text file: could not parse line '" + line + "'.\n");
            }
        } else {
            if(!parse_pair_line<Kernel>(line, x, y)) {
                throw std::runtime_error("Invalid text file: could not parse line '" + line + "'.\n");
            }
        }
        points.emplace_back(x, y);
    }
    return unique_points(std::move(points));
}

template<typename Kernel, int CompressionMethod>
std::vector<typename Kernel::Point_2> read_compressed_instance(std::istream &input) {
    // first try JSON; we cannot seekg in compressed streams,
    // but we can seekg on the original input stream
    std::string json_error;

    auto setup_instream = [&](boost::iostreams::filtering_istream &in) {
        if constexpr(CompressionMethod == 0) {
            in.push(boost::iostreams::gzip_decompressor());
        } else if constexpr(CompressionMethod == 1) {
            in.push(boost::iostreams::bzip2_decompressor());
        } else {
            static_assert(CompressionMethod == 0, "Invalid/unknown compression method");
        }
        in.push(input);
        in.set_auto_close(false);
    };

    {
        // first try reading as compressed JSON
        boost::iostreams::filtering_istream in;
        setup_instream(in);
        auto json_result = read_instance_json<Kernel>(in);
        if(std::holds_alternative<std::string>(json_result)) {
            json_error = std::get<std::string>(json_result);
        } else {
            return std::move(std::get<std::vector<typename Kernel::Point_2>>(json_result));
        }
    }

    input.seekg(0, std::ios::beg);

    {
        // then try reading as compressed text
        boost::iostreams::filtering_istream in;
        setup_instream(in);
        try {
            return read_instance_text<Kernel>(in);
        } catch(std::exception &ex) {
            throw std::runtime_error("Reading as text file failed: " + std::string(ex.what()) +
                                     "\nParsing as JSON also failed: " + json_error);
        }
    }
}

template<typename Kernel> std::vector<typename Kernel::Point_2> read_instance(std::istream &input) {
    if(is_probably_gzip_compressed(input)) {
        return read_compressed_instance<Kernel, 0>(input);
    }

    if(is_probably_bzip2_compressed(input)) {
        return read_compressed_instance<Kernel, 1>(input);
    }

    auto json_result = read_instance_json<Kernel>(input);
    if(std::holds_alternative<std::string>(json_result)) {
        try {
            input.seekg(0, std::ios::beg);
            return read_instance_text<Kernel>(input);
        } catch(const std::runtime_error &ex) {
            throw std::runtime_error(std::string("Reading as text file failed: ") + ex.what() +
                                     "\nParsing as JSON also failed: " + std::get<std::string>(json_result));
        }
    } else {
        return std::move(std::get<std::vector<typename Kernel::Point_2>>(json_result));
    }
}

template<typename Kernel> std::vector<typename Kernel::Point_2> read_instance(const std::string &filename) {
    std::fstream input;
    input.exceptions(std::ios::failbit | std::ios::badbit);
    input.open(filename.c_str(), std::ios::in | std::ios::binary);
    return read_instance<Kernel>(input);
}

template<typename Kernel>
std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>>
solution_from_json_result(const nlohmann::json &json_result) {
    std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>> result;
    std::vector<typename Kernel::Point_2> &points = result.first;
    std::vector<std::array<std::uint64_t, 2>> &edges = result.second;
    for(const auto &point : json_result["points"]) {
        points.emplace_back(interpret_json_array2<Kernel>(point));
    }
    edges = json_result["edges"].get<std::vector<std::array<std::uint64_t, 2>>>();
    return result;
}

template<typename Kernel>
std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>>
read_solution(std::istream &input) {
    if(is_probably_gzip_compressed(input)) {
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(input);
        return read_solution<Kernel>(in);
    }

    if(is_probably_bzip2_compressed(input)) {
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::bzip2_decompressor());
        in.push(input);
        return read_solution<Kernel>(in);
    }

    auto json_result = nlohmann::json::parse(input);
    return solution_from_json_result<Kernel>(json_result);
}

template<typename Kernel>
std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>>
transform_text_solution(const std::vector<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>> &solution) {
    using Point = typename Kernel::Point_2;
    std::pair<std::vector<Point>, std::vector<std::array<std::uint64_t, 2>>> result;
    std::vector<Point> &points = result.first;
    std::vector<std::array<std::uint64_t, 2>> &edges = result.second;
    std::vector<Point> all_points;
    for(const auto &edge : solution) {
        all_points.emplace_back(edge.first);
        all_points.emplace_back(edge.second);
    }
    std::sort(all_points.begin(), all_points.end());
    all_points.erase(std::unique(all_points.begin(), all_points.end()), all_points.end());
    auto point_index = [&all_points](const Point &p) {
        return std::size_t(std::lower_bound(all_points.begin(), all_points.end(), p) - all_points.begin());
    };
    for(const auto &edge : solution) {
        edges.push_back({point_index(edge.first), point_index(edge.second)});
    }
    points = std::move(all_points);
    points.shrink_to_fit();
    edges.shrink_to_fit();
    return result;
}

template<typename Kernel>
std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>>
read_text_solution(std::istream &input) {
    using Point = typename Kernel::Point_2;
    using PPair = std::pair<Point, Point>;
    if(is_probably_gzip_compressed(input)) {
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(input);
        return read_text_solution<Kernel>(in);
    }

    if(is_probably_bzip2_compressed(input)) {
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::bzip2_decompressor());
        in.push(input);
        return read_text_solution<Kernel>(in);
    }

    std::regex line_split_regex("^\\s*[(]([^()]+)[)]\\s+[(]([^()]+)[)]\\s*$");
    std::string line;
    std::vector<PPair> result;
    input.exceptions(std::ios::badbit);
    while(std::getline(input, line)) {
        if(line.empty()) {
            continue;
        }
        std::smatch match;
        if(!std::regex_match(line, match, line_split_regex)) {
            throw std::runtime_error("Invalid text file: could not parse line '" + line + "'.\n");
        }
        std::istringstream input1(match.str(1));
        std::istringstream input2(match.str(2));
        PPair ps;
        if(!(input1 >> ps.first) || !(input2 >> ps.second)) {
            throw std::runtime_error("Invalid text file: could not parse line '" + line + "'.\n");
        }
        result.emplace_back(std::move(ps));
    }
    return transform_text_solution<Kernel>(result);
}

template<typename Kernel>
std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>>
read_solution(const std::string &filename) {
    std::fstream input;
    input.exceptions(std::ios::failbit | std::ios::badbit);
    input.open(filename.c_str(), std::ios::in | std::ios::binary);
    return read_solution<Kernel>(input);
}

template<typename Kernel>
std::pair<std::vector<typename Kernel::Point_2>, std::vector<std::array<std::uint64_t, 2>>>
read_text_solution(const std::string &filename) {
    std::fstream input;
    input.exceptions(std::ios::failbit | std::ios::badbit);
    input.open(filename.c_str(), std::ios::in | std::ios::binary);
    return read_text_solution<Kernel>(input);
}

} // namespace mwt

#endif
