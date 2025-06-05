#ifndef CGAL_MWT_RESOURCE_MONITOR_H_INCLUDED_
#define CGAL_MWT_RESOURCE_MONITOR_H_INCLUDED_

#include <condition_variable>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>

namespace mwt {

struct CurrentUsage {
    double process_user_time = -1.0;
    double process_system_time = -1.0;
    std::int64_t max_resident_set_size = -1;
    std::int64_t resident_set_size = -1;
    std::int64_t minor_page_faults = -1;
    std::int64_t major_page_faults = -1;
    std::int64_t voluntary_context_switches = -1;
    std::int64_t involuntary_context_switches = -1;
    std::int64_t swap_usage = -1;

    void update_maxima(const CurrentUsage &other) {
        process_user_time = std::max(process_user_time, other.process_user_time);
        process_system_time = std::max(process_system_time, other.process_system_time);
        max_resident_set_size = std::max(max_resident_set_size, other.max_resident_set_size);
        resident_set_size = std::max(resident_set_size, other.resident_set_size);
        max_resident_set_size = std::max(max_resident_set_size, resident_set_size);
        minor_page_faults = std::max(minor_page_faults, other.minor_page_faults);
        major_page_faults = std::max(major_page_faults, other.major_page_faults);
        voluntary_context_switches = std::max(voluntary_context_switches, other.voluntary_context_switches);
        involuntary_context_switches = std::max(involuntary_context_switches, other.involuntary_context_switches);
        swap_usage = std::max(swap_usage, other.swap_usage);
    }

    nlohmann::json to_json() const {
        nlohmann::json result;
        result["process_user_time"] = process_user_time;
        result["process_system_time"] = process_system_time;
        result["max_resident_set_size"] = max_resident_set_size;
        result["resident_set_size"] = resident_set_size;
        result["minor_page_faults"] = minor_page_faults;
        result["major_page_faults"] = major_page_faults;
        result["voluntary_context_switches"] = voluntary_context_switches;
        result["involuntary_context_switches"] = involuntary_context_switches;
        result["swap_usage"] = swap_usage;
        return result;
    }
};

} // namespace mwt

#if __has_include(<sys/resource.h>)

#include <sys/resource.h>

namespace mwt {

inline void update_with_rusage(CurrentUsage &usage) {
    struct rusage rusage;
    if(getrusage(RUSAGE_SELF, &rusage) == 0) {
        usage.process_user_time = double(rusage.ru_utime.tv_sec) + rusage.ru_utime.tv_usec * 1.0e-6;
        usage.process_system_time = double(rusage.ru_stime.tv_sec) + rusage.ru_stime.tv_usec * 1.0e-6;
        usage.max_resident_set_size = rusage.ru_maxrss;
#if defined(__linux__)
#if __linux__
        usage.max_resident_set_size *= 1024;
#endif
#endif
        usage.minor_page_faults = rusage.ru_minflt;
        usage.major_page_faults = rusage.ru_majflt;
        usage.voluntary_context_switches = rusage.ru_nvcsw;
        usage.involuntary_context_switches = rusage.ru_nivcsw;
    }
}

inline void update_with_proc(CurrentUsage &usage) {
    std::ifstream read_status("/proc/self/status", std::ios::in);
    if(!read_status)
        return;
    std::string line;
    auto kb_line_into_value = [](const std::string &line, std::int64_t &value) {
        std::istringstream iss(line);
        std::string name;
        iss >> name >> value;
        value *= 1024;
    };
    while(std::getline(read_status, line)) {
        if(line.rfind("VmRSS:", 0) == 0) {
            kb_line_into_value(line, usage.resident_set_size);
        }
        if(line.rfind("VmSwap:", 0) == 0) {
            kb_line_into_value(line, usage.swap_usage);
        }
    }
}

static constexpr bool CGAL_MWT_RESOURCE_MONITOR_SUPPORTED = true;

} // namespace mwt

#else

namespace mwt {

void update_with_rusage(CurrentUsage &usage) {}
void update_with_proc(CurrentUsage &usage) {}

static constexpr bool CGAL_MWT_RESOURCE_MONITOR_SUPPORTED = false;

} // namespace mwt

#endif

namespace mwt {

inline CurrentUsage get_current_usage() {
    CurrentUsage usage;
    if(CGAL_MWT_RESOURCE_MONITOR_SUPPORTED) {
        update_with_rusage(usage);
        update_with_proc(usage);
    }
    return usage;
}

class ResourceMonitor {
  public:
    ~ResourceMonitor() { stop(); }

    void start() {
        if(!CGAL_MWT_RESOURCE_MONITOR_SUPPORTED)
            return;
        if(res_mon_thread.joinable())
            return;
        res_mon_should_stop = false;
        res_mon_thread = std::thread([this]() { this->p_res_mon_thread_main(); });
    }

    void stop() {
        if(!res_mon_thread.joinable())
            return;
        {
            std::unique_lock<std::mutex> lock(res_mon_lock);
            res_mon_should_stop = true;
        }
        res_cond_var.notify_all();
        res_mon_thread.join();
        res_max_usage.update_maxima(get_current_usage());
    }

    CurrentUsage get_max_usage() const {
        CurrentUsage result;
        {
            std::unique_lock<std::mutex> lock(res_mon_lock);
            result = res_max_usage;
        }
        return result;
    }

  private:
    void p_res_mon_thread_main() {
        std::unique_lock<std::mutex> lock(res_mon_lock);
        while(!res_mon_should_stop) {
            res_max_usage.update_maxima(get_current_usage());
            res_cond_var.wait_for(lock, std::chrono::duration<double>(res_mon_sleep_interval));
        }
    }

    mutable std::mutex res_mon_lock;
    std::condition_variable res_cond_var;
    std::thread res_mon_thread;
    bool res_mon_should_stop = false;
    double res_mon_sleep_interval = 0.25;
    CurrentUsage res_max_usage;
};

} // namespace mwt

#endif
