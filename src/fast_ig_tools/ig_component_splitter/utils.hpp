#pragma once

#include <perfcounter.hpp>

#if __DARWIN || __DARWIN_UNIX03
#include <mach/task.h>
#include <mach/mach.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#if __DARWIN || __DARWIN_UNIX03
inline size_t get_max_rss() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS !=
      task_info(mach_task_self(),
                TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    return -1U;

  return t_info.resident_size / 1024;
}
#else
inline size_t get_max_rss() {
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);

  return ru.ru_maxrss;
}
#endif

std::string join_cmd_line(size_t argc, char **argv);

void create_console_logger(std::string log_props_file = "");

std::string running_time_format(const perf_counter &pc);
