#include "misc/time.h"
#include <unistd.h>

static struct timespec difference(struct timespec lhs, struct timespec rhs) {
  struct timespec diff;
  diff.tv_sec = lhs.tv_sec - rhs.tv_sec;
  diff.tv_nsec = lhs.tv_nsec - rhs.tv_nsec;
  if (lhs.tv_nsec < rhs.tv_nsec) {
    diff.tv_sec -= 1;
    diff.tv_nsec += 1e9;
  }
  return diff;
}

// wait for up to one second from the given start time
// returns the current time when it exits
struct timespec wait_until(struct timespec start, long nsec_after_start) {
  struct timespec now;
  get_current_time(&now);

  struct timespec already_elapsed = difference(now, start);
  if (already_elapsed.tv_sec == 0) {
    long wait_for = (nsec_after_start - already_elapsed.tv_nsec) / 1000;
    if (wait_for > 0) {
      usleep(wait_for);
      get_current_time(&now);
    }
  }

  return now;
}

void get_current_time(struct timespec* time) {
  clock_gettime(CLOCK_MONOTONIC, time);
}
