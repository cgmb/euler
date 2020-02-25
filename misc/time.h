#pragma once

#include <time.h>

struct timespec wait_until(struct timespec start, long nsec_after_start);
void get_current_time(struct timespec* time);
