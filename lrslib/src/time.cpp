extern "C" {
#include "sys/time.h"
}

void gettimeofday(struct timeval *tv, void*) {
    if (tv) {
        struct timespec ts;
        timespec_get(&ts, TIME_UTC);
        tv->tv_sec = static_cast<long>(ts.tv_sec);
        tv->tv_usec = static_cast<long>(ts.tv_nsec / 1000);
    }
}