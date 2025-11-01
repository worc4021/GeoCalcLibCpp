#ifndef SYS_TIME_H
#define SYS_TIME_H
#include <time.h>

struct timeval {
    time_t t;
    long  tv_sec;
    long  tv_usec;
};
void gettimeofday(struct timeval *tv, void*);
#endif /* SYS_TIME_H */