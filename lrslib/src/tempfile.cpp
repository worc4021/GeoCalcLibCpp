extern "C" {
#include "tempfile.h"
}

#include <stdio.h>
#include <io.h>
#include <string.h>

extern "C" char *mktemp(char *t){
    return  _mktemp(t);
}

extern "C" int mkstemp(char *t){
    int err = _mktemp_s(t, strlen(t) + 1);
    if (err != 0) {
        return -1;
    }

    FILE *fp;
    fopen_s(&fp, t, "w+");
    if (fp == NULL) {
        return -1;
    }

    int fd = _fileno(fp);
    return fd;
}

extern "C" int close(int fd){
    return _close(fd);
}