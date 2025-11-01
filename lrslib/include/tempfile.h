#ifndef TEMPFILE_H
#define TEMPFILE_H

char *mktemp(char *t);
int mkstemp(char *t);
int close(int fd);

#endif 