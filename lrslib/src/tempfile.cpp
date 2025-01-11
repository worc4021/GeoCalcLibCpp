#include "tempfile.h"
#include <io.h>

char *mktemp(char *t){
    return  _mktemp(t);
}