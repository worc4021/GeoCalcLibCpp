#include <string>

extern "C" {
    char *basename(char *path);
    char *dirname(char *path);
}

std::string
pathname_directory(const std::string &pathname)
{
    std::size_t len = pathname.find_last_of("/\\");
    return len == std::string::npos ? "" : pathname.substr(0,len);
}

std::string
pathname_sans_directory(const std::string &pathname)
{
    std::size_t len = pathname.find_last_of("/\\");
    return len == std::string::npos ? pathname : pathname.substr(len + 1);
}

char *basename(char *path){
    std::string pathname(path);
    std::string basename = pathname_sans_directory(pathname);
    char *cstr = new char[basename.length() + 1];
    strcpy_s(cstr, basename.length() + 1,basename.c_str());
    return cstr;
}

char *dirname(char *path){
    std::string pathname(path);
    std::string dirname = pathname_directory(pathname);
    char *cstr = new char[dirname.length() + 1];
    strcpy_s(cstr, dirname.length() + 1, dirname.c_str());
    return cstr;
}