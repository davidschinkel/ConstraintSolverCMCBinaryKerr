#include "main.h"

//float
template<>
float cast<>(const std::string in)
{
    return atof(in.c_str());
}

template<>
float cast<>(char* in)
{
    return atof(in);
}
//double
template<>
double cast<>(const std::string in)
{
    return atof(in.c_str());
}

template<>
double cast<>(char* in)
{
    return atof(in);
}
//long double
template<>
long double cast<>(const std::string in)
{
    return strtold(in.c_str(),NULL);
}

template<>
long double cast<>(char* in)
{
    return strtold(in,NULL);
}
