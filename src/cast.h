#ifndef CAST_H
#define CAST_H

template <typename X, typename Y>
X cast(Y in)
{
    return X(in);
}
//float
template<>
float cast<>(const std::string in);
template<>
float cast<>(char* in);
//double
template<>
double cast<>(const std::string in);
template<>
double cast<>(char* in);
//long double
template<>
long double cast<>(const std::string in);
template<>
long double cast<>(char* in);

#endif
