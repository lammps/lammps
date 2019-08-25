# ifndef ERF_H
# define ERF_H

# ifdef _WIN32

# ifdef __cplusplus
extern "C" {
# endif

double erf(double x);
double erfc(double x);

# ifdef __cplusplus
}
# endif

# endif

# endif
