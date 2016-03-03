#ifndef SINCOS_LOOKUP_H
#define SINCOS_LOOKUP_H

class sincos_lookup
{
public:
	sincos_lookup( int N );
	~sincos_lookup();
	double sin(double x);
	double cos(double x);
	void sincos(double x, double *s, double *c);

	void test_lookup();

#if __cplusplus < 201103L
        #define my_pi2 1.57079632679489661923
        #define my_pi  3.14159265358979323846
        #define my_2pi 6.28318530717958647692
#else
	static constexpr const double my_pi2 = 1.57079632679489661923;
	static constexpr const double my_pi  = 3.14159265358979323846;
	static constexpr const double my_2pi = 6.28318530717958647692;
#endif

private:
	double *cos_table;
	double dx_tab, inv_dx_tab;
	int size;



	void make_cos_table();
};

#endif // SINCOS_LOOKUP_H
