#ifndef SINCOS_LOOKUP_H
#define SINCOS_LOOKUP_H

namespace LAMMPS_NS {

namespace user_manifold {

class sincos_lookup
{
public:
	sincos_lookup( int N );
	~sincos_lookup();
	double sin(double x);
	double cos(double x);
	void sincos(double x, double *s, double *c);

	void test_lookup();


private:
	double *cos_table;
	double dx_tab, inv_dx_tab;
	int size;

	void make_cos_table();
};

}

}

#endif // SINCOS_LOOKUP_H
