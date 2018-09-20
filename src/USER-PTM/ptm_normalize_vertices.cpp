#include <cmath>

namespace ptm {

void subtract_barycentre(int num, double (*points)[3], double (*normalized)[3])
{
	//calculate barycentre
	double sum[3] = {0, 0, 0};
	for (int i=0;i<num;i++)
	{
		sum[0] += points[i][0];
		sum[1] += points[i][1];
		sum[2] += points[i][2];
	}

	sum[0] /= num;
	sum[1] /= num;
	sum[2] /= num;

	//subtract barycentre
	for (int i=0;i<num;i++)
	{
		normalized[i][0] = points[i][0] - sum[0];
		normalized[i][1] = points[i][1] - sum[1];
		normalized[i][2] = points[i][2] - sum[2];
	}
}

double normalize_vertices(int num, double (*points)[3], double (*normalized)[3])
{
	subtract_barycentre(num, points, normalized);

	//calculate mean length
	double scale = 0.0;
	for (int i=1;i<num;i++)
	{
		double x = normalized[i][0];
		double y = normalized[i][1];
		double z = normalized[i][2];

		double norm = sqrt(x*x + y*y + z*z);
		scale += norm;
	}
	scale /= num;

	//scale vertices such that mean length is 1
	for (int i=0;i<num;i++)
	{
		normalized[i][0] /= scale;
		normalized[i][1] /= scale;
		normalized[i][2] /= scale;
	}

	return scale;
}

}

