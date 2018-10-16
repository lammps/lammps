#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "ptm_constants.h"
#include "ptm_voronoi_cell.h"
#include "ptm_neighbour_ordering.h"


namespace ptm {

typedef struct
{
	double area;
	double dist;
	int index;
} sorthelper_t;

static bool sorthelper_compare(sorthelper_t const& a, sorthelper_t const& b)
{
	if (a.area > b.area)
		return true;

	if (a.area < b.area)
		return false;

	if (a.dist < b.dist)
		return true;

	return false;
}

//todo: change voronoi code to return errors rather than exiting
static int calculate_voronoi_face_areas(int num_points, const double (*_points)[3], double* normsq, double max_norm, ptm_voro::voronoicell_neighbor* v, std::vector<int>& nbr_indices, std::vector<double>& face_areas)
{
	const double k = 1000 * max_norm;	//todo: reduce this constant
	v->init(-k,k,-k,k,-k,k);

	for (int i=1;i<num_points;i++)
	{
		double x = _points[i][0] - _points[0][0];
		double y = _points[i][1] - _points[0][1];
		double z = _points[i][2] - _points[0][2];
		v->nplane(x,y,z,normsq[i],i);
	}

	v->neighbors(nbr_indices);
	v->face_areas(face_areas);
	return 0;
}

int calculate_neighbour_ordering(void* _voronoi_handle, int num_points, const double (*_points)[3], int8_t* ordering)
{
	assert(num_points <= PTM_MAX_INPUT_POINTS);

	ptm_voro::voronoicell_neighbor* voronoi_handle = (ptm_voro::voronoicell_neighbor*)_voronoi_handle;

	double max_norm = 0;
	double points[PTM_MAX_INPUT_POINTS][3];
	double normsq[PTM_MAX_INPUT_POINTS];
	for (int i = 0;i<num_points;i++)
	{
		double x = _points[i][0] - _points[0][0];
		double y = _points[i][1] - _points[0][1];
		double z = _points[i][2] - _points[0][2];
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;

		normsq[i] = x*x + y*y + z*z;
		max_norm = std::max(max_norm, normsq[i]);
#ifdef DEBUG
		printf("point %d: %f\t%f\t%f\t%f\n", i, x, y, z, x*x + y*y + z*z);
#endif
	}

	max_norm = sqrt(max_norm);

	std::vector<int> nbr_indices(num_points + 6);
	std::vector<double> face_areas(num_points + 6);
	int ret = calculate_voronoi_face_areas(num_points, points, normsq, max_norm, voronoi_handle, nbr_indices, face_areas);
	if (ret != 0)
		return ret;

	double areas[PTM_MAX_INPUT_POINTS];
	memset(areas, 0, num_points * sizeof(double));
	areas[0] = INFINITY;
	for (size_t i=0;i<nbr_indices.size();i++)
	{
		int index = nbr_indices[i];
		if (index > 0)
			areas[index] = face_areas[i];
	}

	sorthelper_t data[PTM_MAX_INPUT_POINTS];
	for (int i=0;i<num_points;i++)
	{
		assert(areas[i] == areas[i]);
		data[i].area = areas[i];
		data[i].dist = normsq[i];
		data[i].index = i;
	}

	std::sort(data, data + num_points, &sorthelper_compare);

#ifdef DEBUG
	for (int i=0;i<num_points;i++)
		printf("%d %f\n", data[i].index, data[i].area);
#endif

	for (int i=0;i<num_points;i++)
		ordering[i] = data[i].index;

	return ret;
}

void* voronoi_initialize_local()
{
	ptm_voro::voronoicell_neighbor* ptr = new ptm_voro::voronoicell_neighbor;
	return (void*)ptr;
}

void voronoi_uninitialize_local(void* _ptr)
{
	ptm_voro::voronoicell_neighbor* ptr = (ptm_voro::voronoicell_neighbor*)_ptr;
	delete ptr;
}


typedef struct
{
	double dist;
	int p;
	int index;
} diamond_t;

static bool diamond_compare(diamond_t const& a, diamond_t const& b)
{
	return a.dist < b.dist;
}

int calculate_diamond_neighbour_ordering(	int num_points, double (*unpermuted_points)[3], int32_t* unpermuted_numbers,
						int8_t* ordering, double (*points)[3], int32_t* numbers)
{
	assert(num_points <= PTM_MAX_INPUT_POINTS);

	diamond_t data[4 * (PTM_MAX_INPUT_POINTS - 5)];
	int index = 0;
	for (int i=5;i<num_points;i++)
	{
		for (int j=1;j<5;j++)
		{
			double dx = unpermuted_points[i][0] - unpermuted_points[j][0];
			double dy = unpermuted_points[i][1] - unpermuted_points[j][1];
			double dz = unpermuted_points[i][2] - unpermuted_points[j][2];

			double d = dx*dx + dy*dy + dz*dz;

			data[index].p = j - 1;
			data[index].index = i;
			data[index].dist = d;
			index++;
		}
	}
	int n = index;

	std::sort(data, data + n, &diamond_compare);

	for (index=0;index<5;index++)
		ordering[index] = index;

	int num_found = 0;
	bool hit[PTM_MAX_INPUT_POINTS] = {0};
	int counts[4] = {0};
	for (int i=0;i<n;i++)
	{
		int p = data[i].p;
		int q = data[i].index;
		if (hit[q] || counts[p] >= 3)
			continue;

		ordering[1 + 4 + 3 * p + counts[p]] = q;
		counts[p]++;
		index++;
		num_found++;
		if (num_found >= 12)
			break;
	}

	if (num_found != 12)
		return -1;

	for (int i=0;i<PTM_NUM_NBRS_DCUB+1;i++)
	{
		memcpy(points[i], &unpermuted_points[ordering[i]], 3 * sizeof(double));

		if (unpermuted_numbers != NULL)
			numbers[i] = unpermuted_numbers[ordering[i]];
	}

	return 0;
}

}

