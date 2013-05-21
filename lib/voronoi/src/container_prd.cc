// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file container_prd.cc
 * \brief Function implementations for the container_periodic_base and
 * related classes. */

#include "container_prd.hh"

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (bx_) The x coordinate of the first unit vector.
 * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
 * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
 *                            vector.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] init_mem_ the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_periodic_base::container_periodic_base(double bx_,double bxy_,double by_,
		double bxz_,double byz_,double bz_,int nx_,int ny_,int nz_,int init_mem_,int ps_)
	: unitcell(bx_,bxy_,by_,bxz_,byz_,bz_), voro_base(nx_,ny_,nz_,bx_/nx_,by_/ny_,bz_/nz_),
	ey(int(max_uv_y*ysp+1)), ez(int(max_uv_z*zsp+1)), wy(ny+ey), wz(nz+ez),
	oy(ny+2*ey), oz(nz+2*ez), oxyz(nx*oy*oz), id(new int*[oxyz]), p(new double*[oxyz]),
	co(new int[oxyz]), mem(new int[oxyz]), img(new char[oxyz]), init_mem(init_mem_), ps(ps_) {
	int i,j,k,l;

	// Clear the global arrays
	int *pp=co;while(pp<co+oxyz) *(pp++)=0;
	pp=mem;while(pp<mem+oxyz) *(pp++)=0;
	char *cp=img;while(cp<img+oxyz) *(cp++)=0;

	// Set up memory for the blocks in the primary domain
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		mem[l]=init_mem;
		id[l]=new int[init_mem];
		p[l]=new double[ps*init_mem];
	}
}

/** The container destructor frees the dynamically allocated memory. */
container_periodic_base::~container_periodic_base() {
	for(int l=oxyz-1;l>=0;l--) if(mem[l]>0) {
		delete [] p[l];
		delete [] id[l];
	}
	delete [] img;
	delete [] mem;
	delete [] co;
	delete [] id;
	delete [] p;
}

/** The class constructor sets up the geometry of container.
 * \param[in] (bx_) The x coordinate of the first unit vector.
 * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
 * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
 *                            vector.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *			    coordinate directions.
 * \param[in] init_mem_ the initial memory allocation for each block. */
container_periodic::container_periodic(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
	int nx_,int ny_,int nz_,int init_mem_)
	: container_periodic_base(bx_,bxy_,by_,bxz_,byz_,bz_,nx_,ny_,nz_,init_mem_,3),
	vc(*this,2*nx_+1,2*ey+1,2*ez+1) {}

/** The class constructor sets up the geometry of container.
 * \param[in] (bx_) The x coordinate of the first unit vector.
 * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
 * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
 *                            vector.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *			    coordinate directions.
 * \param[in] init_mem_ the initial memory allocation for each block. */
container_periodic_poly::container_periodic_poly(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
	int nx_,int ny_,int nz_,int init_mem_)
	: container_periodic_base(bx_,bxy_,by_,bxz_,byz_,bz_,nx_,ny_,nz_,init_mem_,4),
	vc(*this,2*nx_+1,2*ey+1,2*ez+1) {ppr=p;}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_periodic::put(int n,double x,double y,double z) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	double *pp=p[ijk]+3*co[ijk]++;
	*(pp++)=x;*(pp++)=y;*pp=z;
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_periodic_poly::put(int n,double x,double y,double z,double r) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	double *pp=p[ijk]+4*co[ijk]++;
	*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
	if(max_radius<r) max_radius=r;
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
 * 			  in, with (0,0,0) corresponding to the primary domain.
 */
void container_periodic::put(int n,double x,double y,double z,int &ai,int &aj,int &ak) {
	int ijk;
	put_locate_block(ijk,x,y,z,ai,aj,ak);
	id[ijk][co[ijk]]=n;
	double *pp=p[ijk]+3*co[ijk]++;
	*(pp++)=x;*(pp++)=y;*pp=z;
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle.
 * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
 * 			  in, with (0,0,0) corresponding to the primary domain.
 */
void container_periodic_poly::put(int n,double x,double y,double z,double r,int &ai,int &aj,int &ak) {
	int ijk;
	put_locate_block(ijk,x,y,z,ai,aj,ak);
	id[ijk][co[ijk]]=n;
	double *pp=p[ijk]+4*co[ijk]++;
	*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
	if(max_radius<r) max_radius=r;
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_periodic::put(particle_order &vo,int n,double x,double y,double z) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	vo.add(ijk,co[ijk]);
	double *pp=p[ijk]+3*co[ijk]++;
	*(pp++)=x;*(pp++)=y;*pp=z;
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_periodic_poly::put(particle_order &vo,int n,double x,double y,double z,double r) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	vo.add(ijk,co[ijk]);
	double *pp=p[ijk]+4*co[ijk]++;
	*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
	if(max_radius<r) max_radius=r;
}

/** Takes a particle position vector and computes the region index into which
 * it should be stored. If the container is periodic, then the routine also
 * maps the particle position to ensure it is in the primary domain. If the
 * container is not periodic, the routine bails out.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
void container_periodic_base::put_locate_block(int &ijk,double &x,double &y,double &z) {

	// Remap particle in the z direction if necessary
	int k=step_int(z*zsp);
	if(k<0||k>=nz) {
		int ak=step_div(k,nz);
		z-=ak*bz;y-=ak*byz;x-=ak*bxz;k-=ak*nz;
	}

	// Remap particle in the y direction if necessary
	int j=step_int(y*ysp);
	if(j<0||j>=ny) {
		int aj=step_div(j,ny);
		y-=aj*by;x-=aj*bxy;j-=aj*ny;
	}

	// Remap particle in the x direction if necessary
	ijk=step_int(x*xsp);
	if(ijk<0||ijk>=nx) {
		int ai=step_div(ijk,nx);
		x-=ai*bx;ijk-=ai*nx;
	}

	// Compute the block index and check memory allocation
	j+=ey;k+=ez;
	ijk+=nx*(j+oy*k);
	if(co[ijk]==mem[ijk]) add_particle_memory(ijk);
}

/** Takes a particle position vector and computes the region index into which
 * it should be stored. If the container is periodic, then the routine also
 * maps the particle position to ensure it is in the primary domain. If the
 * container is not periodic, the routine bails out.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
 *                        in, with (0,0,0) corresponding to the primary domain.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
void container_periodic_base::put_locate_block(int &ijk,double &x,double &y,double &z,int &ai,int &aj,int &ak) {

	// Remap particle in the z direction if necessary
	int k=step_int(z*zsp);
	if(k<0||k>=nz) {
		ak=step_div(k,nz);
		z-=ak*bz;y-=ak*byz;x-=ak*bxz;k-=ak*nz;
	} else ak=0;

	// Remap particle in the y direction if necessary
	int j=step_int(y*ysp);
	if(j<0||j>=ny) {
		aj=step_div(j,ny);
		y-=aj*by;x-=aj*bxy;j-=aj*ny;
	} else aj=0;

	// Remap particle in the x direction if necessary
	ijk=step_int(x*xsp);
	if(ijk<0||ijk>=nx) {
		ai=step_div(ijk,nx);
		x-=ai*bx;ijk-=ai*nx;
	} else ai=0;

	// Compute the block index and check memory allocation
	j+=ey;k+=ez;
	ijk+=nx*(j+oy*k);
	if(co[ijk]==mem[ijk]) add_particle_memory(ijk);
}

/** Takes a position vector and remaps it into the primary domain.
 * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
 *                        with (0,0,0) corresponding to the primary domain.
 * \param[out] (ci,cj,ck) the index of the block that the position vector is
 *                        within, once it has been remapped.
 * \param[in,out] (x,y,z) the position vector to consider, which is remapped
 *                        into the primary domain during the routine.
 * \param[out] ijk the block index that the vector is within. */
inline void container_periodic_base::remap(int &ai,int &aj,int &ak,int &ci,int &cj,int &ck,double &x,double &y,double &z,int &ijk) {

	// Remap particle in the z direction if necessary
	ck=step_int(z*zsp);
	if(ck<0||ck>=nz) {
		ak=step_div(ck,nz);
		z-=ak*bz;y-=ak*byz;x-=ak*bxz;ck-=ak*nz;
	} else ak=0;

	// Remap particle in the y direction if necessary
	cj=step_int(y*ysp);
	if(cj<0||cj>=ny) {
		aj=step_div(cj,ny);
		y-=aj*by;x-=aj*bxy;cj-=aj*ny;
	} else aj=0;

	// Remap particle in the x direction if necessary
	ci=step_int(x*xsp);
	if(ci<0||ci>=nx) {
		ai=step_div(ci,nx);
		x-=ai*bx;ci-=ai*nx;
	} else ai=0;

	cj+=ey;ck+=ez;
	ijk=ci+nx*(cj+oy*ck);
}

/** Takes a vector and finds the particle whose Voronoi cell contains that
 * vector. This is equivalent to finding the particle which is nearest to the
 * vector.
 * \param[in] (x,y,z) the vector to test.
 * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
 *                        contains the vector. This may point to a particle in
 *                        a periodic image of the primary domain.
 * \param[out] pid the ID of the particle.
 * \return True if a particle was found. If the container has no particles,
 * then the search will not find a Voronoi cell and false is returned. */
bool container_periodic::find_voronoi_cell(double x,double y,double z,double &rx,double &ry,double &rz,int &pid) {
	int ai,aj,ak,ci,cj,ck,ijk;
	particle_record w;
	double mrs;

	// Remap the vector into the primary domain and then search for the
	// Voronoi cell that it is within
	remap(ai,aj,ak,ci,cj,ck,x,y,z,ijk);
	vc.find_voronoi_cell(x,y,z,ci,cj,ck,ijk,w,mrs);

	if(w.ijk!=-1) {

		// Assemble the position vector of the particle to be returned,
		// applying a periodic remapping if necessary
		ci+=w.di;if(ci<0||ci>=nx) ai+=step_div(ci,nx);
		rx=p[w.ijk][3*w.l]+ak*bxz+aj*bxy+ai*bx;
		ry=p[w.ijk][3*w.l+1]+ak*byz+aj*by;
		rz=p[w.ijk][3*w.l+2]+ak*bz;
		pid=id[w.ijk][w.l];
		return true;
	}
	return false;
}

/** Takes a vector and finds the particle whose Voronoi cell contains that
 * vector. Additional wall classes are not considered by this routine.
 * \param[in] (x,y,z) the vector to test.
 * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
 *                        contains the vector. If the container is periodic,
 *                        this may point to a particle in a periodic image of
 *                        the primary domain.
 * \param[out] pid the ID of the particle.
 * \return True if a particle was found. If the container has no particles,
 * then the search will not find a Voronoi cell and false is returned. */
bool container_periodic_poly::find_voronoi_cell(double x,double y,double z,double &rx,double &ry,double &rz,int &pid) {
	int ai,aj,ak,ci,cj,ck,ijk;
	particle_record w;
	double mrs;

	// Remap the vector into the primary domain and then search for the
	// Voronoi cell that it is within
	remap(ai,aj,ak,ci,cj,ck,x,y,z,ijk);
	vc.find_voronoi_cell(x,y,z,ci,cj,ck,ijk,w,mrs);

	if(w.ijk!=-1) {

		// Assemble the position vector of the particle to be returned,
		// applying a periodic remapping if necessary
		ci+=w.di;if(ci<0||ci>=nx) ai+=step_div(ci,nx);
		rx=p[w.ijk][4*w.l]+ak*bxz+aj*bxy+ai*bx;
		ry=p[w.ijk][4*w.l+1]+ak*byz+aj*by;
		rz=p[w.ijk][4*w.l+2]+ak*bz;
		pid=id[w.ijk][w.l];
		return true;
	}
	return false;
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_periodic_base::add_particle_memory(int i) {

	// Handle the case when no memory has been allocated for this block
	if(mem[i]==0) {
		mem[i]=init_mem;
		id[i]=new int[init_mem];
		p[i]=new double[ps*init_mem];
		return;
	}

	// Otherwise, double the memory allocation for this block. Carry out a
	// check on the memory allocation size, and print a status message if
	// requested.
	int l,nmem(mem[i]<<1);
	if(nmem>max_particle_memory)
		voro_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=3
	fprintf(stderr,"Particle memory in region %d scaled up to %d\n",i,nmem);
#endif

	// Allocate new memory and copy in the contents of the old arrays
	int *idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	double *pp=new double[ps*nmem];
	for(l=0;l<ps*co[i];l++) pp[l]=p[i][l];

	// Update pointers and delete old arrays
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}

/** Import a list of particles from an open file stream into the container.
 * Entries of four numbers (Particle ID, x position, y position, z position)
 * are searched for. If the file cannot be successfully read, then the routine
 * causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_periodic::import(FILE *fp) {
	int i,j;
	double x,y,z;
	while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(i,x,y,z);
	if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position) are searched for. If the file cannot be
 * successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_periodic::import(particle_order &vo,FILE *fp) {
	int i,j;
	double x,y,z;
	while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(vo,i,x,y,z);
	if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream into the container.
 * Entries of five numbers (Particle ID, x position, y position, z position,
 * radius) are searched for. If the file cannot be successfully read, then the
 * routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_periodic_poly::import(FILE *fp) {
	int i,j;
	double x,y,z,r;
	while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(i,x,y,z,r);
	if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position, radius) are searched for. If the file
 * cannot be successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_periodic_poly::import(particle_order &vo,FILE *fp) {
	int i,j;
	double x,y,z,r;
	while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(vo,i,x,y,z,r);
	if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Outputs the a list of all the container regions along with the number of
 * particles stored within each. */
void container_periodic_base::region_count() {
	int i,j,k,*cop=co;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
		printf("Region (%d,%d,%d): %d particles\n",i,j,k,*(cop++));
}

/** Clears a container of particles. */
void container_periodic::clear() {
	for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
}

/** Clears a container of particles, also clearing resetting the maximum radius
 * to zero. */
void container_periodic_poly::clear() {
	for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
	max_radius=0;
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_periodic::print_custom(const char *format,FILE *fp) {
	c_loop_all_periodic vl(*this);
	print_custom(vl,format,fp);
}

/** Computes all the Voronoi cells and saves customized
 * information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_periodic_poly::print_custom(const char *format,FILE *fp) {
	c_loop_all_periodic vl(*this);
	print_custom(vl,format,fp);
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] filename the name of the file to write to. */
void container_periodic::print_custom(const char *format,const char *filename) {
	FILE *fp=safe_fopen(filename,"w");
	print_custom(format,fp);
	fclose(fp);
}

/** Computes all the Voronoi cells and saves customized
 * information about them
 * \param[in] format the custom output string to use.
 * \param[in] filename the name of the file to write to. */
void container_periodic_poly::print_custom(const char *format,const char *filename) {
	FILE *fp=safe_fopen(filename,"w");
	print_custom(format,fp);
	fclose(fp);
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
void container_periodic::compute_all_cells() {
	voronoicell c;
	c_loop_all_periodic vl(*this);
	if(vl.start()) do compute_cell(c,vl);
	while(vl.inc());
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
void container_periodic_poly::compute_all_cells() {
	voronoicell c;
	c_loop_all_periodic vl(*this);
	if(vl.start()) do compute_cell(c,vl);while(vl.inc());
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_periodic::sum_cell_volumes() {
	voronoicell c;
	double vol=0;
	c_loop_all_periodic vl(*this);
	if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
	return vol;
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_periodic_poly::sum_cell_volumes() {
	voronoicell c;
	double vol=0;
	c_loop_all_periodic vl(*this);
	if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
	return vol;
}

/** This routine creates all periodic images of the particles. It is meant for
 * diagnostic purposes only, since usually periodic images are dynamically
 * created in when they are referenced. */
void container_periodic_base::create_all_images() {
	int i,j,k;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) create_periodic_image(i,j,k);
}

/** Checks that the particles within each block lie within that block's bounds.
 * This is useful for diagnosing problems with periodic image computation. */
void container_periodic_base::check_compartmentalized() {
	int c,l,i,j,k;
	double mix,miy,miz,max,may,maz,*pp;
	for(k=l=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++,l++) if(mem[l]>0) {

		// Compute the block's bounds, adding in a small tolerance
		mix=i*boxx-tolerance;max=mix+boxx+tolerance;
		miy=(j-ey)*boxy-tolerance;may=miy+boxy+tolerance;
		miz=(k-ez)*boxz-tolerance;maz=miz+boxz+tolerance;

		// Print entries for any particles that lie outside the block's
		// bounds
		for(pp=p[l],c=0;c<co[l];c++,pp+=ps) if(*pp<mix||*pp>max||pp[1]<miy||pp[1]>may||pp[2]<miz||pp[2]>maz)
			printf("%d %d %d %d %f %f %f %f %f %f %f %f %f\n",
			       id[l][c],i,j,k,*pp,pp[1],pp[2],mix,max,miy,may,miz,maz);
	}
}

/** Creates particles within an image block that is aligned with the primary
 * domain in the z axis. In this case, the image block may be comprised of
 * particles from two primary blocks. The routine considers these two primary
 * blocks, and adds the needed particles to the image. The remaining particles
 * from the primary blocks are also filled into the neighboring images.
 * \param[in] (di,dj,dk) the index of the block to consider. The z index must
 *			 satisfy ez<=dk<wz. */
void container_periodic_base::create_side_image(int di,int dj,int dk) {
	int l,dijk=di+nx*(dj+oy*dk),odijk,ima=step_div(dj-ey,ny);
	int qua=di+step_int(-ima*bxy*xsp),quadiv=step_div(qua,nx);
	int fi=qua-quadiv*nx,fijk=fi+nx*(dj-ima*ny+oy*dk);
	double dis=ima*bxy+quadiv*bx,switchx=di*boxx-ima*bxy-quadiv*bx,adis;

	// Left image computation
	if((img[dijk]&1)==0) {
		if(di>0) {
			odijk=dijk-1;adis=dis;
		} else {
			odijk=dijk+nx-1;adis=dis+bx;
		}
		img[odijk]|=2;
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l]>switchx) put_image(dijk,fijk,l,dis,by*ima,0);
			else put_image(odijk,fijk,l,adis,by*ima,0);
		}
	}

	// Right image computation
	if((img[dijk]&2)==0) {
		if(fi==nx-1) {
			fijk+=1-nx;switchx+=(1-nx)*boxx;dis+=bx;
		} else {
			fijk++;switchx+=boxx;
		}
		if(di==nx-1) {
			odijk=dijk-nx+1;adis=dis-bx;
		} else {
			odijk=dijk+1;adis=dis;
		}
		img[odijk]|=1;
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l]<switchx) put_image(dijk,fijk,l,dis,by*ima,0);
			else put_image(odijk,fijk,l,adis,by*ima,0);
		}
	}

	// All contributions to the block now added, so set both two bits of
	// the image information
	img[dijk]=3;
}

/** Creates particles within an image block that is not aligned with the
 * primary domain in the z axis. In this case, the image block may be comprised
 * of particles from four primary blocks. The routine considers these four
 * primary blocks, and adds the needed particles to the image. The remaining
 * particles from the primary blocks are also filled into the neighboring
 * images.
 * \param[in] (di,dj,dk) the index of the block to consider. The z index must
 *			 satisfy dk<ez or dk>=wz. */
void container_periodic_base::create_vertical_image(int di,int dj,int dk) {
	int l,dijk=di+nx*(dj+oy*dk),dijkl,dijkr,ima=step_div(dk-ez,nz);
	int qj=dj+step_int(-ima*byz*ysp),qjdiv=step_div(qj-ey,ny);
	int qi=di+step_int((-ima*bxz-qjdiv*bxy)*xsp),qidiv=step_div(qi,nx);
	int fi=qi-qidiv*nx,fj=qj-qjdiv*ny,fijk=fi+nx*(fj+oy*(dk-ima*nz)),fijk2;
	double disy=ima*byz+qjdiv*by,switchy=(dj-ey)*boxy-ima*byz-qjdiv*by;
	double disx=ima*bxz+qjdiv*bxy+qidiv*bx,switchx=di*boxx-ima*bxz-qjdiv*bxy-qidiv*bx;
	double switchx2,disxl,disxr,disx2,disxr2;

	if(di==0) {dijkl=dijk+nx-1;disxl=disx+bx;}
	else {dijkl=dijk-1;disxl=disx;}

	if(di==nx-1) {dijkr=dijk-nx+1;disxr=disx-bx;}
	else {dijkr=dijk+1;disxr=disx;}

	// Down-left image computation
	bool y_exist=dj!=0;
	if((img[dijk]&1)==0) {
		img[dijkl]|=2;
		if(y_exist) {
			img[dijkl-nx]|=8;
			img[dijk-nx]|=4;
		}
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l+1]>switchy) {
				if(p[fijk][ps*l]>switchx) put_image(dijk,fijk,l,disx,disy,bz*ima);
				else put_image(dijkl,fijk,l,disxl,disy,bz*ima);
			} else {
				if(!y_exist) continue;
				if(p[fijk][ps*l]>switchx) put_image(dijk-nx,fijk,l,disx,disy,bz*ima);
				else put_image(dijkl-nx,fijk,l,disxl,disy,bz*ima);
			}
		}
	}

	// Down-right image computation
	if((img[dijk]&2)==0) {
		if(fi==nx-1) {
			fijk2=fijk+1-nx;switchx2=switchx+(1-nx)*boxx;disx2=disx+bx;disxr2=disxr+bx;
		} else {
			fijk2=fijk+1;switchx2=switchx+boxx;disx2=disx;disxr2=disxr;
		}
		img[dijkr]|=1;
		if(y_exist) {
			img[dijkr-nx]|=4;
			img[dijk-nx]|=8;
		}
		for(l=0;l<co[fijk2];l++) {
			if(p[fijk2][ps*l+1]>switchy) {
				if(p[fijk2][ps*l]>switchx2) put_image(dijkr,fijk2,l,disxr2,disy,bz*ima);
				else put_image(dijk,fijk2,l,disx2,disy,bz*ima);
			} else {
				if(!y_exist) continue;
				if(p[fijk2][ps*l]>switchx2) put_image(dijkr-nx,fijk2,l,disxr2,disy,bz*ima);
				else put_image(dijk-nx,fijk2,l,disx2,disy,bz*ima);
			}
		}
	}

	// Recomputation of some intermediate quantities for boundary cases
	if(fj==wy-1) {
		fijk+=nx*(1-ny)-fi;
		switchy+=(1-ny)*boxy;
		disy+=by;
		qi=di+step_int(-(ima*bxz+(qjdiv+1)*bxy)*xsp);
		int dqidiv=step_div(qi,nx)-qidiv;qidiv+=dqidiv;
		fi=qi-qidiv*nx;
		fijk+=fi;
		disx+=bxy+bx*dqidiv;
		disxl+=bxy+bx*dqidiv;
		disxr+=bxy+bx*dqidiv;
		switchx-=bxy+bx*dqidiv;
	} else {
		fijk+=nx;switchy+=boxy;
	}

	// Up-left image computation
	y_exist=dj!=oy-1;
	if((img[dijk]&4)==0) {
		img[dijkl]|=8;
		if(y_exist) {
			img[dijkl+nx]|=2;
			img[dijk+nx]|=1;
		}
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l+1]>switchy) {
				if(!y_exist) continue;
				if(p[fijk][ps*l]>switchx) put_image(dijk+nx,fijk,l,disx,disy,bz*ima);
				else put_image(dijkl+nx,fijk,l,disxl,disy,bz*ima);
			} else {
				if(p[fijk][ps*l]>switchx) put_image(dijk,fijk,l,disx,disy,bz*ima);
				else put_image(dijkl,fijk,l,disxl,disy,bz*ima);
			}
		}
	}

	// Up-right image computation
	if((img[dijk]&8)==0) {
		if(fi==nx-1) {
			fijk2=fijk+1-nx;switchx2=switchx+(1-nx)*boxx;disx2=disx+bx;disxr2=disxr+bx;
		} else {
			fijk2=fijk+1;switchx2=switchx+boxx;disx2=disx;disxr2=disxr;
		}
		img[dijkr]|=4;
		if(y_exist) {
			img[dijkr+nx]|=1;
			img[dijk+nx]|=2;
		}
		for(l=0;l<co[fijk2];l++) {
			if(p[fijk2][ps*l+1]>switchy) {
				if(!y_exist) continue;
				if(p[fijk2][ps*l]>switchx2) put_image(dijkr+nx,fijk2,l,disxr2,disy,bz*ima);
				else put_image(dijk+nx,fijk2,l,disx2,disy,bz*ima);
			} else {
				if(p[fijk2][ps*l]>switchx2) put_image(dijkr,fijk2,l,disxr2,disy,bz*ima);
				else put_image(dijk,fijk2,l,disx2,disy,bz*ima);
			}
		}
	}

	// All contributions to the block now added, so set all four bits of
	// the image information
	img[dijk]=15;
}

/** Copies a particle position from the primary domain into an image block.
 * \param[in] reg the block index within the primary domain that the particle
 *                is within.
 * \param[in] fijk the index of the image block.
 * \param[in] l the index of the particle entry within the primary block.
 * \param[in] (dx,dy,dz) the displacement vector to add to the particle. */
void container_periodic_base::put_image(int reg,int fijk,int l,double dx,double dy,double dz) {
	if(co[reg]==mem[reg]) add_particle_memory(reg);
	double *p1=p[reg]+ps*co[reg],*p2=p[fijk]+ps*l;
	*(p1++)=*(p2++)+dx;
	*(p1++)=*(p2++)+dy;
	*p1=*p2+dz;
	if(ps==4) *(++p1)=*(++p2);
	id[reg][co[reg]++]=id[fijk][l];
}

}
