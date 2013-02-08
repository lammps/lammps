#!/bin/sed -f
# sed edit script to convert files in USER-OMP for better compiler optimization

# replace shortcuts to per atom arrays
s/const double \* const \* const x = atom->x;/const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];/
s/const double \* const \* const v = atom->v;/const dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];/
s/const double \* const \* const omega = atom->omega;/const dbl3_t * _noalias const omega = (dbl3_t *) atom->omega[0];/
s/const double \* const \* const mu = atom->mu;/const dbl4_t * _noalias const mu = (dbl4_t *) atom->mu[0];/
s/double \* const \* const f = thr->get_f();/dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];/
s/double \* const \* const tor = thr->get_torque();/dbl3_t * _noalias const tor = (dbl3_t *) thr->get_torque()[0];/
s/const double \* const q = atom->q;/const double * _noalias const q = atom->q;/
s/const double \* const radius = atom->radius;/const double * _noalias const radius = atom->radius;/
s/const double \* const rmass = atom->rmass;/const double * _noalias const rmass = atom->rmass;/
s/const int \* const type = atom->type;/const int * _noalias const type = atom->type;/
s/const int \* const tag = atom->tag;/const int * _noalias const tag = atom->tag;/
s/const double \* const special_lj/const double * _noalias const special_lj/
s/const double \* const special_coul/const double * _noalias const special_coul/
s/const int \* const ellipsoid = atom->ellipsoid;/const int * _noalias const ellipsoid = atom->ellipsoid;/
s/AtomVecEllipsoid::Bonus *bonus = avec->bonus;/const AtomVecEllipsoid::Bonus * _noalias const bonus = avec->bonus;/
s/const int \* const \* const bondlist = neighbor->bondlist;/const int3_t * _noalias const bondlist = (int3_t *) neighbor->bondlist[0];/
s/const int \* const \* const anglelist = neighbor->anglelist;/const int4_t * _noalias const anglelist = (int4_t *) neighbor->anglelist[0];/
s/const int \* const \* const dihedrallist = neighbor->dihedrallist;/const int5_t * _noalias const dihedrallist = (int5_t *) neighbor->dihedrallist[0];/
s/const int \* const \* const improperlist = neighbor->improperlist;/const int5_t * _noalias const improperlist = (int5_t *) neighbor->improperlist[0];/


# change refereces to coordinate elements to use dbl3_t members instead
s/x\[\(.\)\]\[0\]/x[\1].x/g
s/x\[\(.\)\]\[1\]/x[\1].y/g
s/x\[\(.\)\]\[2\]/x[\1].z/g
s/v\[\(.\)\]\[0\]/v[\1].x/g
s/v\[\(.\)\]\[1\]/v[\1].y/g
s/v\[\(.\)\]\[2\]/v[\1].z/g
s/f\[\(.\)\]\[0\]/f[\1].x/g
s/f\[\(.\)\]\[1\]/f[\1].y/g
s/f\[\(.\)\]\[2\]/f[\1].z/g
s/tor\[\(.\)\]\[0\]/tor[\1].x/g
s/tor\[\(.\)\]\[1\]/tor[\1].y/g
s/tor\[\(.\)\]\[2\]/tor[\1].z/g
s/omega\[\(.\)\]\[0\]/omega[\1].x/g
s/omega\[\(.\)\]\[1\]/omega[\1].y/g
s/omega\[\(.\)\]\[2\]/omega[\1].z/g
s/mu\[\(.\)\]\[0\]/mu[\1].x/g
s/mu\[\(.\)\]\[1\]/mu[\1].y/g
s/mu\[\(.\)\]\[2\]/mu[\1].z/g
s/mu\[\(.\)\]\[3\]/mu[\1].w/g
s/x\[\(..\)\]\[0\]/x[\1].x/g
s/x\[\(..\)\]\[1\]/x[\1].y/g
s/x\[\(..\)\]\[2\]/x[\1].z/g
s/v\[\(..\)\]\[0\]/v[\1].x/g
s/v\[\(..\)\]\[1\]/v[\1].y/g
s/v\[\(..\)\]\[2\]/v[\1].z/g
s/f\[\(..\)\]\[0\]/f[\1].x/g
s/f\[\(..\)\]\[1\]/f[\1].y/g
s/f\[\(..\)\]\[2\]/f[\1].z/g
s/tor\[\(..\)\]\[0\]/tor[\1].x/g
s/tor\[\(..\)\]\[1\]/tor[\1].y/g
s/tor\[\(..\)\]\[2\]/tor[\1].z/g
s/omega\[\(..\)\]\[0\]/omega[\1].x/g
s/omega\[\(..\)\]\[1\]/omega[\1].y/g
s/omega\[\(..\)\]\[2\]/omega[\1].z/g
s/mu\[\(..\)\]\[0\]/mu[\1].x/g
s/mu\[\(..\)\]\[1\]/mu[\1].y/g
s/mu\[\(..\)\]\[2\]/mu[\1].z/g
s/mu\[\(..\)\]\[3\]/mu[\1].w/g
s/bondlist\[\(.\)\]\[0\]/bondlist[\1].a/g
s/bondlist\[\(.\)\]\[1\]/bondlist[\1].b/g
s/bondlist\[\(.\)\]\[2\]/bondlist[\1].t/g
s/anglelist\[\(.\)\]\[0\]/anglelist[\1].a/g
s/anglelist\[\(.\)\]\[1\]/anglelist[\1].b/g
s/anglelist\[\(.\)\]\[2\]/anglelist[\1].c/g
s/anglelist\[\(.\)\]\[3\]/anglelist[\1].t/g
s/dihedrallist\[\(.\)\]\[0\]/dihedrallist[\1].a/g
s/dihedrallist\[\(.\)\]\[1\]/dihedrallist[\1].b/g
s/dihedrallist\[\(.\)\]\[2\]/dihedrallist[\1].c/g
s/dihedrallist\[\(.\)\]\[3\]/dihedrallist[\1].d/g
s/dihedrallist\[\(.\)\]\[4\]/dihedrallist[\1].t/g
s/improperlist\[\(.\)\]\[0\]/improperlist[\1].a/g
s/improperlist\[\(.\)\]\[1\]/improperlist[\1].b/g
s/improperlist\[\(.\)\]\[2\]/improperlist[\1].c/g
s/improperlist\[\(.\)\]\[3\]/improperlist[\1].d/g
s/improperlist\[\(.\)\]\[4\]/improperlist[\1].t/g
