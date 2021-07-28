/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <array>
#include <regex>
#include <string.h>
#include <cmath>
//#include <filesystem>

static const std::string data_file0 = "system.init";
static const std::string data_dump0 = "config.dump";
static const std::string out_dir0 = "out";

struct Particle {
	double x, y, z, vx, vy, vz, Es, Eb, Et, Ep, Ek;
	char type, nx, ny, nz;
};

class Lamps_base {
public:
	Lamps_base() = default;
	virtual ~Lamps_base() = default;

	int open(const std::string& filename);
	int next();											//get next snapshot from the opened file	
	virtual int write(const std::string& filename) const = 0;

	inline double get_X1() const { return X1; };
	inline double get_X2() const { return X2; };
	inline double get_Y1() const { return Y1; };
	inline double get_Y2() const { return Y2; };
	inline double get_Z1() const { return Z1; };
	inline double get_Z2() const { return Z2; };
	inline int get_Natoms() const { return Natoms; };
	inline int get_Nsteps() const { return Nsteps; };
	inline int is_open() const { return open_stat; };
	inline const Particle& get(int i) const { return particles[i]; };
	inline Particle& get(int i) { return particles[i]; };

protected:
	virtual int load() = 0;

	int Nsteps, Natoms, open_stat;
	double X1, X2, Y1, Y2, Z1, Z2;
	std::vector<Particle> particles;
	std::ifstream in;
};

class Lamps_dump : public Lamps_base {
public:
	Lamps_dump() = default;
	~Lamps_dump() = default;

	virtual int write(const std::string& filename) const override;
private:
	virtual int load() override;
};

int Lamps_base::open(const std::string& filename) {
	in.open(filename); if (!in.is_open()) return EXIT_FAILURE;
	return load();
}

int Lamps_base::next() {
	return load();
}

int Lamps_dump::write(const std::string& filename) const {
	return EXIT_FAILURE;
}

int Lamps_dump::load() {
	std::string inbuf; char* tmp_cptr;
	open_stat = 0;

	if (!getline(in, inbuf)) return EXIT_FAILURE;
	if (!getline(in, inbuf)) return EXIT_FAILURE;
	Nsteps = std::stoi(inbuf);
	if (!getline(in, inbuf)) return EXIT_FAILURE;
	if (!getline(in, inbuf)) return EXIT_FAILURE;
	Natoms = std::stoi(inbuf);
	particles.resize(Natoms);
	if (!getline(in, inbuf)) return EXIT_FAILURE;

	if (!getline(in, inbuf)) return EXIT_FAILURE;
	X1 = strtof(inbuf.c_str(), &tmp_cptr);
	X2 = strtof(tmp_cptr + 1, &tmp_cptr);

	if (!getline(in, inbuf)) return EXIT_FAILURE;
	Y1 = strtof(inbuf.c_str(), &tmp_cptr);
	Y2 = strtof(tmp_cptr + 1, &tmp_cptr);

	if (!getline(in, inbuf)) return EXIT_FAILURE;
	Z1 = strtof(inbuf.c_str(), &tmp_cptr);
	Z2 = strtof(tmp_cptr + 1, &tmp_cptr);

	if (!getline(in, inbuf)) return EXIT_FAILURE;
	for (int i = 0; i < Natoms; i++) {
		if (!getline(in, inbuf)) return EXIT_FAILURE;
		int id = strtol(inbuf.c_str(), &tmp_cptr, 10) - 1;		// modify based on a particular file format
		particles[id].type = static_cast<char>(strtol(tmp_cptr + 1, &tmp_cptr, 10));
		particles[id].x = strtof(tmp_cptr + 1, &tmp_cptr);
		particles[id].y = strtof(tmp_cptr + 1, &tmp_cptr);
		particles[id].z = strtof(tmp_cptr + 1, &tmp_cptr);
		particles[id].Es = strtof(tmp_cptr + 1, &tmp_cptr);
		particles[id].Eb = strtof(tmp_cptr + 1, &tmp_cptr);
		particles[id].Et = strtof(tmp_cptr + 1, &tmp_cptr);
		particles[id].Ep = particles[id].Es + particles[id].Eb + particles[id].Et;
		particles[id].Ek = strtof(tmp_cptr + 1, &tmp_cptr);
	}
	open_stat = true;

	return EXIT_SUCCESS;
}


int main(int argc, char* argv[]) {
	std::string data_file = (argc > 1) ? argv[1] : data_file0;
	std::string data_dump = (argc > 2) ? argv[2] : data_dump0;
	std::string out_dir  = (argc > 3) ? argv[3] : out_dir0;
	//std::filesystem::remove_all(out_dir);
	//std::filesystem::create_directories(out_dir);

	//list of bonds
	std::ifstream in(data_file); 
	if (!in.is_open()) {
		std::cout << "cannot open " << data_file << std::endl;
		return EXIT_FAILURE;
	}
	std::string buf;
	std::string atoms_l = "Atoms";
	while (std::getline(in, buf)){
		if (buf == atoms_l) break;
		if (in.eof()) return EXIT_FAILURE;
	}
	std::getline(in, buf);
	char* tmp_cptr;
	std::vector<std::array<int, 2>> bonds;
	while (std::getline(in, buf)) {
		if (in.eof() || buf.size() == 0) break;
		int idx = strtol(buf.c_str(), &tmp_cptr, 10);
		int m_idx = strtol(tmp_cptr + 1, &tmp_cptr, 10);
		int type = strtol(tmp_cptr + 1, &tmp_cptr, 10);
		int id1 = strtol(tmp_cptr + 1, &tmp_cptr, 10);
		int id2 = strtol(tmp_cptr + 1, &tmp_cptr, 10);
		if(id1 >= 0 && id2 >= 0) bonds.push_back({id1 - 1, id2 - 1});
	}

	//dump
	Lamps_dump dump;
	dump.open(data_dump);
	if (!dump.is_open()) {
		std::cout << "cannot open " << data_dump << std::endl;
		return EXIT_FAILURE;
	}
	double Lx = dump.get_X2() - dump.get_X1();
	double Ly = dump.get_Y2() - dump.get_Y1();
	double Lz = dump.get_Z2() - dump.get_Z1();
	while (1) {
		std::ofstream out(out_dir + "/cnt" + std::to_string(dump.get_Nsteps()) + ".vtk");
		if (!out.is_open()) {
			std::cout << "cannot create " << out_dir + "/cnt" + std::to_string(dump.get_Nsteps()) + ".vtk" << std::endl;
			std::cout << "create the output directory \"" << out_dir << "\" manually" << std::endl;
			return EXIT_FAILURE;
		}
		out << "# vtk DataFile Version 3.0\n# \nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
		out << "POINTS " << dump.get_Natoms() << " float\n";
		for (int i = 0; i < dump.get_Natoms(); i++) {
			out << dump.get(i).x << " " << dump.get(i).y << " " << dump.get(i).z << " " << "\n";
		}

		int bond_count = 0;
		for (int i = 0; i < bonds.size(); i++) {
			double f1 = std::fabs(dump.get(bonds[i][0]).x - dump.get(bonds[i][1]).x);
			double f2 = std::fabs(dump.get(bonds[i][0]).y - dump.get(bonds[i][1]).y);
			double f3 = std::fabs(dump.get(bonds[i][0]).z - dump.get(bonds[i][1]).z);
			if ((std::fabs(dump.get(bonds[i][0]).x - dump.get(bonds[i][1]).x) < 0.5*Lx)
			 && (std::fabs(dump.get(bonds[i][0]).y - dump.get(bonds[i][1]).y) < 0.5*Ly)
			 && (std::fabs(dump.get(bonds[i][0]).z - dump.get(bonds[i][1]).z) < 0.5*Lz))
				bond_count++;
		}
		out << "\nCELLS " << bond_count << " " << 3*bond_count << "\n";
		for (int i = 0; i < bonds.size(); i++) {
			if ((std::fabs(dump.get(bonds[i][0]).x - dump.get(bonds[i][1]).x) < 0.5 * Lx)
				&& (std::fabs(dump.get(bonds[i][0]).y - dump.get(bonds[i][1]).y) < 0.5 * Ly)
				&& (std::fabs(dump.get(bonds[i][0]).z - dump.get(bonds[i][1]).z) < 0.5 * Lz))
					out << "2 " << bonds[i][0] << " " << bonds[i][1] << " " << "\n";
		}

		out << "\nCELL_TYPES " << bond_count << "\n";
		for (int i = 0; i < bond_count; i++) {
			out << "4\n";
		}

		out << "\nPOINT_DATA " << dump.get_Natoms() << "\n";
		out << "SCALARS Ep float 1\n";
		out << "LOOKUP_TABLE default\n";
		for (int i = 0; i < dump.get_Natoms(); i++) {
			out << dump.get(i).Ep << "\n";
		}

		out << "\nSCALARS Ek float 1\n";
		out << "LOOKUP_TABLE default\n";
		for (int i = 0; i < dump.get_Natoms(); i++) {
			out << dump.get(i).Ek << "\n";
		}

		out << "\nSCALARS Es float 1\n";
		out << "LOOKUP_TABLE default\n";
		for (int i = 0; i < dump.get_Natoms(); i++) {
			out << dump.get(i).Es << "\n";
		}

		out << "\nSCALARS Eb float 1\n";
		out << "LOOKUP_TABLE default\n";
		for (int i = 0; i < dump.get_Natoms(); i++) {
			out << dump.get(i).Eb << "\n";
		}

		out << "\nSCALARS Et float 1\n";
		out << "LOOKUP_TABLE default\n";
		for (int i = 0; i < dump.get_Natoms(); i++) {
			out << dump.get(i).Et << "\n";
		}

		if (dump.next() != EXIT_SUCCESS) break;
	}
	return EXIT_SUCCESS;
}
