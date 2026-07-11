/*
 * convert a LAMMPS / Smooth-Mach Dynamics triangle output file into a VTK unstructured grid file.
 * LAMMPS output file contains STL triangle vertices
 *
 * ASSUMPTIONS:
 *
 * A compute exists for and STL triangulated surface, which computes the triangl vertices with a command like this:
 * compute  F tri_group smd/triangle_vertices
 * c_F[1-9] now hold the triangle vertices.
 *
 * LAMMPS dump file has the following entries per atom line:
 * id type mol x y z c_F[1] c_F[2] c_F[3] c_F[4] c_F[5] c_F[6] c_F[7] c_F[8] c_F[9]
 *
 * Author: Georg Ganzenmuller at the 
 * Fraunhofer-Institute for High-Speed Dynamics,
 * Ernst Mach Institute in Germany
 * (georg.ganzenmueller at emi.fhg.de)
 */

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

#define RENDER 0 // switch on / off rendering of snapshots

void WriteVTU(std::string filename, int ntime, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	// create filename with the current timestep
					std::ostringstream oss;
					oss << ntime;
					filename += "." + oss.str() + ".vtu";

					writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
					writer->SetInput(unstructuredGrid);
	#else
					writer->SetInputData(unstructuredGrid);
	#endif
					writer->Write();
}

void WriteVTK(std::string filename, int ntime, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	// create filename with the current timestep
					std::ostringstream oss;
					oss << ntime;
					filename += "_" + oss.str() + ".vtk";

					writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
					writer->SetInput(unstructuredGrid);
	#else
					writer->SetInputData(unstructuredGrid);
	#endif
					writer->Write();
}

int main(int argc, char *argv[]) {

	string line;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkCellArray> vertices;
	vtkSmartPointer<vtkCellArray> cellArray;
	vtkIdType id_v1, id_v2, id_v3;
	vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New(); // create an empty tri;
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;


	unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Parse command line arguments
	if (argc != 3) {
		std::cout << "Required arguments: InputFilename OutputFile_Basename"
				<< std::endl;
		return EXIT_FAILURE;
	}

	std::string lammps_dump = argv[1];
	std::string filename = argv[2];

	ifstream myfile(lammps_dump.c_str());

	int nlines_read = 0;
	int count = -1;
	int ntime;
	int N; // number of atoms in this timestep
	int id, type, mol;
	double x, y, z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;
	if (myfile.is_open()) {

		while (getline(myfile, line)) {
			std::size_t found = line.find("TIMESTEP");
			if (found != std::string::npos) {
				cout << line << endl;
				count = 0;
				nlines_read = 0;

				// reinitialize VTK arrays so they are empty
				points = vtkSmartPointer<vtkPoints>::New();
				vertices = vtkSmartPointer<vtkCellArray>::New();
				cellArray = vtkSmartPointer<vtkCellArray>::New();
			}
			if (count >= 0)
				count++;

			if (count == 2) {
				istringstream convert(line);
				convert >> ntime;
				cout << "timestep is: " << ntime << endl;
			}

			if (count == 4) {
				istringstream convert(line);
				convert >> N;
				cout << "number of particles is: " << N << endl;
			}

			if (count > 9) {
				//cout << line << endl;
				//const vector<string> words = split(line, " ");

				std::istringstream is(line);
				std::istream_iterator<float> eos;
				std::vector<float> out(std::istream_iterator<float>(is), eos);

				/*
				 * we expect the following entries in a data line:
				 * id type mol x y z v1x v1y v1z v2x v2y v2z v3x v3y v3z
				 */

				id = out[0];
				type = out[1];
				mol = out[2];
				x = out[3];  // triangle center
				y = out[4];
				z = out[5];
				v1x = out[6]; // triangle vertex 1
				v1y = out[7];
				v1z = out[8];
				v2x = out[9]; // triangle vertex 2
				v2y = out[10];
				v2z = out[11];
				v3x = out[12]; // triangle vertex 3
				v3y = out[13];
				v3z = out[14];

				id_v1 = points->InsertNextPoint(v1x, v1y, v1z); // store tri vertices in points array and retain indices
				id_v2 = points->InsertNextPoint(v2x, v2y, v2z);
				id_v3 = points->InsertNextPoint(v3x, v3y, v3z);



				tri->GetPointIds()->SetId(0, id_v1); // fill triangle with retained indices
				tri->GetPointIds()->SetId(1, id_v2);
				tri->GetPointIds()->SetId(2, id_v3);

				cellArray->InsertNextCell(tri); // insert tri into array of cells

				nlines_read += 1;

			}

			if (count % (9 + N) == 0) { // we have reached the end of this snapshot
				printf(
						"----------------------------------------------------\n");
				printf("number of lines read from lammps dump file: %d\n",
						nlines_read);


				unstructuredGrid->SetPoints(points);
				unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);

				WriteVTK(filename, ntime, unstructuredGrid);


//				if (ntime > 100) {
//					return 0;
//				}

#if RENDER == 1
				// Read and display file for verification that it was written correclty
				vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
						vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
				reader->SetFileName(filename.c_str());
				reader->Update();

				vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<
						vtkDataSetMapper>::New();
				mapper->SetInputConnection(reader->GetOutputPort());

				vtkSmartPointer<vtkActor> actor =
						vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);

				vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<
						vtkRenderer>::New();
				vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<
						vtkRenderWindow>::New();
				renderWindow->AddRenderer(renderer);
				vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
						vtkSmartPointer<vtkRenderWindowInteractor>::New();
				renderWindowInteractor->SetRenderWindow(renderWindow);

				renderer->AddActor(actor);
				renderer->SetBackground(.3, .6, .3); // Background color green

				renderWindow->Render();
				renderWindowInteractor->Start();
#endif
			}

		}
		myfile.close();
	} else {
		cout << "Unable to open file";
	}

	return EXIT_SUCCESS;
}
