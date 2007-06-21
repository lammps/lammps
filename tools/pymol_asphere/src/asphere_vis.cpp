/***************************************************************************
                               asphere_vis.cpp
                             -------------------

  Convert a Lammps ellipsoid trajectory to a Pymol CGO trajectory

 __________________________________________________________________________
    This file is part of the Graphics Utilities package for command-line
    access to Graphics Library functions
__________________________________________________________________________

    begin                : Fri Jan 12 2007
    copyright            : (C) 2007 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#include "commandline.h"
#include "glsurface.h"
#include <limits>

// Describe the program parameters
void Describe(CommandLine &cl,ostream &out);
// Parse the command line parameters
void HandleArgs(CommandLine &cl, int argc, char *argv[], Error *error);
// Parse past an ITEM line in the lammps dump file
bool parse_to(const char *token,ifstream &in);

int main(int argc, char *argv[]) {
  CommandLine cl;
  Error error;
  Colors colors;
  FileIterator fi;

  // Parse the command line
  HandleArgs(cl,argc,argv,&error);

  // ----------------- Set up file names
  if (cl['s']) {
    fi.set_file_header(a::namewoext(cl.argstring(' ',2))+".");  
    fi.set_file_extensions("."+a::extension(cl.argstring(' ',2)));
    fi.set_lead_zeros(4);
  }
  
  // ----------------- Get the frame interval
  unsigned start_frame=0;
  unsigned skip=0;
  unsigned end_frame=std::numeric_limits<unsigned>::max();
  if (cl['i']) {
    start_frame=cl.argint('i',0);
    skip=cl.argint('i',1);
    end_frame=cl.argint('i',2);
  }
  
  // ----------------- Get the atom type info from a data file
  unsigned atom_types=0;
  vector<cPt> shapes;
  ifstream in;
  if (!cl['d']) {
    a::fileopen(in,cl.argstring(' ',0),error);
    while (!in.eof()) {
      char iline[500];
      in.getline(iline,500);
      vector<string> tokens;
      a::get_tokens(iline,tokens);
      if (tokens.size()>2)
        if (tokens[1]=="atom" && tokens[2]=="types") {
          atom_types=atoi(tokens[0].c_str());
          break;
        }
      if (!in)
        break;
    }
    if (atom_types==0)
      error.generate_error(0,"asphere_vis",
          "Could not find number of atom types in data file.");
    // ----------------- Get the atom type shapes
    if (parse_to("Shapes",in))
      for (unsigned i=0; i<atom_types; i++) {
        double tempd;
        cPt shape;
        in >> tempd >> shape;
        shape.x *= 0.5;
        shape.y *= 0.5;
        shape.z *= 0.5;
        if (in.eof() || !in)
          break;
        shapes.push_back(shape);
        char iline[500];
        in.getline(iline,500);
      }
    if (shapes.size()!=atom_types) {
      error.buffer() << "Error reading shapes from Pair Coeffs section of "
                    << "data file. Read " << shapes.size() << " valid shapes, "
                    << "but expected " << atom_types;
      error.addbuf(0,"asphere_vis");
    }
    a::fileclose(in,error);
  } else {
    // ----------------- Get the atom type info from a input file
    a::fileopen(in,cl.argstring(' ',0),error);
    while (!in.eof()) {
      string token;
      in >> token;
      if (token=="create_box") {
        in >> atom_types;
        shapes.assign(atom_types,cPt(0.5,0.5,0.5));
      } else if (token=="shape") {
        unsigned type;
        cPt shape;
        in >> type >> shape;
        shape.x *= 0.5;
        shape.y *= 0.5;
        shape.z *= 0.5;
        if (type>atom_types) {
          error.buffer() << "Error reading shapes from LAMMPS input file. "
                         << "I thought there were " << atom_types
                         << " atom types. But found an shape command for "
                         << "atom type: " << type;
          error.addbuf(0,"asphere_vis");
        }
        shapes[type-1]=shape;
      } else {
        char iline[500];
        in.getline(iline,500);
      }
      
      if (!in && !in.eof())
        error.generate_error(0,"asphere_vis",
                             "Error reading from LAMMPS input file");
    }
    a::fileclose(in,error);
  }
  if (atom_types==0)
    error.generate_error(0,"asphere_vis","Found 0 atom types!");
    
  // ----------------- Get the colors and alpha values for atom types
  vector<colorPt> color_list;
  vector<double> alpha_list;
  if (cl['c']) {
    a::fileopen(in,cl.argstring('c',0),error);
    for (unsigned i=0; i<atom_types; i++) {
      double alpha;
      string color;
      in >> color >> alpha;
      if (in.eof() || !in)
        error.generate_error(0,"asphere_vis",
                             "Improperly formatted color file.");
      color_list.push_back(colors[color]);
      alpha_list.push_back(alpha);
    }
    a::fileclose(in,error);
  }
  
  a::fileopen(in,cl.argstring(' ',1),error);
  ofstream out;
  if (!cl['s']) {
    a::fileopen(out,cl.argstring(' ',2),error);
  }
  

  // ----------------- Get the atom count
  unsigned atom_count;
 if (parse_to("ITEM: NUMBER OF ATOMS",in))
   in >> atom_count;
 else
   error.generate_error(0,"asphere_vis",
        "Could not find ITEM: NUMBER OF ATOMS in input file.");
 if (!in)
   error.generate_error(0,"asphere_vis",
        "Error reading ITEM: NUMBER OF ATOMS in input file.");

  // ----------------- Get the triangles per ellipsoid
  unsigned ellip_res=10;
  if (cl['r'])
    ellip_res=cl.argint('r',0);
  if (ellip_res==0) {
    error.addwarning(0,9,"asphere_vis","Cannot use -r 0. Setting to 10.");
    ellip_res=10;
  }
  
  // ----------------- Get the bounding box
  bool bound_found=false;
  if (!cl['s']) {
    if (parse_to("ITEM: BOX BOUNDS",in)) {
      bound_found=true;
      cPt bound[2];
     in >> bound[0].x >> bound[1].x;
     in >> bound[0].y >> bound[1].y;
     in >> bound[0].z >> bound[1].z;
      GLSurface gls;
      Vertex v;
      v.transparency=1;
      v.valid_normal=false;
      v.color=colors["white"];
      for (unsigned i=0; i<2; i++)
        for (unsigned j=0; j<2; j++)
          for (unsigned k=0; k<2; k++) {
           v.cpt=cPt(bound[i].x,bound[j].y,bound[k].z);
            gls.addvertex(v);
          }
      gls.addline(0,1);
      gls.addline(0,2);
      gls.addline(0,4);
      gls.addline(1,3);
      gls.addline(1,5);
      gls.addline(2,3);
      gls.addline(2,6);
      gls.addline(3,7);
      gls.addline(4,5);
      gls.addline(4,6);
      gls.addline(5,7);
      gls.addline(6,7);
      gls.writelines(out,"gridb");
      out << "cmd.set('cgo_dot_width',8)\n";
   } else
     error.addwarning(0,9,"asphere_vis",
          "Could not find ITEM: BOX BOUNDS in input file. No box output.");
  }
  
 if (!in)
   error.generate_error(0,"asphere_vis",
        "Error reading ITEM: BOX BOUNDS.");
  
  // ----------------- Generate the frames
  unsigned frame=0;
  unsigned max_frame=std::numeric_limits<unsigned>::max();
  if (cl['f'])
    max_frame=cl.argint('f',0);
    
  colorPt color=colors["marine"];
  double alpha=1.0;
  Vertex v;
  v.color=color;
  v.transparency=1.0;

  // ----------------- Get to the start frame
  while (frame<start_frame)
    if (!parse_to("ITEM: ATOMS",in))
      error.generate_error(0,"asphere_vis",
                      "Could not find first frame in interval in dump file.");
    else
      frame++;

  unsigned wrote=0;      
  while (true) {
    if (frame>end_frame)
       break;
    if (!parse_to("ITEM: ATOMS",in))
      break;
    GLSurface gls;
   for (unsigned i=0; i<atom_count; i++) {
      unsigned id;
      unsigned atom_type;
    in >> id >> atom_type;
    cPt atom_center;
    in >> atom_center;
      Quaternion q;
      in >> q;
    if (!in) {
      error.addwarning(0,9,"asphere_vis","Error reading frame: "+
                      a::itoa(frame));
     break;
    }
      cPt radius(shapes[atom_type-1]);
      if (radius.x == radius.y && radius.y == radius.z) {
        v.cpt=atom_center;
        if (cl['c']) {
          v.color=color_list[atom_type-1];
          v.transparency=alpha_list[atom_type-1];
        }
        gls.addvertex(v);
        gls.add_sphere(gls.size_vertices()-1,radius.x);
      } else {
        if (cl['c']) {
          color=color_list[atom_type-1];
          alpha=alpha_list[atom_type-1];
        }
        gls.add_ellipsoid(atom_center,radius,q,
                        color,alpha,ellip_res);
      }
   }
    if (!in)
      break;
    if (cl['s']) {
      if (cl['b'])
        fi.set_file_num(frame);
      a::fileopen(out,fi.nextfilename(),error);
    }
    gls.writetris(out,"ellipse");
    if (gls.size_spheres()!=0)
      gls.writespheres(out,"spheres");
    if (cl['s'])
      a::fileclose(out,error);
    wrote++;
    frame++;
    if (frame==max_frame)
      break;
    for (unsigned i=0; i<skip; i++)
      if (!parse_to("ITEM: ATOMS",in))
        break;
      else
        frame++;
  }
  if (frame==0)
    error.addwarning(0,9,"asphere_vis",
                     "Could not find any frams in input_file!");
  if (cl['i'] && frame<end_frame) {
    error.buffer() << "Only found " << frame << " frames in input file.";
    error.addbuf(0,9,"asphere_vis");
  }

  if (bound_found)
    out << "cmd.zoom(\"gridb\",animate=-1)\n";
    
  cout << "Wrote " << wrote << " frames to output file.\n";
  a::fileclose(out,error);

 return 0;
}

void Describe(CommandLine &cl,ostream &out) {
 string name=cl.program_name();
  string progname=a::strcenter(name,70);
  string gridversion=a::strcenter("Graphics Library Version 0.1",70);
 out << endl << progname << endl << gridversion << endl
   << "______________________________________________________________________\n"
   << a::strcenter("W. Michael Brown",70) << endl
   << a::strcenter("1/12/2007",70) << endl
   << "______________________________________________________________________\n"
   << "Tool for LAMMPS aspherical trajectory visualization in pymol.\n\n"
   << cl.format_synopsis("","","") << endl << endl
   << "Use '" << name << " -h > " << name
   << ".1' to generate a man page for this\n"
   << "program and type 'man ./" << name << ".1' for help\n"
   << "______________________________________________________________________\n";
 return;
}

void HandleArgs(CommandLine &cl, int argc, char *argv[], Error *error) {
  // Arguments
  cl.addmanditory(' ',3);
  cl.addargname(' ',"input_data_file");
 cl.addargname(' ',"input_dump_file");
 cl.addargname(' ',"output_py_file");
  cl.add('d',0);
  cl.adddescription('d',"Use a LAMMPS input file rather than a data file for extracting atom type information. The input filename is specified as input_data_file.");
  cl.add('f',1);
  cl.addargname('f',"max_frame");
  cl.adddescription('f',"Do not write more than max_frame frames to the output file.");
  cl.add('r',1);
  cl.addargname('r',"ellip_res");
  cl.adddescription('r',"Resolution of ellipsoids in trajectory. The number of triangles per ellipsoid is equal to 2*(ellip_res^2). Default is 10.");
  cl.add('c',1);
  cl.addargname('c',"color_file");
  cl.adddescription('c',"Color the atom_types and set transparency based on the input file. The input file contains a space delimited set sequence of the color for an atom followed by the alpha. The color should be the string name and the alpha should be between 0 and 1.");
  cl.add('s',0);
  cl.adddescription('s',"Output the results into separate .py files. The filename and extension for the output files is taken from output_py_file.");
  cl.add('i',3);
  cl.addargname('i',"start_frame");
  cl.addargname('i',"skip");
  cl.addargname('i',"end_frame");
  cl.adddescription('i',"Render the specified frame interval inclusive between start_frame and end_frame. skip gives the number of frames to skip between each rendered frame. A value of 0 outputs every frame between start_frame and end_frame. The first frame in the dump file is frame 0.");
  cl.add('b',0);
  cl.adddescription('b',"When used with -s, the option will number the filenames based on the frame number. By default, they are numbered consequtively from zero.");
        
  // Stuff for every executable
  cl.addhelp('h',0);
  cl.adddescription('h',"Print out the man page for help");
  cl.add('n',1);
  cl.addargname('n',"notice_level");
  cl.adddescription('n',"Set the degree of program output.  Use: \n\n\t-n  0\tNo output\n\t-n 10\tNormal program output\n\t-n 20\tParameters useful for reproducing the results\n\t-n 30\tAll output");

  // Short Description
  cl.addtoman_chapter("NAME","Tools for LAMMPS ellipsoid trajectory visualization in PyMol.");

  // Version
  cl.addtoman_chapter("VERSION","Version 0.1");

  // Full Description
  const string desc[5]={
  "Tool for LAMMPS trajectory visualization in PyMol. The input is a LAMMPS ",
    "'data' file or a 'in' file with ellipsoid semi-axes specified using the ",
    "ellipsoid command. The trajectory is input from a 'dump' file that must ",
    "be generated using a custom style with the following arguments in order:\n",
    ".TP\\fItag type x y z quatw quati quatj quatk\\fR\n"
 };
 cl.addtoman_chapter("DESCRIPTION",5,desc);

  Colors colors;
  cl.addtoman_chapter("AVAILABLE COLORS",colors.colorlist());
  
  // bugs
  const string bugs[3]={
  "Comments are not allowed at any point between a section header and ",
    "the end of the contents for a section in either the data file or ",
    "the input file."
 };
 cl.addtoman_chapter("KNOWN BUGS",3,bugs);


  // Authors
 cl.addtoman_chapter("AUTHORS","W. Michael Brown");

  // Parse the commandline
  if (!cl.parse(argc,argv,error)) {
  Describe(cl,cout);
  error->generate_error(0,a::filenameonly(argv[0]),"Bad Command Line\n");
 }

  // Set the notice level
  if (cl['n'])
  error->note.set_notice_level(cl.argint('n',0));

  // Generate a notice with the command line for records purposes
  string cm=cl.program_name();
  for (int j=1; j<argc; j++)
  cm+=' '+string(argv[j]);
 cm+="\n";
  error->note.notice(19,"CommandLine",cm);

  // Output the help
  if (cl['h']) {
  cl.write_man_page(cout,"0.1","Graphics Utilities");
  exit(0);
 }
}

// Parse past an ITEM line in the lammps dump file
bool parse_to(const char * token,ifstream &in) {
  char iline[5000];
  while (true) {
    in.getline(iline,5000);
    if (in.eof() || !in)
      return false;
    if (strcmp(token,iline)==0)
      return true;
  }
}
