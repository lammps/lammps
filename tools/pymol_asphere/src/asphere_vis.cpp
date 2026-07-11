/***************************************************************************
                               asphere_vis.cpp
                             -------------------

  Convert a Lammps trajectory to a Pymol CGO trajectory

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
#include <cstring>

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
  
  // ----------------- Assign default atom type info
  unsigned max_atom_types=1000;
  vector<cPt> shapes;
  vector<colorPt> color_list;
  vector<double> alpha_list;
  vector<bool> spherical;
  shapes.assign(max_atom_types,cPt(0.5,0.5,0.5));
  color_list.assign(max_atom_types,colors["blue"]);
  alpha_list.assign(max_atom_types,1.0);
  spherical.assign(max_atom_types,true);

  // ----------------- Get the atom type info from a flavor file
  unsigned atom_type;
  cPt shape;
  string color;
  double alpha,diameter_x,diameter_y,diameter_z;
  unsigned num_types=0;
  char iline[5000];

  ifstream in;
  a::fileopen(in,cl.argstring(' ',0),error);
  while (!in.eof()) {
    in.getline(iline,500);
    if (!in || in.eof())
      break;
    istringstream line_in(iline);
    line_in >> atom_type;
    if (!line_in || line_in.eof())
      break;
    if (atom_type>max_atom_types)
      error.generate_error(0,"asphere_vis",
                    "asphere_vis will not handle atom types greater than "+
                     a::itoa(max_atom_types));
    line_in >> color >> alpha >> diameter_x;
    num_types++;
    if (!line_in)
      error.generate_error(0,"asphere_vis",
                           "Improperly formatted flavor_file at line:\n  "+
                           string(iline));
    color_list[atom_type]=colors[color];
    alpha_list[atom_type]=alpha;
    line_in >> diameter_y;
    if (!line_in) {
      shapes[atom_type].x=diameter_x/2.0;
      continue;
    }
    line_in >> diameter_z;
    if (!line_in)
      error.generate_error(0,"asphere_vis",
                           "Improperly formatted flavor_file at line:\n  "+
                           string(iline));
    shapes[atom_type]=cPt(diameter_x/2.0,diameter_y/2.0,diameter_z/2.0);
    spherical[atom_type]=false;    
  }
  error.note[9] << "Read in " << num_types << " atom types from flavor_file.\n";
  a::fileclose(in,error);
  
  // ----------------- Get the triangles per ellipsoid
  unsigned ellip_res=10;
  if (cl['r'])
    ellip_res=cl.argint('r',0);
  if (ellip_res==0) {
    error.addwarning(0,9,"asphere_vis","Cannot use -r 0. Setting to 10.");
    ellip_res=10;
  }
  
  a::fileopen(in,cl.argstring(' ',1),error);
  ofstream out;
  if (!cl['s']) {
    a::fileopen(out,cl.argstring(' ',2),error);
  }
  // ----------------- Get the bounding box
  GLSurface glb;
  bool bound_found=false;
  if (!cl['o']) {
    if (parse_to("ITEM: BOX BOUNDS",in)) {
      bound_found=true;
      cPt bound[2];
      in >> bound[0].x >> bound[1].x;
      in >> bound[0].y >> bound[1].y;
      in >> bound[0].z >> bound[1].z;
      Vertex v;
      v.transparency=1;
      v.valid_normal=false;
      v.color=colors["white"];
      for (unsigned i=0; i<2; i++)
        for (unsigned j=0; j<2; j++)
          for (unsigned k=0; k<2; k++) {
            v.cpt=cPt(bound[i].x,bound[j].y,bound[k].z);
            glb.addvertex(v);
          }
      glb.addline(0,1);
      glb.addline(0,2);
      glb.addline(0,4);
      glb.addline(1,3);
      glb.addline(1,5);
      glb.addline(2,3);
      glb.addline(2,6);
      glb.addline(3,7);
      glb.addline(4,5);
      glb.addline(4,6);
      glb.addline(5,7);
      glb.addline(6,7);
      if (!cl['s'])
        glb.writelines(out,"gridb");
      out << "cmd.set('cgo_dot_width',8)\n";
   } else
     error.addwarning(0,9,"asphere_vis",
          "Could not find ITEM: BOX BOUNDS in input file. No box output.");
  }
  
  if (!in)
    error.generate_error(0,"asphere_vis",
                         "Error reading ITEM: BOX BOUNDS.");
  a::fileclose(in,error);
  
  // ----------------- Generate the frames
  unsigned frame=0;
  unsigned max_frame=std::numeric_limits<unsigned>::max();
  if (cl['f'])
    max_frame=cl.argint('f',0);
    
  a::fileopen(in,cl.argstring(' ',1),error);

  // ----------------- Get to the start frame
  while (frame<start_frame)
    if (!parse_to("ITEM: ATOMS",in))
      error.generate_error(0,"asphere_vis",
                      "Could not find first frame in interval in dump file.");
    else
      frame++;

  unsigned wrote=0;
  Vertex v;      
  while (true) {
    // ----------------- Get the atom count
    unsigned atom_count;
    if (parse_to("ITEM: NUMBER OF ATOMS",in))
      in >> atom_count;
    else
      break;
    if (!in)
      error.generate_error(0,"asphere_vis",
                         "Error reading ITEM: NUMBER OF ATOMS in input file.");

    if (frame>end_frame)
       break;
    if (!parse_to("ITEM: ATOMS",in))
      break;
    GLSurface gls;
    unsigned id, atom_type;
    cPt atom_center;
    Quaternion q;
    for (unsigned i=0; i<atom_count; i++) {
      in.getline(iline,5000);
      istringstream line_in(iline);
      line_in >> id >> atom_type;
      if (atom_type>max_atom_types)
        error.generate_error(0,"asphere_vis",
                      "asphere_vis will not handle atom types greater than "+
                       a::itoa(max_atom_types));
      
      line_in >> atom_center;
      if (spherical[atom_type]==false)
        line_in >> q;
      if (!line_in) {
        error.addwarning(0,9,"asphere_vis","Error reading frame: "+
                         a::itoa(frame));
        break;
      }
      if (spherical[atom_type]) {
        v.cpt=atom_center;
        v.color=color_list[atom_type];
        v.transparency=alpha_list[atom_type];
        gls.addvertex(v);
        gls.add_sphere(gls.size_vertices()-1,shapes[atom_type].x);
      } else {
        gls.add_ellipsoid(atom_center,shapes[atom_type],q,color_list[atom_type],
                          alpha_list[atom_type],ellip_res);
      }
    }
    if (!in)
      break;
    if (cl['s']) {
      if (cl['b'])
        fi.set_file_num(frame);
      a::fileopen(out,fi.nextfilename(),error);
      if (!cl['o']) {
        glb.writelines(out,"gridb");
        out << "cmd.set('cgo_dot_width',8)\n";
      }
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
  string gridversion=a::strcenter("Graphics Library Version 0.2",70);
 out << endl << progname << endl << gridversion << endl
   << "______________________________________________________________________\n"
   << a::strcenter("W. Michael Brown",70) << endl
   << a::strcenter("1/12/2007",70) << endl
   << "______________________________________________________________________\n"
   << "Tool for LAMMPS trajectory visualization in Pymol.\n\n"
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
  cl.addargname(' ',"flavor_file");
  cl.addargname(' ',"dump_file");
  cl.addargname(' ',"output_py_file");
  cl.add('f',1);
  cl.addargname('f',"max_frame");
  cl.adddescription('f',"Do not write more than max_frame frames to the output file.");
  cl.add('r',1);
  cl.addargname('r',"ellip_res");
  cl.adddescription('r',"Resolution of ellipsoids in PyMol. The number of triangles per ellipsoid is equal to 2*(ellip_res^2). Default is 10.");
  cl.add('s',0);
  cl.adddescription('s',"Output the results into separate .py files. The filename and extension for the output files is taken from output_py_file.");
  cl.add('i',3);
  cl.addargname('i',"start_frame");
  cl.addargname('i',"skip");
  cl.addargname('i',"end_frame");
  cl.adddescription('i',"Render the specified frame interval inclusive between start_frame and end_frame. skip gives the number of frames to skip between each rendered frame. A value of 0 outputs every frame between start_frame and end_frame. The first frame in the dump file is frame 0.");
  cl.add('b',0);
  cl.adddescription('b',"When used with -s, the option will number the filenames based on the frame number. By default, they are numbered consequtively from zero.");
  cl.add('o',0);
  cl.adddescription('o',"Do not output the outline for the simulation box.");
        
  // Stuff for every executable
  cl.addhelp('h',0);
  cl.adddescription('h',"Print out the man page for help");
  cl.add('n',1);
  cl.addargname('n',"notice_level");
  cl.adddescription('n',"Set the degree of program output.  Use: \n\n\t-n  0\tNo output\n\t-n 10\tNormal program output\n\t-n 20\tParameters useful for reproducing the results\n\t-n 30\tAll output");

  // Short Description
  cl.addtoman_chapter("NAME","Tools for ellipsoid visualization in PyMol of a LAMMPS trajectory.");

  // Version
  cl.addtoman_chapter("VERSION","Version 0.2");

  // Full Description
  const string desc[22]={
    "Tool for converting LAMMPS trajectories into compiled graphics objects for ",
    "visualization in PyMol. The flavor_file is an input file that describes ",
    "the color, transparency, and size/shape of each atom type. The flavor_file ",
    "consists of two possible line formats. For spherical particles, the format ",
    "is:\n\n",
    "\t\\fIatom_type color alpha diameter\\fR\n\n",
    "where alpha is used to adjust the transparency of the particle. For ",
    "ellipsoidal particles, the format is:\n\n"
    "\t\\fIatom_type color alpha diameter_x diameter_y diameter_z\\fR\n\n",
    "Ellipsoidal and spherical line formats can be mixed in the same flavor_file ",
    "For any atom type not listed in the flavor_file a blue sphere of size 1 is ",
    "assumed.\n\n",
    "The dump_file is a LAMMPS trajectory. For atom types specified as spherical ",
    "in the flavor_file, the dump_file must contain \\fItag type x y z\\fR as ",
    "the first columns. For atom types specified as ellipsoidal in the ",
    "flavor_file, the columns are \\fItag type x y z quatw quati quatj quatk\\fR.",
    "The latter can be gerenated, for example, with the ",
    "LAMMPS dump_style custom command with the following arguments in order:\n\n",
    "\t\\fItag type x y z quatw quati quatj quatk\\fR\n\n",
    "The output file is a python file for input to Pymol. This can be viewed ",
    "from the command line using \\fIpymol output.py\\fR or by using the \\fIrun\\fR ",
    "command from within Pymol."
  };
  cl.addtoman_chapter("DESCRIPTION",22,desc);

  Colors colors;
  cl.addtoman_chapter("AVAILABLE COLORS",colors.colorlist());
  
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
  cl.write_man_page(cout,"0.2","Graphics Utilities");
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
    if (strncmp(token,iline,strlen(token))==0)
      return true;
  }
}
