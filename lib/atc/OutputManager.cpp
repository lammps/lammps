#include <string>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include "OutputManager.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"

namespace ATC {

static const int kFieldPrecison = 12;
static const int kFieldWidth = kFieldPrecison + 6;

static const int kFileNameSize = 26; // HERE <<<<
static string tensor_component_names[9] = {"11","12","13",
                                           "21","22","23", 
                                           "31","32","33"}; 
static string sym_tensor_component_names[6] = {"11","22","33","12","13","23"}; 
static string vector_component_names[3] = {"_X","_Y","_Z"}; 
static string list_component_names[26] = {"_a","_b","_c","_d","_e","_f","_g","_h","_i","_j","_k","_l","_m","_n","_o","_p","_q","_r","_s","_t","_u","_v","_w","_x","_y","_z"}; 

string* get_component_names(int type) { // HERE <<<<
 string* componentNames = list_component_names;
 if (type==VECTOR_OUTPUT)
   componentNames = vector_component_names;
 else if (type == SYM_TENSOR_OUTPUT)
   componentNames = sym_tensor_component_names;
 else if (type == TENSOR_OUTPUT)
   componentNames = tensor_component_names;
 return componentNames;
}

             
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
OutputManager::OutputManager(string outputPrefix, set<int> & otypes)
  : initialized_(false),
    firstStep_(true),
    firstGlobalsWrite_(true),
    writeGlobalsHeader_(true),
    coordinates_(NULL),
    connectivities_(NULL),
    dataType_(POINT),
    outputPrefix_(outputPrefix),
    ensightOutput_(otypes.count(ENSIGHT)),
    textOutput_(otypes.count(GNUPLOT)),
    fullTextOutput_(otypes.count(FULL_GNUPLOT)),
    vtkOutput_(otypes.count(VTK)),
    tensorToComponents_(false), // paraview does not support tensors
    vectorToComponents_(false)
    
{}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
OutputManager::OutputManager()
  : initialized_(false),
    firstStep_(true),
    firstGlobalsWrite_(true),
    writeGlobalsHeader_(true),
    coordinates_(NULL),
    connectivities_(NULL),
    dataType_(POINT),
    outputPrefix_("NULL"),
    ensightOutput_(true),
    textOutput_(false),
    fullTextOutput_(false),
    vtkOutput_(false),
    tensorToComponents_(false), // paraview does not support tensors
    vectorToComponents_(false)
{}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
OutputManager::~OutputManager() {}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
void OutputManager::set_option(OutputOption option, bool value) {
  if      (option == OUTPUT_VECTOR_COMPONENTS) vectorToComponents_ = value;
  else if (option == OUTPUT_TENSOR_COMPONENTS) tensorToComponents_ = value;
  else throw ATC_Error("unsupported output option");
};

//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
void OutputManager::initialize(string outputPrefix, set<int> & otypes) 
{ 
  if (outputPrefix_ !=  outputPrefix ) { // new stream with existing object
    outputPrefix_ =  outputPrefix;
    initialized_ = false;
  }
  outputTimes_.clear();
  if (otypes.count(ENSIGHT) > 0) ensightOutput_ = true;
  else                           ensightOutput_ = false;
  if (otypes.count(GNUPLOT) > 0) textOutput_ = true;
  if (otypes.count(FULL_GNUPLOT) > 0) fullTextOutput_ = true;
  if (otypes.count(VTK) > 0) vtkOutput_ = true;
  firstStep_ = true;
  firstGlobalsWrite_ = true;
  writeGlobalsHeader_ = true;
}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
void OutputManager::print_custom_names() {
  map<string,vector<string> >::const_iterator itr;
  string msg = "output custom names:\n";
  for (itr = fieldNames_.begin(); itr != fieldNames_.end(); itr++) {
    string stem = itr->first;
    vector<string> names = itr->second;
    for (unsigned int i = 0; i < names.size(); i++) {
      msg+= stem+" : "+names[i]+"\n";
    }
  }
  ATC::LammpsInterface::instance()->print_msg_once(msg);
}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
// Dump text-based fields to disk for later restart
void OutputManager::write_restart_file(string fileName, RESTART_LIST *data)
{
  FILE * fp=NULL;
  fp=fopen(fileName.c_str(),"wb"); // open
  RESTART_LIST::iterator iter;
  for (iter = data->begin(); iter != data->end(); iter++) {
    const MATRIX* field_data = iter->second;
    for (int i = 0; i < field_data->nRows(); ++i) {
      for (int j = 0; j < field_data->nCols(); ++j) {
        double x = (*field_data)(i,j);
        fwrite(&x,sizeof(double),1,fp);
      }
    }
  }
  fclose(fp);
}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
// Read a file corresponding to a write by write_restart_file

void OutputManager::read_restart_file(string fileName, RESTART_LIST *data)
{
  FILE * fp=NULL;
  fp=fopen(fileName.c_str(),"rb"); // open
  RESTART_LIST::iterator iter;
  for (iter = data->begin(); iter != data->end(); iter++) {
    MATRIX* field_data = iter->second;
    for (int i = 0; i < field_data->nRows(); ++i) {
      for (int j = 0; j < field_data->nCols(); ++j) {
        double myVal;
        fread(&myVal,sizeof(double),1,fp);
        (*field_data)(i,j) = myVal;  
    
      }
    }
  }
  fclose(fp);
}
//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
void OutputManager::write_globals(void) 
{
  if ( outputPrefix_ == "NULL") return;
  string file = outputPrefix_ + ".GLOBALS";
  ofstream text;
  if ( firstGlobalsWrite_ ) text.open(file.c_str(),ios_base::out);
  else                      text.open(file.c_str(),ios_base::app);
  firstGlobalsWrite_ = false;

  map<string, double>::iterator iter;
  // header
  if ( firstStep_ || writeGlobalsHeader_) {
    text << "# time:1 ";
    int index = 2;
    for (iter = globalData_.begin(); iter != globalData_.end(); iter++) 
    {
      string name  = iter->first;
      string str; stringstream out; out << ":" << index++; 
      str = out.str();
      name.append(str);
      text.width(kFieldWidth); text << name << " ";
    }
    text << '\n';
  }
  writeGlobalsHeader_ = false;
  // data
  text.width(kFieldWidth); text << outputTimes_[outputTimes_.size()-1] << " ";
  for (iter = globalData_.begin(); 
       iter != globalData_.end(); iter++) {
    double value = iter->second;
    text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << value << " ";
  }
  text << "\n";
}
//-----------------------------------------------------------------------------
//* 
//-----------------------------------------------------------------------------
void OutputManager::write_geometry(const MATRIX *coordinates, 
                                   const Array2D<int> *connectivities)
{
  if ( outputPrefix_ == "NULL") throw ATC_Error( "No outputPrefix given."); 
 
  number_of_nodes_ = coordinates->nCols();
  coordinates_ = coordinates;
  connectivities_ = connectivities;
  if (ensightOutput_) write_geometry_ensight();
  if (textOutput_)    write_geometry_text();
  initialized_ = true;
}
//-----------------------------------------------------------------------------
//* 
//-----------------------------------------------------------------------------
void OutputManager::write_geometry_ensight(void)
{
  // geometry based on a reference configuration
  string geom_file_name = outputPrefix_ + ".geo";

  // open file
  FILE * fp=NULL;
  char buffer[80];
  if ( ! initialized_ ) {
    fp=fopen(geom_file_name.c_str(),"wb"); // open
    strcpy(buffer,"C Binary");
    fwrite(buffer,sizeof(char),80,fp);
  }
  else {
    fp=fopen(geom_file_name.c_str(),"ab"); // append
  }
  if (fp == NULL) {
    throw ATC_Error("can not create Ensight geometry file");
  }

  // write preamble
  strcpy(buffer,"BEGIN TIME STEP");
  fwrite(buffer,sizeof(char),80,fp);
  strcpy(buffer,"Ensight geometry file");
  fwrite(buffer,sizeof(char),80,fp);
  strcpy(buffer,"description");
  fwrite(buffer,sizeof(char),80,fp);
  strcpy(buffer,"node id assign");
  fwrite(buffer,sizeof(char),80,fp);
  strcpy(buffer,"element id assign");
  fwrite(buffer,sizeof(char),80,fp);

  // per part
  strcpy(buffer,"part");
  fwrite(buffer,sizeof(char),80,fp);
  int part_number=1;
  fwrite(&part_number,sizeof(int),1,fp);
  strcpy(buffer,"description");
  fwrite(buffer,sizeof(char),80,fp);
  const MATRIX & coordinates = *coordinates_;
  // write coordinates
  strcpy(buffer,"coordinates");
  fwrite(buffer,sizeof(char),80,fp);
  fwrite(&number_of_nodes_,sizeof(int),1,fp);
   
  int number_of_spatial_dimensions = coordinates.nRows();
  if (number_of_spatial_dimensions != 3) 
    throw ATC_Error("Ensight writer needs a 3D geometry");

  for (int i = 0; i < number_of_spatial_dimensions; ++i) 
  {
    for (int j = 0; j < number_of_nodes_; ++j) 
    {
      float x = (float) coordinates(i,j);
      fwrite(&x,sizeof(float),1,fp);
    }
  }

  // write mesh connectivities or point "connectivities"
  if (connectivities_) 
  { 
    dataType_ = MESH;
    int nodes_per_element = connectivities_->nRows();
    if      (nodes_per_element == 4) { strcpy(buffer,"tetra4"); }
    else if (nodes_per_element == 8) { strcpy(buffer,"hexa8"); }
    else if (nodes_per_element == 20) { strcpy(buffer,"hexa20"); }
    else if (nodes_per_element == 27) { strcpy(buffer,"hexa27"); }
    else 
      throw ATC_Error("Ensight writer does not recoginize element type");
    fwrite(buffer,sizeof(char),80,fp);
    int number_of_elements = connectivities_->nCols();
    fwrite(&number_of_elements,sizeof(int),1,fp);
    int number_of_nodes_per_element = connectivities_->nRows();
    for (int j = 0; j < number_of_elements; ++j) 
    {
      for (int i = 0; i < number_of_nodes_per_element; ++i) 
      {
        int inode = (*connectivities_)(i,j) +1; // 1 based numbering
        fwrite(&inode,sizeof(int),1,fp);
      }
    }
  } 
  else 
  {
    strcpy(buffer,"point");
    fwrite(buffer,sizeof(char),80,fp);
    int number_of_elements = number_of_nodes_;
    fwrite(&number_of_elements,sizeof(int),1,fp);
    for (int j = 0; j < number_of_elements; ++j)  
    {
      int inode = j +1; // 1 based numbering
      fwrite(&inode,sizeof(int),1,fp);
    }
  }

  // end per part
  strcpy(buffer,"END TIME STEP");
  fwrite(buffer,sizeof(char),80,fp);
  fclose(fp);
}
//-----------------------------------------------------------------------------
//* 
//-----------------------------------------------------------------------------
void OutputManager::write_geometry_text(void)
{
  if ( outputPrefix_ == "NULL") return;
  // geometry based on a reference configuration
  string geom_file_text = outputPrefix_ + ".XYZ";

  // open file
  ofstream text;
  text.open(geom_file_text.c_str(),ios_base::out);
  if (connectivities_) 
  { 
    int number_of_elements = connectivities_->nCols();
    int number_of_nodes_per_element = connectivities_->nRows();
    for (int j = 0; j < number_of_elements; ++j) 
    {
      text << "#";
      for (int i = 0; i < number_of_nodes_per_element; ++i) 
      {
        int inode = (*connectivities_)(i,j) +1; // 1 based numbering
          text << setw(6) << inode;
      }
      text << "\n";
    }
  }

  const MATRIX & coordinates = *coordinates_;
  int number_of_spatial_dimensions = coordinates.nRows();
  for (int j = 0; j < number_of_nodes_; ++j) 
  {
    text << setw(6) << j+1 << " ";
    for (int i = 0; i < number_of_spatial_dimensions; ++i) 
    {
      text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << coordinates(i,j) << " ";
    }
    text << "\n";
  }
  text << "\n";
}
//-----------------------------------------------------------------------------
/** pack "soln" into data */
//-----------------------------------------------------------------------------
void OutputManager::write_data(double time, FIELDS *soln, OUTPUT_LIST *data, const int *node_map) 
{
  // pack
  OUTPUT_LIST combined_data;  
  if (soln) 
  {
    FIELDS::iterator iter;
    for (iter = soln->begin(); iter != soln->end(); iter++) 
    {
      FieldName field_index = iter->first;
      MATRIX* field_data = &((iter->second).set_quantity());  
      string field_name = field_to_string(field_index);
      combined_data[field_name] = field_data;
    }
  }
  if (data) 
  {
    OUTPUT_LIST::iterator iter;
    for (iter = data->begin(); iter != data->end(); iter++) 
    {
      string field_name = iter->first;
      const MATRIX* field_data = iter->second;
      combined_data[field_name] = field_data;
    }
  }
  write_data(time, &(combined_data), node_map);
};
//-----------------------------------------------------------------------------
/** write data  */
//-----------------------------------------------------------------------------
void OutputManager::write_data(double time, OUTPUT_LIST *data, const int *node_map)
{
  if (! initialized_) {
    throw ATC_Error("must write geometry before data");
  }

  // store the time step value
  outputTimes_.push_back(time);

  if (ensightOutput_) {
    // write data
    OUTPUT_LIST::iterator iter;
    for (iter = data->begin(); iter != data->end(); iter++) 
    {
      string field_name = iter->first;
      const MATRIX* field_data = iter->second;
      write_data_ensight(field_name, field_data, node_map);
    }
    // write dictionary
    write_dictionary(time,data);
  }

  // write text dump
  if (textOutput_) {
    write_data_text(data); 
    if (firstStep_ && node_map) {
      string map_file_text = outputPrefix_ + ".MAP";
      ofstream text;
      text.open(map_file_text.c_str(),ios_base::out);
      for (int i=0; i< number_of_nodes_ ; i++) {
        text << node_map[i] << "\n";
      }
      text.close();
    }
  }
  else if (fullTextOutput_) {
    write_data_text(data,node_map); 
  }
  if (vtkOutput_) {
    write_data_vtk(data); 
  }

  // global variables
  if (! globalData_.empty()) write_globals();

  if (firstStep_) firstStep_ = false;
}
//-----------------------------------------------------------------------------
/** write (ensight gold format "C" binary) data  */
// use "ens_checker" to check binary format
//-----------------------------------------------------------------------------
void OutputManager::write_data_ensight(string field_name, const MATRIX *field_data, const int *node_map)
{
  int ndof = field_data->nCols();
  int col_start = 0;
  int col_end   = ndof;
  string filenames[kFileNameSize];
  int nfiles = 1;
  filenames[0] = outputPrefix_ + "." + field_name;
  int type = data_type(ndof);
  if (use_component_names(type)){
    nfiles = ndof;
    if (nfiles > kFileNameSize) {
      if (ATC::LammpsInterface::instance()->rank_zero()) {
        stringstream ss;
        ss << " only writing " << kFileNameSize
           << " components of " <<  field_name << " which has " << ndof;
        ATC::LammpsInterface::instance()->print_msg(ss.str());
      }
      nfiles =  kFileNameSize;
    }
    string* component_names = get_component_names(type);
    for (int ifile = 0; ifile < nfiles; ++ifile) 
    {
      string comp_name;
      if (! custom_name(field_name,ifile,comp_name))
        comp_name = field_name + component_names[ifile];
      filenames[ifile] = outputPrefix_ + "." + comp_name;
    }
  }

  for (int ifile = 0; ifile < nfiles; ++ifile) 
  {
    // for vector/tensor to components
    if ( nfiles > 1 ) 
    { 
      col_start = ifile;   
      col_end   = ifile+1;   
    }
 
    // open or append data file
    string data_file_name = filenames[ifile];
    FILE * fp=NULL;
    if ( outputTimes_.size() == 1 ) {
      fp=fopen(data_file_name.c_str(),"wb"); // open
    }
    else {
      fp=fopen(data_file_name.c_str(),"ab"); // append
    }
    if (fp == NULL) {
      throw ATC_Error("can not create Ensight data file: "+data_file_name);
    }
  
    // write data
    char buffer[80];
    strcpy(buffer,"BEGIN TIME STEP");
    fwrite(buffer,sizeof(char),80,fp);
    strcpy(buffer,"field name"); 
    fwrite(buffer,sizeof(char),80,fp);

    // per part
    strcpy(buffer,"part");
    fwrite(buffer,sizeof(char),80,fp);
    int part_number = 1;
    fwrite(&part_number,sizeof(int),1,fp);
    strcpy(buffer,"coordinates");
    fwrite(buffer,sizeof(char),80,fp);
    if (node_map) 
    {
      for (int j = col_start; j < col_end; ++j) 
      {
        for (int i = 0; i < number_of_nodes_; ++i) 
        {
          int inode = node_map[i];
          float x = (float) (*field_data)(inode,j);
          fwrite(&x,sizeof(float),1,fp);
        }
      }
    } 
    else 
    {
      for (int j = col_start; j < col_end; ++j)  
      {
        for (int i = 0; i < field_data->nRows(); ++i) 
        {
          float x = (float) (*field_data)(i,j);
          fwrite(&x,sizeof(float),1,fp);
        }
      }
    }
    // end per part

    strcpy(buffer,"END TIME STEP");
    fwrite(buffer,sizeof(char),80,fp);
    fclose(fp);
  }
}

//-----------------------------------------------------------------------------
/** write data dict for both text & full_text */
//-----------------------------------------------------------------------------
void OutputManager::write_text_data_header(OUTPUT_LIST *data, ofstream & text, int k)
{
    if (data) 
    {
      OUTPUT_LIST::iterator iter;
      for (iter = data->begin(); iter != data->end(); iter++) 
      {
        string field_name = iter->first;
        int nrows = iter->second->nRows();
        if (!(nrows>0)) {
          string msg = field_name + " does not have data for output";
          throw ATC_Error(msg);
        }
        int ncols = iter->second->nCols();
        if (ncols > kFileNameSize) {
          if (ATC::LammpsInterface::instance()->rank_zero()) {
            stringstream ss;
            ss << " only writing " << kFileNameSize
               << " components of " <<  field_name << " which has " << ncols;
            ATC::LammpsInterface::instance()->print_msg(ss.str());
          }
          ncols =  kFileNameSize;
        }
        if (ncols == 1) {
          string name = field_name;
          custom_name(field_name,0,name);
          string str; stringstream out; out <<":"<<k; str = out.str();
          name.append(str);
          text.width(kFieldWidth); text << name << " ";
          k++;
        }
        else {
          for (int i = 1; i <= ncols; i++) {
            string name = field_name;
            string str; stringstream out; 
            if (! custom_name(field_name,i-1,name)) { out <<"_"<<i; }
            out <<":"<<k; str = out.str();
            name.append(str);
            text.width(kFieldWidth); text << name << " ";
            k++;
          }
        }
      }
    } else { throw ATC_Error(" data missing from output");}
    text << "\n";
}

//-----------------------------------------------------------------------------
/** write data  in text format */
//-----------------------------------------------------------------------------
void OutputManager::write_data_text(OUTPUT_LIST *data)
{
  string data_file_text = outputPrefix_ + ".DATA";
  ofstream text;
  if  (firstStep_) text.open(data_file_text.c_str(),ios_base::out);
  else text.open(data_file_text.c_str(),ios_base::app);

  // write data label header
  if (firstStep_) 
  {
    text.width(6); text << "# index:1" << " "; // give an ordinate for gnuplot
    text.width(10); text << " step:2" << " "; 
    write_text_data_header(data,text,3);
  }
  text << "# timestep " << outputTimes_.size() << " : "
       << outputTimes_[outputTimes_.size()-1]  << "\n";
  
  int nrows = 0;
  OUTPUT_LIST::iterator iter;
  iter = data->begin();
  if (iter == data->end()) { throw ATC_Error(" no data in output");}
  const MATRIX* field_data = iter->second;
  nrows = field_data->nRows();

  for (int i = 0; i < nrows; ++i) 
  {
    text.width(6); text << i << " "; // give an ordinate for gnuplot
    text.width(10); text << outputTimes_.size() << " "; 
    OUTPUT_LIST::iterator iter;
    for (iter = data->begin(); iter != data->end(); iter++) 
    {
      const MATRIX* field_data = iter->second;
      int ncols = field_data->nCols();
      if (ncols > kFileNameSize) { ncols = kFileNameSize;}
      for (int j = 0; j < ncols; ++j) 
      {
        text.width(kFieldWidth); 
        text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << (*field_data)(i,j) << " ";
      }
    }
    text <<"\n";
  }
  text <<"\n";
}

//-----------------------------------------------------------------------------
/** write data  in full text format */
//-----------------------------------------------------------------------------
void OutputManager::write_data_text(OUTPUT_LIST *data, const int *node_map)
{
  string data_file_text = outputPrefix_ + ".DATA";
  ofstream text;
  if  (firstStep_) text.open(data_file_text.c_str(),ios_base::out);
  else text.open(data_file_text.c_str(),ios_base::app);

  // write data label header
  if (firstStep_) 
  {
    text.width(6); text << "# index:1" << " "; 
    text.width(6); text << " id:2" << " "; 
    text.width(10); text << " step:3" << " "; 
    text.width(4); text << " x:4" << " "; 
    text.width(4); text << " y:5" << " "; 
    text.width(4); text << " z:6" << " "; 
    write_text_data_header(data,text,7);

    if (connectivities_) 
    { 
      int number_of_elements = connectivities_->nCols();
      int number_of_nodes_per_element = connectivities_->nRows();
      text << "# connectivities  number_of_elements: " << number_of_elements 
           << " nodes_per_element: " << number_of_nodes_per_element << "\n"; 
      for (int j = 0; j < number_of_elements; ++j) 
      {
        text << "#";
        for (int i = 0; i < number_of_nodes_per_element; ++i) 
        {
          int inode = (*connectivities_)(i,j) +1; // 1 based numbering
              text << setw(6) << inode;
        }
        text << "\n";
      }
    }
  }
  text << "# timestep " << outputTimes_.size() << " : "
       << outputTimes_[outputTimes_.size()-1]  << "\n";
  
  OUTPUT_LIST::iterator iter;
  iter = data->begin();
  if (iter == data->end()) { throw ATC_Error(" no data in output");}
  int nnodes = coordinates_->nCols();

  for (int i = 0; i < nnodes; ++i) 
  {
    int unode = i;
    if (node_map) unode = node_map[i];
    text.width(6); text << i << " "; 
    text.width(6); text << unode << " "; 
    text.width(10); text << outputTimes_.size() << " "; 
    // coordinates
    for (int j = 0; j < coordinates_->nRows(); ++j) {
      text.width(kFieldWidth); 
      text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << (*coordinates_)(j,i) << " ";
    }
    // data
    OUTPUT_LIST::iterator iter;
    for (iter = data->begin(); iter != data->end(); iter++) 
    {
      const MATRIX* field_data = iter->second;
      int ncols = field_data->nCols();
      if (ncols > kFileNameSize) { ncols = kFileNameSize; }
      for (int j = 0; j < ncols; ++j) 
      {
        text.width(kFieldWidth); 
        text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << (*field_data)(unode,j) << " ";
      }
    }
    text <<"\n";
  }
  text <<"\n";
}

//-----------------------------------------------------------------------------
/** write data  in vtk text format */
//-----------------------------------------------------------------------------
void OutputManager::write_data_vtk(OUTPUT_LIST *data)
{
  string data_file_text = outputPrefix_ + ".vtk";
  ofstream text;
  if  (firstStep_) text.open(data_file_text.c_str(),ios_base::out);
  else throw ATC_Error(" vtk format can not handle multiple steps");
  text << "# vtk DataFile Version 3.0\n";
  text << "# " << outputPrefix_ << "\n";
  text << "ASCII\n";
  text << "DATASET UNSTRUCTURED_GRID\n";
  // geometry
  int nnodes = coordinates_->nCols();
  text << "POINTS " << nnodes << " float\n";
  for (int i = 0; i < nnodes; ++i) {
    for (int j = 0; j < coordinates_->nRows(); ++j) {
      text.width(kFieldWidth); 
      text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << (*coordinates_)(j,i) << " ";
    }
    text << "\n";
  }
  text << "\n";
  int nelems = connectivities_->nCols();
  int nodes_per_element = connectivities_->nRows();
  text << "CELLS " << nelems << " " << nelems*(nodes_per_element+1) << "\n";
  for (int j = 0; j < nelems; ++j) {
    text << setw(6) << nodes_per_element;
    for (int i = 0; i < nodes_per_element; ++i) {
       int inode = (*connectivities_)(i,j); // 0 based numbering
       text << setw(6) << inode;
    }
    text << "\n";
  }
  text << "\n";
  int cell_type = 4 ; 
  text << "CELL_TYPES " << nelems << "\n";
  for (int j = 0; j < nelems; ++j) {
    text << cell_type << "\n";
  }
  text << "\n";
  // data
  text << "POINT_DATA " << nnodes << "\n";
  text << "\n";
  OUTPUT_LIST::iterator iter;
  for (iter = data->begin(); iter != data->end(); iter++) 
  {
    string field_name = iter->first;
    const MATRIX* field_data = iter->second;
    int ncols = field_data->nCols();
    if (ncols == 1) { 
      text << "SCALARS " << field_name << " float\n";
      text << "LOOKUP_TABLE default\n";
    }
    else {
      text << "VECTORS " << field_name << " float\n";
    }
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < ncols; ++j) {
        text.width(kFieldWidth); 
        text << setw(kFieldWidth) << std::scientific << std::setprecision(kFieldPrecison) << (*field_data)(i,j) << " ";
      }
      text <<"\n";
    }
  }
  text <<"\n";
}

/** write (ensight gold : ASCII "C" format) dictionary */
void OutputManager::write_dictionary(double time, OUTPUT_LIST *data)
{
  // file names
  string dict_file_name = outputPrefix_ + ".case";
  string geom_file_name = outputPrefix_ + ".geo";
  
  // open file
  FILE * fp=NULL;
  if ((fp=fopen(dict_file_name.c_str(),"w")) == NULL) 
  {
    throw ATC_Error("can not create Ensight case file");
  }

  // write file
  fprintf(fp,"FORMAT\n");
  fprintf(fp,"type: ensight gold\n");
  fprintf(fp,"GEOMETRY\n");
  if ( dataType_ == POINT) {
    fprintf(fp,"model: 1 1 %s change_coords_only\n", geom_file_name.c_str());
  } else {
    fprintf(fp,"model: %s\n", geom_file_name.c_str());
  }
  fprintf(fp,"VARIABLE\n");

  // data types
  if (!data) throw ATC_Error("no data for output");
  OUTPUT_LIST::iterator iter;
  int ncols = 0;
  for (iter = data->begin(); iter != data->end(); iter++) {
    string field_name = iter->first;
    string field_file = outputPrefix_ + "." + field_name;
    const MATRIX* field_data = iter->second;
    int fieldCols = field_data->nCols();
    ncols += fieldCols;
    int type = data_type(fieldCols);
    if (use_component_names(type)){
      string* component_names = get_component_names(type);
      int ndof = fieldCols;
      if (ndof > kFileNameSize) ndof = kFileNameSize;
      for (int j = 0; j < ndof; ++j) 
      {
        string comp_name;
        if (! custom_name(field_name,j,comp_name))
        comp_name = field_name + component_names[j];
        string comp_file = outputPrefix_ + "." + comp_name;
        fprintf(fp,"scalar per node: 1 1 %s %s\n",
          comp_name.c_str(),comp_file.c_str());
      }
    }
    else if (type == VECTOR_OUTPUT) {
      fprintf(fp,"vector per node: 1 1 %s %s\n",
        field_name.c_str(),field_file.c_str());
    } 
    else if (type == SYM_TENSOR_OUTPUT) {
      fprintf(fp,"tensor symm per node: 1 1 %s %s\n",
        field_name.c_str(),field_file.c_str());
    } 
    else if (type == TENSOR_OUTPUT) {
      fprintf(fp,"tensor asymm per node: 1 1 %s %s\n",
        field_name.c_str(),field_file.c_str());
    } 
    else {
      fprintf(fp,"scalar per node: 1 1 %s %s\n",
        field_name.c_str(),field_file.c_str());
    }
  }

  if (!firstStep_ && ncols != nDataCols_) {
    throw ATC_Error("number of columns of data has changed: start new output");
  }
  nDataCols_ = ncols;

  int nsteps = outputTimes_.size();
  fprintf(fp,"TIME\n");
  fprintf(fp,"time set: 1\n");
  fprintf(fp,"number of steps: %10d\n",nsteps);
  if ( dataType_ == POINT) {
    fprintf(fp,"filename start number: 0\n");
    fprintf(fp,"filename increment: 1\n");
  }
  fprintf(fp,"time values:\n");
  for (int j = 0; j < nsteps; ++j) {
    double t = outputTimes_[j];
    fprintf(fp,"%12.5e",t);
    if ((j+1)%6 == 0) fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  fprintf(fp,"FILE\n");
  fprintf(fp,"file set: 1\n");
  fprintf(fp,"number of steps: %10d\n",nsteps);
  fclose(fp);
};

} // end ATC namespace


