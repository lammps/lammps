#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include "ATC_TypeDefs.h"

#include <map>
#include <string>


// 1 -> scalar
// 3 -> vector  x,y,z
// NOT 6 -> tensor  xx,xy,xz,yy,yz,zz
// 6 -> tensor  xx,yy,zz,xy,zx,yz
// 9 -> tensor  xx,xy,xz,yx,yy,yz,zx,zy,zz

using namespace std;

namespace ATC {

  enum OutputType     { ENSIGHT=0, GNUPLOT, FULL_GNUPLOT, VTK };
  enum OutputDataType { POINT=0, MESH };
  enum OutputDataCardinality { SCALAR_OUTPUT=0, VECTOR_OUTPUT, TENSOR_OUTPUT, 
    SYM_TENSOR_OUTPUT, LIST_OUTPUT };
  enum OutputOption   { OUTPUT_VECTOR_COMPONENTS=0, OUTPUT_TENSOR_COMPONENTS};

  /**
   *  @class  OutputManager 
   *  @brief  Base class for handling output desired from an AtC computation 
   */

  class OutputManager{

  public:
    OutputManager(void);
    OutputManager(string outputPrefix, set<int> &otypes);
    ~OutputManager(void);

    /** initialize output */
    void initialize(string outputPrefix, set<int> &otypes);

    /** set output options */
    void set_option(OutputOption option, bool value); 

    // Dump text-based field info to disk for later restart
    void write_restart_file(string fileName, RESTART_LIST *data);

    // Read text-based field file written from write_restart_file
    void read_restart_file(string fileName, RESTART_LIST *data);

    /** write initial/reference geometry
        default is to write point data, 
        if connectivities are given then mesh data will be output 
        coordinates : num _total_ points/nodes X num spatial dim 
        connectivities : num elements X num nodes per element*/
    void write_geometry(const MATRIX *coordinates,
                        const Array2D<int> *connectivity=NULL);

    /** write data from a time step
        specify node_map to handle periodic soln & data */
    void write_data(double time, OUTPUT_LIST *data, const int *node_map=NULL);
    void write_data(double time, FIELDS *soln, OUTPUT_LIST *data, 
      const int *node_map=NULL);

    /** add custom names for any field */
    void add_field_names(const string& name, const vector<string>& list) {
      fieldNames_[name] = list; }
    /** add a scalar to a text output file */
    void add_global(const string& name, const double& value) {
      globalData_[name] = value; }

    /** delete a scalar from the output */
    void delete_global(const string& name) { globalData_.erase(name); }

    /** reset the stored output scalars */
    void reset_globals() { globalData_.clear(); writeGlobalsHeader_=true; }

    /** return data type: scalar, vector, tensor, list */
    int data_type(const DENS_MAT & data) const {
      return data_type(data.nCols());
    }
    int data_type(int cols) const {
      if      (cols == 1) return SCALAR_OUTPUT;
      else if (cols == 3) return VECTOR_OUTPUT;
      else if (cols == 6) return SYM_TENSOR_OUTPUT;
      else if (cols == 9) return TENSOR_OUTPUT;
      else                return LIST_OUTPUT;
    }
    bool use_component_names(int type) const {
      if ( (type==LIST_OUTPUT) ||     
       ((type==SYM_TENSOR_OUTPUT || type==TENSOR_OUTPUT) && tensorToComponents_)
      || (type==VECTOR_OUTPUT && vectorToComponents_) ) 
        return true;
      else 
        return false;
    }
    bool custom_name(const string field, const int index, string & name) const {
      map<string,vector<string> >::const_iterator itr = fieldNames_.find(field);
      if (itr == fieldNames_.end()) return false;
      vector<string>  names = itr->second;
      name = names[index];
      return true;
    }
    void print_custom_names(); 

  private:

    void write_geometry_ensight(void);
    void write_geometry_text(void);
    void write_data_ensight(string name, const MATRIX *data, const int *node_map);
    void write_text_data_header(OUTPUT_LIST *data, ofstream & text, int k);
    void write_data_text(OUTPUT_LIST *data);
    void write_data_text(OUTPUT_LIST *data, const int *node_map);
    void write_data_vtk(OUTPUT_LIST *data);
    void write_dictionary(double time, OUTPUT_LIST *data);
    void write_globals();

    /** status flags */
    bool initialized_, firstStep_, firstGlobalsWrite_, writeGlobalsHeader_;

    /** custom field names */
    map<string,vector<string> > fieldNames_;

    /**  pointers to mesh data */
    const MATRIX * coordinates_;
    const Array2D<int> * connectivities_;
    /** number of columns of data */
    int nDataCols_;
    /** number of nodes */ 
    int number_of_nodes_; 
    /** data type */
    int dataType_;
    /** base name for output files */
    string outputPrefix_;
    /** list of output timesteps */
    vector<double> outputTimes_;
    /** output type flags */
    bool ensightOutput_,textOutput_,fullTextOutput_,vtkOutput_;
    /** output tensor as its components */
    bool tensorToComponents_;
    /** output vector as its components */
    bool vectorToComponents_;
    /** global variables */
    map<string,double> globalData_;
  };
}
#endif
