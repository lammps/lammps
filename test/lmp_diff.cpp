/**************************************************************************
                                  mdiff.cpp
                             -------------------

    Numerical diff for regression testing

  __________________________________________________________________________
    Miscellaneous Utilities
  __________________________________________________________________________

    begin      : Thu Feb 14 2008
    authors    : W. Michael Brown
    email      : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc!=4) {
    std::cerr << "lmp_diff numbers file1 file2\n";
    std::cerr << "lmp_diff all file1 file2\n";
    std::cerr << "lmp_diff time file1 file2\n";
    return 1;
  }    
	
  bool time=false;
  if (strcmp(argv[1],"time")==0)
    time=true;
  else if (strcmp(argv[1],"numbers")!=0) {
    std::cerr << "lmp_diff numbers file1 file2\n";
    std::cerr << "lmp_diff time file1 file2\n";
    return 1;
  }    

  ifstream in1, in2;
  in1.open(argv[2]);
  in2.open(argv[3]);
  if (!in1 || !in2) {
    cout << "ERROR " << argv[2] << " " << argv[3] << " "
         << "Files Does Not Exist\n";
    return 1;
  }
	
  bool ok=true;
  double max_error=0;
  double mean_error=0;
  unsigned number_tokens=0;
  double fastest=0,slowest=0;
  
  char cline[256];
  while (true) {
    string token1,token2;
    in1 >> token1;
    in2 >> token2;

    if (in1.eof() || !in1) {
      if (in2.eof() || !in2)
        break;
      ok=false;
      break;
    }
    
    if (in2.eof() || !in2) {
      ok=false;
      break;
    }
    
    if (token1=="LAMMPS" || token1=="by" || token1=="Memory" ||
        token1=="FFT" || (token1=="time" && time==false)) {
      in1.getline(cline,256);
      in2.getline(cline,256);
      continue;
    }

    if (token1=="Nlocal:" || token2=="Nlocal:") {
      while (token1!="Dangerous") {
        in1 >> token1;
        if (in1.eof() || !in1) {
          ok=false;
          break;
        }
      }
      while (token2!="Dangerous") {
        in2 >> token2;
        if (in2.eof() || !in2) {
          ok=false;
          break;
        }
      }
      continue;
    }
        
    bool fix1=false,fix2=false;
    if (token1=="fix") {
      fix1=true;
      in1.getline(cline,256);
    }
    if (token2=="fix") {
      fix2=true;
      in2.getline(cline,256);
    }
    if (fix1 && fix2)
      continue;
    else if (fix1)
      in1 >> token1;
    else if (fix2)
      in2 >> token2;
    
    bool package1=false,package2=false;
    if (token1=="using") {
      package1=true;
      in1.getline(cline,256);
    }
    if (token2=="using") {
      package2=true;
      in2.getline(cline,256);
    }
    if (package1 && package2)
      continue;
    else if (package1)
      in1 >> token1;
    else if (package2)
      in2 >> token2;

    if (token1=="package") {
      package1=true;
      in1.getline(cline,256);
    }
    if (token2=="package") {
      package2=true;
      in2.getline(cline,256);
    }
    if (package1 && package2)
      continue;
    else if (package1)
      in1 >> token1;
    else if (package2)
      in2 >> token2;
    
    if (time) {
      if (token1=="time") {
        in1 >> token1;
        in2 >> token2;
        if (token1=="of") {
          in1 >> token1;
          in2 >> token2;
        } else {
          continue;
        }
      } else
        continue;
    }
      
    if (token1==token2)
      continue;
      
    if (strstr(token1.c_str(),"gpu")!=NULL || 
        strstr(token2.c_str(),"gpu")!=NULL)
      continue;
      
    istringstream sin1, sin2;
    double num1, num2;
    sin1.str(token1);
    sin2.str(token2);
    sin1 >> num1;
    sin2 >> num2;
    if (!sin1 || !sin2) {
      ok=false;
      break;
    }
    double denr=num1;
    if (num1==0 || num2==0)
      denr=1;
    double errorv=fabs((num1-num2)/denr);
    max_error=max(errorv,max_error);
    mean_error+=errorv;
    number_tokens++;
    if (time) {
      if (num2<num1) {
        double ft=(num1-num2)/num1*100.0;
        if (ft>fastest) fastest=ft;
      } else {
        double ft=(num2-num1)/num1*100.0;
        if (ft>slowest) slowest=ft;
      }
    }
  }
  
  if (time) {
    cout << "TIME " << argv[2] << " " << argv[3] << " FASTER: "
         << fastest << "% SLOWER: " << slowest << "%\n";
    return 0;
  }

  if (ok==false) {
    cout << "ERROR " << argv[2] << " " << argv[3] << " "
         << "Files Differ\n";
    return 0;
  }
  
  if (number_tokens>0)
    mean_error/=number_tokens;
    
  if (max_error<0.01)
    cout << "OK ";
  else
    cout << "ERROR ";
  
  cout << argv[2] << " " << argv[3] << " MAX: "
       << max_error*100.0 << "% MEAN: " << mean_error*100.0
       << "%\n";

  return 0;
}
