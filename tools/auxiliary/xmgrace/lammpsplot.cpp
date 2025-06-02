// By Vikas Varshney @ vv0210 at gmail dot com

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

typedef vector<string> sentence;
typedef vector<sentence> filedata;
typedef string::size_type s_type;

vector<filedata> f;

class fd {
public:
  string filename;
  string property;
  vector <float> value;
};


/* This function splits the sentence into array of strings separated by " "*/
sentence split(const string& s) {
  sentence a;
  s_type i=0;
  while (i!=s.size()) {
    while (i!=s.size() && isspace(s[i])) {
      ++i;
    }
    
    s_type j=i;
    while (j!=s.size() && !isspace(s[j]))
      ++j;
    
    if (i!=j) {
      a.push_back(s.substr(i,j-i));
      i=j;
    }
  }
  return a;
}

/* This function add array of strings (passed as string) to a sentence */
sentence add_strings(sentence& temp_entry, const string l) {
  sentence temp=split(l); 
  s_type mm;
  if (temp.size()!=0)
    for (mm=0;mm!=temp.size();mm++)
      temp_entry.push_back(temp[mm]);
  else 
    return temp_entry;
  return temp_entry;
}

/* This function reads the logfile produced by lammps */
void readlogfile_multi(filedata& atemp, ifstream& input, int leave_index) {
  const string ll("----------------");
  const string step("Step");
  const string Mem("Memory");
  const string Loop("Loop");
  int leave=0;
  string line;
  s_type kk=0,mm,kkk;
  s_type length;
  int style;
  int halffile=0;
  while (getline(input,line)) { 
    sentence entry=split(line);
    if (entry.size()!=0) 
      // Find keyword Memory
      if (!entry[0].compare(Mem)) {
	getline(input,line);
	entry=split(line);
	// Style = 1 means multi line format and style = 2 means single line format
	if (!entry[0].compare(ll)) style=1; else if (!entry[0].compare(step)) style=2; else cout << "style not matching" << endl;
 	switch (style) {
 	case 1: {
 	  // Read the line till ll comes up next and store it in a variable.
	  sentence zero;
  	  sentence one;
	  
 	  zero.push_back(entry[1]);
	  one.push_back(entry[2]);
	  getline(input,line);
 	  entry=split(line);
 	  do {
 	    if (entry.size()!=0)
 	      for(int i=0;i<entry.size();i+=3) {
 		zero.push_back(entry[i]);
 		one.push_back(entry[i+2]);
 	      } 
 	    getline(input,line);
	    entry=split(line);
	    if (!entry[0].compare(Loop)) 
	      break;
 	  } 
 	  while (entry[0].compare(ll));
 	  if (!atemp.size())
 	    atemp.push_back(zero);
 	  atemp.push_back(one);
 	  int length=zero.size();
 	  int t=0;

 	  // Read all the lines till Loop comes up.
 	  while (entry[0].compare(Loop)) {
	    t++;
	    sentence a1;
	    a1.push_back(entry[2]);
	    for(int j=0;j<(length+1)/3;j++) {
	      if (!getline(input,line)) {
		halffile=1; 
		break;
	      }
	      entry=split(line);
	      for(int i=0;i<entry.size();i+=3) 
		a1.push_back(entry[i+2]);
	    }
	    if (a1.size()==length)
	      atemp.push_back(a1);
	    if (halffile) 
	      break;
	    if (!getline(input,line)) 
	      break;
	    else
	      entry=split(line);
	  }
 	  break;
	}
 	case 2:
	  {
	    sentence zero;
 	    for(int i=0;i!=entry.size();i++)
	      zero.push_back(entry[i]);
	    getline(input,line);
	    entry=split(line);
	    if (!atemp.size())
	      atemp.push_back(zero);
	    int length=zero.size();
	    while(entry[0].compare(Loop)) {
	      sentence a1;
	      for(int i=0;i!=entry.size();i++) 
		a1.push_back(entry[i]);
	      atemp.push_back(a1);
	      getline(input,line);
	      entry=split(line);
	      if (entry.size()<length) break;
 	    }
	    break;
	  } 	
	}
	
	if (leave_index>leave) {
	  atemp.clear();
	  leave++;
	}
      }
  }
}

/* This function1 searches for a string in a sentence and return true/false*/
bool search_string(string str, sentence sen, int& j) {
  j=0;
  while (j!=sen.size()) {
    if (!sen[j].compare(str)) {
      break;
    }
    j++;
  }  
  if (j==sen.size()) return false; else return true; 
}

/* This function is used to ask user to input parameters to be plotted*/
sentence get_userinput() {
  string temp;
  sentence t;
  cout << "Enter the parameter(s) of which you want the plot for files, read successfully:" << endl; 
  getline(cin,temp);
  t=split(temp);
  if (!t[1].compare("cross"))
    if (t.size()!=4) {
      cout << "Currently only two cross parameters can be plotted \nThe syntax is s/p cross param1 param2\n" << endl;
      exit(1);
    }
  return t;
}

/*This function replaces certain input keywords to their appropriate words found in lammps file*/
void transform_userinput(sentence& temp) {
  if (temp.size()) {
    for(int i=1;i!=temp.size();i++) {
      if (!temp[i].compare("t")) temp[i].replace(0,temp[i].length(),"Temp");
      if (!temp[i].compare("p")) temp[i].replace(0,temp[i].length(),"Press");
      if (!temp[i].compare("v")) temp[i].replace(0,temp[i].length(),"Volume");
      if (!temp[i].compare("vd")) temp[i].replace(0,temp[i].length(),"E_vdwl");
      if (!temp[i].compare("te")) temp[i].replace(0,temp[i].length(),"TotEng");
      if (!temp[i].compare("ke")) temp[i].replace(0,temp[i].length(),"KinEng");
      if (!temp[i].compare("pe")) temp[i].replace(0,temp[i].length(),"PotEng");
      if (!temp[i].compare("be")) temp[i].replace(0,temp[i].length(),"E_bond");
      if (!temp[i].compare("ae")) temp[i].replace(0,temp[i].length(),"E_angle");
      if (!temp[i].compare("de")) temp[i].replace(0,temp[i].length(),"E_dihed");
      if (!temp[i].compare("ce")) temp[i].replace(0,temp[i].length(),"E_coul");
      if (!temp[i].compare("le")) temp[i].replace(0,temp[i].length(),"E_long");
    }
  }
}

/* This function searches if the parameter is present in the lammps file and return a array of 
   associated with the parameter stored in string */
map<string, vector<int> >  update2_userinput(sentence temp) {
  if ((!temp[0].compare("s")) || (!temp[0].compare("p"))) {
    s_type j;
    int jtemp;
    int notfind;
    int notfindall=1;
    map<string, vector<int> > ret;
    sentence::iterator tb=temp.begin();
    for(j=0;j!=f.size();j++) 
      ret[*tb].push_back(0);
    tb++;
     
    for(sentence::iterator aa=tb; aa!=temp.end(); aa++) {
      notfind=1;
      for(j=0;j!=f.size();j++) {
	if (search_string(*aa,f[j][0],jtemp)) {
	  ret[*aa].push_back(jtemp);
	  notfind=0;
	  notfindall=0;
	}
	else
	  ret[*aa].push_back(0);
      }
      if (notfind) {
	temp.erase(aa);
	aa--;
      }
    }
    
    if (notfindall) {
      cerr << "\nDid not find any matching parameters in any input files, read successfully. Please check!!!\n\n";
      exit(1);
    }
    return ret;
  }
  else {
    cerr << "\nFirst string should be s (for series) or p (for parallel)\n\n" << endl; 
    exit(1);
  }
}

/* This function just print out a line with string str on standard output*/
void print_line(string str) {
  for(int i=0;i<100;i++)
    cout << str;
  cout << endl;
}

/*This function reads the data from the input files*/
void read_input(int argc, char** argv, sentence& sf,string& analysisfile) {
  filedata a; 
  int success_filecounter=0; 
  int leave_index=0;
  // Read the data into the file.
  for (int i=1; i<argc; ++i) {	
    if (!strcmp(argv[i],"-leave")) {
      leave_index=atoi(argv[i+1]);
      i++;
      continue;
    }

     if (!strcmp(argv[i],"-analysis")) {
      analysisfile=argv[i+1];
      i++;
      continue;
    }

    ifstream in(argv[i]);
    if (in) {
      sf.push_back(argv[i]);
      cout << "Reading file    : " << argv[i] << endl;
      if (a.size())
	a.clear();
      readlogfile_multi(a,in,leave_index);
      f.push_back(a);
      in.close();
      success_filecounter++;
    } else {
      cerr << "Cannot open file: " << argv[i] << endl;  
    }
  }

  if(!success_filecounter)
    exit(1);
  else {
    sentence aa;
    int exists;
    cout << endl;
    cout << "Select parameters to plot from following: --- " << endl;
    for (int i=0;i<success_filecounter;i++) {
      for(int j=0;j<f[i][0].size();j++) {
	exists=0;
	for(int k=0;k!=aa.size();k++)
	  if (!f[i][0][j].compare(aa[k])) {
	    exists=1;
	    break;
	  }       
	if (!exists)
	  aa.push_back(f[i][0][j]);
      }
    }

    for(int k=0;k!=aa.size();k++) {
      if (!aa[k].compare("Temp")) aa[k].replace(0,aa[k].length(),"Temp(t)");
      if (!aa[k].compare("Press")) aa[k].replace(0,aa[k].length(),"Press(p)");
      if (!aa[k].compare("Volume")) aa[k].replace(0,aa[k].length(),"Volume(v)");     
      if (!aa[k].compare("E_vdwl")) aa[k].replace(0,aa[k].length(),"E_vdwl(vd)");     
      if (!aa[k].compare("TotEng")) aa[k].replace(0,aa[k].length(),"TotEng(te)");     
      if (!aa[k].compare("KinEng")) aa[k].replace(0,aa[k].length(),"KinEng(ke)");     
      if (!aa[k].compare("PotEng")) aa[k].replace(0,aa[k].length(),"PotEng(pe)");     
      if (!aa[k].compare("E_bond")) aa[k].replace(0,aa[k].length(),"E_bond(be)");     
      if (!aa[k].compare("E_angle")) aa[k].replace(0,aa[k].length(),"E_angle(ae)");     
      if (!aa[k].compare("E_dihed")) aa[k].replace(0,aa[k].length(),"E_dihed(de)");     
      if (!aa[k].compare("E_coul")) aa[k].replace(0,aa[k].length(),"E_coul(ce)");     
      if (!aa[k].compare("E_long")) aa[k].replace(0,aa[k].length(),"E_long(le)");
      cout << aa[k];
      if (aa[k].size()<8) cout << "\t";
      cout << "\t";
      if (k%6==5) cout << endl;
    }
    cout << "\n--------------\n\n";
    cout << "Please enter your parameters in the following format: --- " << endl;
    cout << "1. First string should be 'p' or 's' for plotting parameters from different files in parallel or series" << endl;
    cout << "2. Please enter one or more parameters listed above after entering the first string." << endl;
    cout << "3. For certain parameters, short forms may be used as listed in bracket above" << endl;
    cout << "   ----- Example: To plot total energy, please type 'p te'." << endl;
    cout << "4. In order to plot one parameter with respect to other, please type 'cross' after p/s followed by two parameters" << endl;
    cout << "5. The independent parameter could be 'Step' as mentioned above." << endl;
    cout << "6. Currently, only two parameters can be plotted with respect to each other." << endl;
    cout << "   ----- Example: To plot temperature vs. pressure p cross t p" << endl ;
    cout << "--------------\n\n";
    cout << "For analysis, please exit and use -analysis flag with filename. The file should contain atleast one of these commands: --- " << endl;
    cout << "1. average     all/single_filename startvalue(start) endvalue(end)" << endl;
    cout << "2. scale       filename property factor" << endl;
    cout << "3. inverse     filename property" << endl;
    cout << "4. write       Under Coding" << endl;
    cout << "5. subtract    property1 property2 filename outputpropertyname" << endl;
    cout << "6. add         Under Coding" << endl;
    cout << "--------------\n\n";
  }
}

string print_space(int t) {
  string f;
  for (int i=0;i!=t;i++)
    f.push_back(' ');
  return f;
}


/*This function prints out the parameter associated index in successfully read file */
void print_file_index_output(sentence s, map<string, vector<int> > x, sentence sf) {
  size_t maxlen=0;
  for(int i=0;i!=sf.size();i++)
    maxlen=max(maxlen,sf[i].length());
  cout << print_space(maxlen+11);
  for(sentence::iterator aa=++s.begin(); aa!=s.end(); aa++) 
    cout << *aa << "\t";
  cout << endl;

  for(int j=0;j!=f.size();j++) {
    cout << "File-> " << sf[j] << print_space(maxlen-sf[j].length()+1) <<":  ";
    for(sentence::iterator aa=++s.begin(); aa!=s.end(); aa++) {
      cout << x[*aa][j] << "\t";
    }
    cout << endl;
  }
}

void plot_xmgrace(sentence s, vector <fd> data) {
  ofstream out("temp");
  // Some string definitions
  string s1 = "@ s";
  string s2 = "   comment ";
  string s3 = "   legend  ";
  int plotiterator=0;  int set;  int style;
  int half_data=data.size()/2;

  // Printing out X and Y axis properties
  out << "@    xaxis  label char size 1.500000\n";
  out << "@    xaxis  ticklabel char size 1.250000\n";
  if (s[1].compare("cross")) 
    out << "@    xaxis label " << "\"time\""<< endl;
  else
    out << "@    xaxis label " << "\"" << s[2] << "\""<< endl;
  out << "@    yaxis  label char size 1.500000\n";
  out << "@    yaxis  ticklabel char size 1.250000\n";

  // Printing out legend properties
  out << "@    legend box linestyle 0\n";
  out << "@    legend 0.25, 0.8\n"; 

  // Selecting style for series or parallel plotting.
  if (!s[0].compare("s")) style=1; else style=2;

  if (s[1].compare("cross")) {
    if (s.size()==2) {
      out << "@ yaxis label " << "\"";
      for(int i=1;i!=s.size();i++) 
	out << s[i] << " ";
      out << "\""<< endl;
    }
  }
  else
    out << "@    yaxis label " << "\"" << s[3] << "\""<< endl; 
  int length=data.size();
  if (!s[1].compare("cross")) length/=2;
  // Storing data in file "temp" as asked (series or parallel)
  for(int i=0;i!=length;i++) {
    switch (style) {
    case 1: { // Series Case
      if (s[1].compare("cross")) { // This is executed if not cross
	out << s1 << plotiterator << s2 << "\""<< data[i].property << "\"" << endl;
	out << s1 << plotiterator << s3 << "\""<< data[i].property << "\"" << endl;
      }
      else
	{
	  out << s1 << plotiterator << s2 << "\""<< data[i+half_data].property << " vs " << data[i].property << "\"" << endl;
	  out << s1 << plotiterator << s3 << "\""<< data[i+half_data].property << " vs " << data[i].property << "\"" << endl;
	}
      out << s1 << plotiterator << " symbol 1" << endl;
      out << s1 << plotiterator << " symbol size 0.50000" << endl;
      out << s1 << plotiterator << " symbol color 1" << endl;
      out << s1 << plotiterator << " symbol pattern 1" << endl;
      out << s1 << plotiterator << " symbol fill color "<< plotiterator+1 << endl;
      out << s1 << plotiterator << " symbol fill pattern 1" << endl;
      out << s1 << plotiterator << " symbol linewidth 1.0" << endl;
      out << s1 << plotiterator << " symbol linestyle 1" << endl;
      out << s1 << plotiterator << " symbol char 65" << endl;
      out << s1 << plotiterator << " symbol char font 0" << endl;
      plotiterator++;
      for(int j=0;j!=data[i].value.size(); j++)  
	if (s[1].compare("cross"))  
	  out << data[i].value[j] << endl; // This is executed if not cross
	else
	  out << data[i].value[j] << " " << data[i+half_data].value[j]<< endl; // This is executed if cross
      out << "&" << endl;
      break;
    }
      
    case 2: { // Parallel Case 
      if (s[1].compare("cross")) {
	out << s1 << plotiterator << s2 << "\""<< data[i].property << ": " << data[i].filename << "\"" << endl;
	if (s.size()<=2)
	  out << s1 << plotiterator << s3 << "\"" << data[i].filename << "\"" << endl;
	else
	  out << s1 << plotiterator << s3 << "\""<< data[i].property << ": " << data[i].filename << "\"" << endl;
      }
      else {
	out << s1 << plotiterator << s2 << "\""<< data[i+half_data].property << " vs " << data[i].property << ": " << data[i].filename << "\"" << endl;
	out << s1 << plotiterator << s3 << "\""<< data[i+half_data].property << " vs " << data[i].property << ": " << data[i].filename << "\"" << endl;
      }
      out << s1 << plotiterator << " symbol 0" << endl;
      out << s1 << plotiterator  << " symbol size 0.50000" << endl;
      out << s1 << plotiterator  << " symbol color 1" << endl;
      out << s1 << plotiterator  << " symbol pattern 1" << endl;
      out << s1 << plotiterator  << " symbol fill color "<< plotiterator+1 << endl;
      out << s1 <<  plotiterator  << " symbol fill pattern 1" << endl;
      out << s1 <<  plotiterator  << " symbol linewidth 1.0" << endl;
      out << s1 <<  plotiterator  << " symbol linestyle 1" << endl;
      out << s1 <<  plotiterator  << " symbol char 65" << endl;
      out << s1 <<  plotiterator  << " symbol char font 0" << endl;
      plotiterator++;
      for(int j=0;j!=data[i].value.size(); j++)  
	if (s[1].compare("cross"))  {
	  out << data[i].value[j] << endl; // This is executed if not cross
        }
	else
	  out << data[i].value[j] << " " << data[i+half_data].value[j]<< endl; // This is executed if cross
      out << "&" << endl;
      break;
    }
    }
  }
  int rv=system("xmgrace temp &");
}
void plot_xmgrace_cross(sentence s, vector <fd> data) {
  ofstream out("temp");
  // Some string definitions
  string s1 = "@ s";
  string s2 = "   comment ";
  string s3 = "   legend  ";
  int plotiterator=0;  int set;  int style;
  int half_data=data.size()/2;
  // Printing out X and Y axis properties
  out << "@    xaxis  label char size 1.500000\n";
  out << "@    xaxis  ticklabel char size 1.250000\n";
  out << "@    xaxis label " << "\"" << s[2] << "\""<< endl;
  out << "@    yaxis  label char size 1.500000\n";
  out << "@    yaxis  ticklabel char size 1.250000\n";
  out << "@    yaxis label " << "\"" << s[3] << "\""<< endl;
 
  // Printing out legend properties
  out << "@    legend box linestyle 0\n";
  out << "@    legend 0.25, 0.8\n"; 

  // Selecting style for series or parallel plotting.
  if (!s[0].compare("s")) style=1; else style=2;

  // Storing data in file "temp" as asked (series or parallel)
  for(int i=0;i!=data.size()/2;i++) {
    switch (style) {
    case 1: { // Series Case
      out << s1 << plotiterator << s2 << "\""<< data[i+half_data].property << " vs " << data[i].property << "\"" << endl;
      out << s1 << plotiterator << s3 << "\""<< data[i+half_data].property << " vs " << data[i].property << "\"" << endl;
      out << s1 << plotiterator << " symbol 1" << endl;
      out << s1 << plotiterator << " symbol size 0.50000" << endl;
      out << s1 << plotiterator << " symbol color 1" << endl;
      out << s1 << plotiterator << " symbol pattern 1" << endl;
      out << s1 << plotiterator << " symbol fill color "<< plotiterator+1 << endl;
      out << s1 << plotiterator << " symbol fill pattern 1" << endl;
      out << s1 << plotiterator << " symbol linewidth 1.0" << endl;
      out << s1 << plotiterator << " symbol linestyle 1" << endl;
      out << s1 << plotiterator << " symbol char 65" << endl;
      out << s1 << plotiterator << " symbol char font 0" << endl;
      plotiterator++;
      for(int j=0;j!=data[i].value.size(); j++) 	
	out << data[i].value[j] << " " << data[i+half_data].value[j]<< endl;
      out << "&" << endl;
      break;
    }
      
    case 2: { // Parallel Case 
      out << s1 << plotiterator << s2 << "\""<< data[i+half_data].property << " vs " << data[i].property << ": " << data[i].filename << "\"" << endl;
      out << s1 << plotiterator << s3 << "\""<< data[i+half_data].property << " vs " << data[i].property << ": " << data[i].filename << "\"" << endl;
      out << s1 << plotiterator << " symbol 1" << endl;
      out << s1 << plotiterator  << " symbol size 0.50000" << endl;
      out << s1 << plotiterator  << " symbol color 1" << endl;
      out << s1 << plotiterator  << " symbol pattern 1" << endl;
      out << s1 << plotiterator  << " symbol fill color "<< plotiterator+1 << endl;
      out << s1 <<  plotiterator  << " symbol fill pattern 1" << endl;
      out << s1 <<  plotiterator  << " symbol linewidth 1.0" << endl;
      out << s1 <<  plotiterator  << " symbol linestyle 1" << endl;
      out << s1 <<  plotiterator  << " symbol char 65" << endl;
      out << s1 <<  plotiterator  << " symbol char font 0" << endl;
      plotiterator++;
      
      for(int j=0;j!=data[i].value.size(); j++) 	
	out << data[i].value[j] << " " << data[i+half_data].value[j]<< endl;
      out << "&" << endl;
      break;
    }
    }
  }
  int rv=system("xmgrace temp &");
}

vector<fd> convert_to_float(sentence s, map<string, vector<int> > x, sentence sf) {
  vector<fd> datatemp;
  int style,set;
  
  if (!s[0].compare("s")) style=1; else style=2;
  for(int i=1;i!=s.size();i++) {
    fd temp;
    switch (style) {
    case 1: {  // This is series case
      set=0;
      for(int j=0;j!=f.size(); j++) {
	if (x[s[i]][j]) {
	  if (!set) {
	    temp.filename="Allfiles";
	    temp.property=s[i];
	    set=1;
	  }
	  for(int tt=1;tt!=f[j].size();tt++) {
	    float t=atof(&f[j][tt][x[s[i]][j]][0]);
	    temp.value.push_back(t);
	  }
	}
	//	datatemp.push_back(temp);
      }
      datatemp.push_back(temp);
      break;
    }
    case 2: { // This is parallel case.
      for(int j=0;j!=f.size(); j++) {
	if (x[s[i]][j]) {
	  temp.filename=sf[j];
	  temp.property=s[i];
	  temp.value.clear();
	  for(int tt=1;tt!=f[j].size();tt++) {
	   float t= 0.0;
            if ((f[j][tt][x[s[i]][j]].compare("nan")))
	      t = atof(&f[j][tt][x[s[i]][j]][0]);
	    temp.value.push_back(t);
	  }
	  datatemp.push_back(temp);
	}
	else if (!s[i].compare("Step")) {
	  temp.filename=sf[j];
	  temp.property=s[i];
	  temp.value.clear();
	  for(int tt=1;tt!=f[j].size();tt++) {
	    float t=atof(&f[j][tt][x[s[i]][j]][0]);
	    temp.value.push_back(t);
	  }
	  datatemp.push_back(temp);
	}
      }
      break;
    }
    }
  }
  return datatemp;
}

void do_averageall(vector<fd> x,int start,int end) {
  for(int i=0;i!=x.size();i++) {
    double average=0;
    double stdev=0;
    int finish;
    if (end==0) 
      finish=x[i].value.size();
    else 
      finish=end;
    
    for (int j=start;j!=finish;j++) 
      average+=x[i].value[j];
    average/=(double)(finish-start);
    for (int j=start;j!=finish;j++) 
      stdev+=pow((x[i].value[j]-average),2);
    stdev/=(double)(finish-1-start);
    stdev=sqrt(stdev);
    cout << x[i].filename <<  ": " << x[i].property << "--- " << average << " " << stdev << " " << start << " " << finish << " Total Size " << x[i].value.size() << endl;
  }
}

void do_averagesingle(vector<fd> x,int start,int end, string str) {
  for(int i=0;i!=x.size();i++) {
    if (!x[i].filename.compare(str)) {
      float average=0;
      float stdev=0;
      int finish;
      if (end==0) finish=x[i].value.size();
      for (int j=start;j!=finish;j++) 
	average+=x[i].value[j];
      average/=(double)(finish-start);
      for (int j=start;j!=finish;j++) 
	stdev+=pow((x[i].value[j]-average),2);
      stdev/=(double)(finish-1-start);
      stdev=sqrt(stdev);
      cout << x[i].filename <<  ": " << x[i].property << "--- " << average << " " << stdev << " " << start << " " << finish << endl;
    }
  }
}

void do_scale_single(vector<fd> &x,float factor, string filename, string property) {
  for(int i=0;i!=x.size();i++) 
    if ((!x[i].filename.compare(filename)) && (!x[i].property.compare(property))) {
      fd temp;
      temp.filename=filename;
      temp.property=property+".scl";
      for (int j=0;j!=x[i].value.size();j++)
	temp.value.push_back(x[i].value[j]*factor);
      x.push_back(temp);
    }
}

void do_inverse_single(vector<fd> &x, string filename, string property) {
  for(int i=0;i!=x.size();i++) 
    if ((!x[i].filename.compare(filename)) && (!x[i].property.compare(property))) {
      fd temp;
      temp.filename=filename;
      temp.property=property+".inv";
      for (int j=0;j!=x[i].value.size();j++)
	temp.value.push_back(1.0/x[i].value[j]);
      x.push_back(temp);
    }
}

void do_write_data(vector<fd> x, sentence line) {
  //   int startvalue,endvalue;
  //   if (!line_data[1].compare("start")) startvalue=0; else startvalue=atoi(&line_data[2][0]);
  //   if (!line_data[2].compare("end")) endvalue=0; else endvalue=atoi(&line_data[3][0]);
  //   for(int j=3;j!=line.size();j++)
  //     for(int i=0;i!=x.size();i++) 
  //       if (!x[i].property.compare(line[i])) {
  // 	// Create a filename
  // 	// open a file.
  // 	// Write the data
  //       }
}

void do_subtract(vector<fd>& x, sentence line) {
  int indexfile1=-1;
  int indexfile2=-1;
  int minlength;
  for(int i=0;i!=x.size();i++) 
    if ((!x[i].property.compare(line[1])) && (!x[i].filename.compare(line[3])))
      indexfile1=i;
  for(int i=0;i!=x.size();i++) 
    if ((!x[i].property.compare(line[2])) && (!x[i].filename.compare(line[3])))
      indexfile2=i;
  if ((indexfile1>=0) && (indexfile2>=0)) {
  fd temp;
  temp.property=line[4];
  temp.filename=line[3];
  if (x[indexfile1].value.size()<x[indexfile2].value.size()) 
    minlength=x[indexfile1].value.size();
  else
    minlength=x[indexfile2].value.size();

  for (int j=0;j!=minlength;j++)
    temp.value.push_back(x[indexfile1].value[j]-x[indexfile2].value[j]);
  x.push_back(temp);    
  }
}

void do_add(vector<fd>& x, sentence line) {
  fd temp;
  int length=line.size();
  temp.property=line[length-1];
  temp.filename=line[length-2];
  vector <int> index;
  int minlength=10000000;

  for(int i=0;i!=length-3;i++)
    index.push_back(-1);
  for(int j=1;j!=length-2;j++)
    for(int i=0;i!=x.size();i++) 
      if ((!x[i].property.compare(line[j])) && (!x[i].filename.compare(line[length-2])))
	index[j-1]=i;

  for(int i=0;i!=length-3;i++) {
    if (index[i]>-1)
      if (minlength>x[index[i]].value.size())
	minlength=x[index[i]].value.size();
  }

  for (int j=0;j!=minlength;j++) {
    float value=0;
    for (int i=0;i!=length-2;i++)
      if (index[i]>0)
	value+=x[index[i]].value[j];
    temp.value.push_back(value);
  }
  x.push_back(temp);
}

void do_analysis(vector<fd> &data,string file) {
 ifstream inp_anal(file.c_str());
  if (inp_anal) {
    string line;
    while(getline(inp_anal,line)) {
      sentence line_data=split(line);
      
      // Average
      if (!line_data[0].compare("average")) {
	int startvalue,endvalue;
	if (!line_data[2].compare("start")) startvalue=0; else startvalue=atoi(&line_data[2][0]);
	if (!line_data[3].compare("end")) endvalue=0; else endvalue=atoi(&line_data[3][0]);
	if (!line_data[1].compare("all")) 
	  do_averageall(data,startvalue,endvalue);
	else 
	  do_averagesingle(data,startvalue,endvalue,line_data[1]);
      }

      // Scale
      if (!line_data[0].compare("scale")) {
	float factor=atof(&line_data[3][0]);
	do_scale_single(data,factor,line_data[1],line_data[2]);
      }

      // Inverse
      if (!line_data[0].compare("inverse")) 
	do_inverse_single(data,line_data[1],line_data[2]);
      
      // Write Data
      if (!line_data[0].compare("write"))
	do_write_data(data,line_data);

    // Add : Think about it and write
      if (!line_data[0].compare("write")) {
	do_add(data,line_data);
      } 
    // Subtract : Think about it and write
      if (!line_data[0].compare("subtract")) {
	cout << "coming here" << endl;
      	      do_subtract(data,line_data);
      }
      
      if (!line_data[0].compare("graph")) {
      }
      
  }
}
}

void print_start() {
int rv=system("clear ");
cout << endl;
cout << "Files Read: --- "  << endl;
}

int main (int argc, char **argv) {
  sentence s; 
  sentence stringfile;
  string anal_file="";
  vector <fd> data;  
  map<string, vector<int> > inp_index;  
  print_start();
  read_input(argc,argv,stringfile,anal_file);       // Read the input file
   s=get_userinput();                                // Read the parameters 
   transform_userinput(s);                           // Transform the parameters to correct names
   inp_index=update2_userinput(s);                   // Assign indices
  print_file_index_output(s,inp_index,stringfile);  // Print file index information
  data=convert_to_float(s,inp_index,stringfile);    // Convert the stored data in float formation
  if (anal_file.size())
    do_analysis(data,anal_file);                    // Do analysis if wanted.
  plot_xmgrace(s,data);                             // Plot the data in xmgrace
  return 0;
}
