// txt2html - written by Steve Plimpton, May 2004
// table formatting by Anna Reese, Jul 2004
// Sandia National Labs, www.cs.sandia.gov/~sjplimp
//
// txt2html converts a text file with simple formatting & markup into HTML
// formatting & markup specification is given in README
//
// Syntax: txt2html options file                read one file, write to stdout
//         txt2html optoins file1 file2 ...     read files, write files.html
//
// options:
//   -b = add a page-break comment to end of each HTML file
//     useful when set of HTML files will be converted to PDF
//   -x file = skip a file even if it appears in file list
//     specify full file name of input file
// input files are first opened as-is
//   if that fails a .txt suffix is added
// output files have an .html suffix added or replaced
//   (unless written to stdout)

#include <string>
#include <vector>
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#define MAXLINE 1024

// function prototypes

static int next_paragraph(FILE *fp, string &paragraph);
static int index_of_first_char_of_last_word(string &paragraph);
static void process_commands(int flag, string &s, string &pre, string &post);
static void substitute(string &s);
static string td_tag(int currentc);
static long find_n(string &s, int nend, int &n1);
static void file_open(int npair, string &infile, FILE **in, FILE **out);

// global variables for links, tables, lists, all command

vector<string> alias1;
vector<string> alias2;
int nlink;

int tableflag;       // makes a table if tb command specified
int rowquit;         // number of cols per row if c=N specified (default = 0)
string dwidth;       // width for all of the columns
string tabledelim;   // speciallized separator
string tablealign;   // alignment for the table as an image
string dataalign;    // alignment for data in table
string rowvalign;    // vertical alignment for table

int ncnum;             // # of columns with specified width
vector<int> cnum;      // column IDs
vector<string> cwidth; // column widths

int ncalign;              // # of columns with specified alignment
vector<int> acolnum;      // column IDs
vector<string> colalign ; // column alignment

int ncvalign;              // # of columns with specified vertical alignment
vector<int> vacolnum;      // column IDs
vector<string> colvalign ; // column vertical alignment

string listflag;
string allflag;

// main program

int main(int narg, char **arg)
{
  int npair;
  size_t n;
  string *infile;
  FILE *in,*out;
  int style,ifirst,ilast;
  string raw,pre,post,body,commands,final;

  // parse command-line options and args
  // setup list of files to process
  // npair = # of files to process
  // infile = input file names

  if (narg == 1) {
    fprintf(stderr,"Syntax: txt2html options file\n");
    fprintf(stderr,"        txt2html options file1 file2 ...\n");
    exit(1);
  }

  int breakflag = 0;
  int nskip = 0;
  char **skipfiles = NULL;

  int iarg = 1;
  while (arg[iarg][0] == '-') {
    if (strcmp(arg[iarg],"-b") == 0) breakflag = 1;
    else if (strcmp(arg[iarg],"-x") == 0) {
      skipfiles = (char **) realloc(skipfiles,(nskip+1)*sizeof(char *));
      n = strlen(arg[iarg+1]) + 1;
      skipfiles[nskip] = new char[n];
      strcpy(skipfiles[nskip],arg[iarg+1]);
      nskip++;
      iarg++;
    } else {
      fprintf(stderr,"Syntax: txt2html options file\n");
      fprintf(stderr,"        txt2html options file1 file2 ...\n");
      exit(1);
    }
    iarg++;
  }

  if (narg-iarg == 1) {
    npair = 1;
    infile = new string[npair];
    infile[0] = arg[narg-1];
  } else {
    npair = narg-iarg;
    infile = new string[npair];
    for (int i = 0; i < npair; i++) infile[i] = arg[i+iarg];
  }

  // loop over files

  for (int ipair = 0; ipair < npair; ipair++) {

    // skip file if matches -x switch

    int flag = 0;
    for (int i = 0; i < nskip; i++)
      if (strcmp(infile[ipair].c_str(),skipfiles[i]) == 0) flag = 1;
    if (flag) continue;

    // clear global variables before processing file

    alias1.clear();
    alias2.clear();
    nlink = 0;
    tableflag = 0;
    listflag = "";
    allflag = "";

    // open files & message to screen

    file_open(0,infile[ipair],&in,&out);
    fprintf(stderr,"Converting %s ...\n",infile[ipair].c_str());

    // scan file for link definitions
    // read file one paragraph at a time
    // process commands, looking only for link definitions

    while ((style = next_paragraph(in,raw))) {

      if (style == 2) {
        int n = index_of_first_char_of_last_word(raw);
        commands = raw.substr(n+1);
        process_commands(0,commands,pre,post);
      }

      raw.erase();
    }

    // close & reopen files

    fclose(in);
    file_open(npair,infile[ipair],&in,&out);

    // write leading <HTML>

    fprintf(out,"<HTML>\n");

    // process entire file
    // read file one paragraph at a time
    // delete newlines when line-continuation char at end-of-line
    // process commands for each paragraph
    // substitute text for each paragraph
    // write HTML to output file

    int rstflag = 0;

    while ((style = next_paragraph(in,raw))) {
      
      if (rstflag && raw.find("END_RST -->") != string::npos) {
        rstflag = 0;
        raw.erase();
        continue;
      } else if (rstflag == 0 && raw.find("<!-- RST") != string::npos) {
        rstflag = 1;
        raw.erase();
        continue;
      } else if (rstflag) {
        raw.erase();
        continue;
      }

      n = raw.find("\\\n");
      while (n < string::npos) {
        raw.erase(n,2);
        n = raw.find("\\\n");
      }

      ifirst = raw.find_first_not_of(" \t\n");
      ilast = raw.find_last_not_of(" \t\n");

      pre.erase();
      post.erase();

      if (raw[ifirst] == '<' && raw[ilast] == '>') {
        body = raw;
      } else if (style == 1) {
        body = raw;
        commands = "p\n";
        process_commands(1,commands,pre,post);
        substitute(body);
      } else {
        int n = index_of_first_char_of_last_word(raw);
        body = raw.substr(0,n) + "\n";
        commands = raw.substr(n+1);
        process_commands(1,commands,pre,post);
        substitute(body);
      }

      final = pre + body + post;
      fprintf(out,"%s\n",final.c_str());

      raw.erase();
    }

    // write trailing </HTML>

    if (breakflag) fprintf(out,"<!-- PAGE BREAK -->\n");
    fprintf(out,"</HTML>\n");

    // close files

    fclose(in);
    if (out != stdout) fclose(out);
  }

  // clean up memory
  
  for (int i = 0; i < nskip; i++) delete [] skipfiles[i];
  if (skipfiles) free(skipfiles);
  delete [] infile;
}

// return next paragraph as string
// discard leading blank lines
// paragraph is terminated by:
//   EOF or blank line or line ending with command that starts with ":"
// return 0 if EOF and no paragraph
// return 1 if no trailing command
// return 2 if trailing command

int next_paragraph(FILE *fp, string &paragraph)
{
  char *ptr;
  char str[MAXLINE];
  int first = 1;
  int len = 0;

  while (1) {
    ptr = fgets(str,MAXLINE,fp);
    if (ptr == NULL && first) return 0;
    if (ptr == NULL) return 1;
    len = strlen(str);
    if (len == MAXLINE-1) {
      fprintf(stderr,"ERROR: File has too-long a string - increase MAXLINE\n");
      exit(1);
    }

    // check for valid 7-bit ascii characters
    bool nonascii = false;
    for (int i=0; i < len; ++i) {
      char c = str[i];
      if (c != '\n' && c != '\t' && (c < ' ' || c > '~'))
        nonascii = true;
    }
    if (nonascii)
      fprintf(stderr,"WARNING: Non-portable characters in line: %s",ptr);

    if (strspn(str," \t\n") == strlen(str) && first) continue;
    if (strspn(str," \t\n") == strlen(str)) return 1;
    first = 0;

    paragraph += str;
    if (paragraph[index_of_first_char_of_last_word(paragraph)] == ':')
      return 2;
  }
}

// return index of first char in last word of paragraph string

int index_of_first_char_of_last_word(string &paragraph)
{
  size_t n = paragraph.find_last_not_of(" \t\n");
  size_t m = paragraph.find_last_of(" \t\n",n);
  if (m == string::npos) return 0;
  else return m+1;
}

// apply commands one after the other to the paragraph

void process_commands(int flag, string &s, string &pre, string &post)
{
  size_t start,stop,last;
  int narg;
  string command;
  vector<string> arg;

  start = 0;
  last = s.find_last_not_of(" \t\n");
  if (last == string::npos) return;

  while (start <= last) {

    // grab a single command with optional arguments
    // command = name of command
    // narg = # of args
    // arg = list of argument strings

    stop = s.find_first_of(",( \t\n",start);
    if (s[stop] == '(') {
      command = s.substr(start,stop-start);
      start = stop+1;
      narg = 0;
      while (1) {
        stop = s.find_first_of(",)",start);
        if (stop == string::npos) {
          fprintf(stderr,"ERROR: No trailing parenthesis in %s\n",s.c_str());
          exit(1);
        }
        arg.resize(narg+1);
        arg[narg] = s.substr(start,stop-start);
        narg++;
        start = stop+1;
        if (s[stop] == ')') {
          start++;
          break;
        }
      }
    } else {
      command = s.substr(start,stop-start);
      start = stop+1;
      narg = 0;
    }

    // if only in scan mode, just operate on link command

    if (flag == 0) {
      if (command == "link" && narg == 2) {
        //      s.erase(s.length()-1,1);
        for (int i = 0; i < nlink; i++)
          if (alias1[i] == arg[0]) {
            fprintf(stderr,"ERROR: Link %s appears more than once\n",
                    arg[0].c_str());
            exit(1);
          }
        alias1.resize(nlink+1);
        alias2.resize(nlink+1);
        alias1[nlink] = arg[0];
        alias2[nlink] = arg[1];
        nlink++;
      } else continue;
    }

    // process the command

    if (command == "line") {
      pre.append("<HR>");
    } else if (command == "p") {
      pre.append("<P>");
      post.insert(0,"</P>");
    } else if (command == "pre") {
      pre.append("<PRE>");
      post.insert(0,"</PRE>");
    } else if (command == "c") {
      pre.append("<CENTER>");
      post.insert(0,"</CENTER>");
    } else if (command == "h1") {
      pre.append("<H1>");
      post.insert(0,"</H1>");
    } else if (command == "h2") {
      pre.append("<H2>");
      post.insert(0,"</H2>");
    } else if (command == "h3") {
      pre.append("<H3>");
      post.insert(0,"</H3>");
    } else if (command == "h4") {
      pre.append("<H4>");
      post.insert(0,"</H4>");
    } else if (command == "h5") {
      pre.append("<H5>");
      post.insert(0,"</H5>");
    } else if (command == "h6") {
      pre.append("<H6>");
      post.insert(0,"</H6>");
    } else if (command == "b") {
      post.insert(0,"<BR>");
    } else if (command == "ulb") {
      pre.append("<UL>");
    } else if (command == "ule") {
      post.insert(0,"</UL>");
    } else if (command == "olb") {
      pre.append("<OL>");
    } else if (command == "ole") {
      post.insert(0,"</OL>");
    } else if (command == "dlb") {
      pre.append("<DL>");
    } else if (command == "dle") {
      post.insert(0,"</DL>");
    } else if (command == "l") {
      pre.append("<LI>");
    } else if (command == "dt") {
      pre.append("<DT>");
    } else if (command == "dd") {
      pre.append("<DD>");
    } else if (command == "ul") {
      listflag = command;
      pre.append("<UL>");
      post.insert(0,"</UL>");
    } else if (command == "ol") {
      listflag = command;
      pre.append("<OL>");
      post.insert(0,"</OL>");
    } else if (command == "dl") {
      listflag = command;
      pre.append("<DL>");
      post.insert(0,"</DL>");
    } else if (command == "link") {
      if (narg == 1) {
        string aname = "<A NAME = \"" + arg[0] + "\"></A>";
        pre.append(aname);
      }
    } else if (command == "image") {
      if (narg == 1) {
        string img = "<IMG SRC = \"" + arg[0] + "\">";
        pre.append(img);
      } else if (narg == 2) {
        string img = "<A HREF = \"" + arg[1] + "\">" +
          "<IMG SRC = \"" + arg[0] + "\">" + "</A>";
        pre.append(img);
      }
    } else if (command == "tb") {   // read the table command and set settings
      tableflag = 1;

      string tableborder = "1";    // these are the table defaults
      rowquit = 0;
      tablealign = "c";
      dataalign = "0";
      rowvalign = "0";

      ncnum = 0;
      ncalign = 0;
      ncvalign = 0;

      cnum.clear();
      acolnum.clear();
      vacolnum.clear();

      cwidth.clear();
      colalign.clear();
      colvalign.clear();

      tabledelim = ",";
      string tw = "";
      dwidth = "0";

      for (int i = 0; i < narg; i++) {     // loop through each tb() arg
        int tbstop;
        string tbcommand;
        tbstop = 0;
        tbstop = arg[i].find("=");
        tbcommand = arg[i].substr(0,tbstop);
        int n = arg[i].length();
        if (tbstop == -1) {
          continue;
        } else if (tbcommand == "c") {
          string collumn= arg[i].substr (tbstop+1,n-(tbstop+1));
          rowquit = atoi(collumn.c_str());
        } else if (tbcommand == "s") {
          tabledelim= arg[i].substr (tbstop+1,n-(tbstop+1));
        } else if (tbcommand == "b") {
          tableborder= arg[i].substr (tbstop+1,n-(tbstop+1));
        } else if (tbcommand == "w") {
          string width = "0";       
          if (arg[i].substr (n-1,1) == "%") {
            string width = arg[i].substr (tbstop+1,n-(tbstop+1));
            tw = " WIDTH=\"" + width + "\"";
          } else
            dwidth = arg[i].substr (tbstop+1,n-(tbstop+1));
        } else if (tbcommand == "ea") {
          dataalign= arg[i].substr (tbstop+1,n-(tbstop+1));
        } else if (tbcommand == "eva") {
          rowvalign= arg[i].substr (tbstop+1,n-(tbstop+1));
        } else if (tbcommand == "a") {
          tablealign= arg[i].substr (tbstop+1,n-(tbstop+1));
        } else if (tbcommand.substr(0,2) == "cw") {
          string cwnum= tbcommand.substr(2,tbstop-1);
          cnum.resize(ncnum+1);
          cnum[ncnum] = atoi(cwnum.c_str());
          cwidth.resize(ncnum+1);
          cwidth[ncnum]= arg[i].substr(tbstop+1,n-(tbstop+1));
          ncnum++;
        } else if (tbcommand.substr(0,2) ==  "ca") {
          string canum= tbcommand.substr(2,tbstop-1);
          acolnum.resize(ncalign+1);
          acolnum[ncalign] = atoi(canum.c_str());
          colalign.resize(ncalign+1);
          colalign[ncalign]= arg[i].substr(tbstop+1,n-(tbstop+1));
          ncalign++;
        } else if (tbcommand.substr(0,3) ==  "cva") {
          string cvanum= tbcommand.substr(2,tbstop-1);
          vacolnum.resize(ncvalign+1);
          vacolnum[ncvalign] = atoi(cvanum.c_str());
          colvalign.resize(ncvalign+1);
          colvalign[ncvalign]= arg[i].substr(tbstop+1,n-(tbstop+1));
          ncvalign++;
        } else {
          fprintf(stderr,
                  "ERROR: Unrecognized table command %s\n",tbcommand.c_str());
          exit(1);
        }
        
        tbstop = s.find("=");
      }
      
      string align;
      if (tablealign=="c") align="center";
      else if (tablealign=="r") align="right ";
      else if (tablealign=="l") align="left  ";
      else align="center";
      string tablea = "<DIV ALIGN=" + align + ">" ;
      pre.append(tablea);
      pre.append("<TABLE ");
      pre.append(tw);
      string border=" BORDER=" + tableborder + " >\n";
      pre.append(border);
      post.insert(0,"</TD></TR></TABLE></DIV>\n"); 
      
    } else if (command == "all") {
      if (narg == 1) allflag = arg[0];
      
    } else {
      fprintf(stderr,"ERROR: Unrecognized command %s\n",command.c_str());
      exit(1);
    }
  }
}
  
// perform substitutions within text of paragraph

void substitute(string &s)
{
  size_t n,m,p;
  char c;
  string text,link,href;
  string punctuation = ".,?!;:()";

  // substitute for bold & italic markers
  // if preceded by \ char, then leave markers in text

  n = s.find_first_of("[]{}");
  while (n != string::npos) {
    c = s[n];
    if (n > 0 && s[n-1] == '\\') s.erase(n-1,1);
    else {
      s.erase(n,1);
      if (c == '[') s.insert(n,"<B>");
      else if (c == ']') s.insert(n,"</B>");
      else if (c == '{') s.insert(n,"<I>");
      else if (c == '}') s.insert(n,"</I>");
    }
    n = s.find_first_of("[]{}",n);
  }

  // substitute for links

  n = s.find("\"_");
  while (n != string::npos) {
    m = s.rfind("\"",n-1);
    if (m == string::npos) {
      fprintf(stderr,"ERROR: Could not find matching \" for \"_ in %s\n",
              s.c_str());
      exit(1);
    }

    p = s.find_first_of(" \t\n",n) - 1;
    if (p == string::npos) {
      fprintf(stderr,"ERROR: Could not find end-of-link in %s\n",s.c_str());
      exit(1);
    }
    while (s.find_first_of(".,?!;:()",p) == p) p--;

    text = s.substr(m+1,n-m-1);
    link = s.substr(n+2,p-n-1);
    for (int i = 0; i < nlink; i++)
      if (alias1[i] == link) {
        link = alias2[i];
        break;
      }

    s.erase(m,p-m+1);
    href = "<A HREF = \"" + link + "\">" + text + "</A>";
    s.insert(m,href);
    n = s.find("\"_");
  }

  // format the paragraph as a table
  
  if (tableflag) {
    tableflag = 0;

    string DT;

    // set up <TR> tag
    // alignment for data in rows

    string tbalign;
    if (dataalign != "0"){
      string align;                                
      if (dataalign=="c") align="\"center\"";
      else if (dataalign=="r") align="\"right\"";
      else if (dataalign=="l") align="\"left\"";
      else {
        fprintf(stderr,
                "ERROR: Unrecognized table alignment argument %s for ea=X\n",
                dataalign.c_str());
        exit(1);
      }
      tbalign = " ALIGN=" + align;
    } else tbalign="";
    
    // set up vertical  alignment for particular columns

    string va;
    if (rowvalign != "0"){
      string valign;
      if (rowvalign == "t") valign= "top";
      else if (rowvalign == "m") valign= "middle";
      else if (rowvalign == "ba") valign= "baseline";
      else if (rowvalign == "bo") valign= "bottom";
      else {
        fprintf(stderr,
                "ERROR: Unrecognized table alignment argument %s for eva=X\n",
                rowvalign.c_str());
        exit(1);
      }
      va = " VALIGN =\"" + valign + "\"";
    } else va="";

    //tr_tag is keyword for data in rows

    string tr_tag= "<TR" + tbalign + va + ">";

    //declare integers to help with counting and finding position
    
    int currentc=0; // current column
    int nend = 0;
    int n1=0;
    long n = find_n(s,nend,n1);    
    
    // if there are no separators, go to the end of the string

    if (n < 0) n = s.length();

    // while n exists:

    while (n != static_cast<long>(string::npos)) {

      // ignore = 0 when pass by \n because looking for delimiters only
      // when ignore==0 do not put in a <tr>
      int ignore=1;

      // For each loop starts nend at n
      nend=n;
      
      // current column is 0, (very first loop), insert first <TR>
      if (currentc == 0){    
        currentc++;
        DT=td_tag(currentc);
        s.insert(0,tr_tag);
        s.insert(tr_tag.length(),DT);
        nend=nend+tr_tag.length()+DT.length();
        n = find_n(s,nend,n1);
        if (n==n1) currentc++;  
        else { 
          // currentc will remain one if rowquit==0
          if (rowquit>0){
            s.erase(n,1);
            n = find_n(s,nend,n1);
            currentc++;
          }
        }
      } else {
      
        // if n is separator
        if (n == n1){
          s.erase(n,tabledelim.length());
          if(currentc==(rowquit+1)&& rowquit!=0){
            s.insert(nend,"</TD></TR>\n");
            nend=nend+11;
            // set current column back to one to start new line
            currentc=1;
          }else{
            DT= td_tag(currentc);
            s.insert (nend,"</TD>");
            nend=nend+5;
            s.insert (nend,DT);
            nend=nend+DT.length();
            // add one so current column is updated
            currentc++;   
            n = find_n(s,nend,n1);
          }

        }
        //if n is newline character
        else{
          s.erase(n,1);
          // if columns == 0 means ARE searching for newlines
          // else erase and ignore insert <tr> later and
          // search for next separator

          if (rowquit==0){
            s.insert(nend,"</TD></TR>\n");
            nend=nend+11;
            // set current column back to one to start new line
            currentc=1;
          }else{
            ignore=0;
            n = find_n(s,nend,n1);
          }
        }

        // if we are at the beginning of the row then insert <TR>

        if (currentc==1&&ignore) {
          DT = td_tag(currentc); // find DT for currentc=1
          s.insert(nend,tr_tag);
          nend=nend+tr_tag.length();
          s.insert(nend,DT);
          n = find_n(s,nend,n1); // search for next separator
          currentc++;
        }
      } // end to else statement
    } // end to while loop
  } // end to if tableflag

  // if listflag is set, put list marker at beginning of every line

  if (listflag != "") {
    string marker;
    int toggle = 0;

    n = s.find('\n');
    while (n != string::npos) {
      m = s.rfind('\n',n-1);
      if (listflag == "dl" && toggle == 0) marker = "<DT>";
      else if (listflag == "dl" && toggle == 1) marker = "<DD>";
      else marker = "<LI>";
      if (m == string::npos) s.insert(0,marker);
      else s.insert(m+1,marker);
      n = s.find('\n',m+1);
      n = s.find('\n',n+1);
      if (toggle) toggle = 0;
      else toggle = 1;
    }

    listflag = "";
  }

  // if allflag is set, add markers to every line

  if (allflag != "") {
    string marker1,marker2;
    if (allflag == "p") {
      marker1 = "<P>";
      marker2 = "</P>";
    } else if (allflag == "c") {
      marker1 = "<CENTER>";
      marker2 = "</CENTER>";
    } else if (allflag == "b") {
      marker1 = "";
      marker2 = "<BR>";
    } else if (allflag == "l") {
      marker1 = "<LI>";
      marker2 = "";
    } else marker1 = marker2 = "";

    n = s.find('\n');
    while (n != string::npos) {
      m = s.rfind('\n',n-1);
      if (m == string::npos) s.insert(0,marker1);
      else s.insert(m+1,marker1);
      n = s.find('\n',m+1);
      s.insert(n,marker2);
      n = s.find('\n',n);
      n = s.find('\n',n+1);
    }

    allflag = "";
  }
}

// open input file as-is or as file.txt
// if npair = 0, don't open output file (is just initial pass thru input)
// if npair = 1, open output file as stdout
// if npair > 1, open output file with .html suffix
//   either replace .txt in input file, or append .html

void file_open(int npair, string &infile, FILE **in, FILE **out)
{
  *in = fopen(infile.c_str(),"r");
  if (*in == NULL) {
    string root = infile;
    infile = infile + ".txt";
    *in = fopen(infile.c_str(),"r");
    if (*in == NULL) {
      fprintf(stderr,"ERROR: Could not open %s or %s\n",
              root.c_str(),infile.c_str());
      exit(1);
    }
  }

  if (npair == 0) return;
  else if (npair == 1) *out = stdout;
  else {
    string outfile;
    size_t pos = infile.rfind(".txt");
    if (pos == infile.length()-4) outfile = infile.substr(0,pos) + ".html";
    else outfile = infile + ".html";
    *out = fopen(outfile.c_str(),"w");
    if (*out == NULL) {
      fprintf(stderr,"ERROR: Could not open %s\n",outfile.c_str());
      exit(1);
    }
  }
}

// for tables:
// build <TD> string (DT) based on current column

string td_tag(int currentc) {

  // eacolumn gives the alignment printout of a specific column
  string eacolumn;
  // va gives vertical alignment to a specific column
  string va;
  // DT is the complete <td> tag, with width and align
  string DT; 
  // dw is the width for tables.  It is also the <dt> tag beginning
  string dw;
  
  // set up alignment for particular columns

  for (int counter=0; counter < ncalign; counter++){
    if (ncalign != 0 && acolnum[counter] == currentc){
      string align;
      if (colalign[counter] == "l") align= "left";
      else if (colalign[counter] == "r") align= "right";
      else if (colalign[counter] == "c") align= "center";
      else {
        fprintf(stderr,
                "ERROR: Unrecognized table alignment argument %s for caM=X\n",
                colalign[counter].c_str());
        exit(1);
      }
      eacolumn= " ALIGN =\"" + align +"\"";
    }else eacolumn= "";
  }

  // set up vertical  alignment for particular columns

  for (int counter=0; counter < ncvalign; counter++){
    if (ncvalign != 0 && vacolnum[counter] == currentc){
      string valign;
      if (colvalign[counter] == "t") valign= "top";
      else if (colvalign[counter] == "m") valign= "middle";
      else if (colvalign[counter] == "ba") valign= "baseline";
      else if (colvalign[counter] == "bo") valign= "bottom";
      else {
        fprintf(stderr,
                "ERROR: Unrecognized table alignment argument %s for cvaM=X\n",
                colvalign[counter].c_str());
        exit(1);
      }
      va = " VALIGN =\"" + valign + "\"";
    } else va = " ";
  }

  // put in special width if specified
  // new code
  // if dwidth has not been set, dw is blank
  // if dwidth has been set, dw has that... unless

  if (dwidth=="0") dw = " ";       
  else dw =" WIDTH=\""+ dwidth + "\"";

  for (int counter = 0; counter < ncnum; counter++){ 
    // if it is the right column, dw = cwidth property
    if (cnum[counter] == currentc) dw= " WIDTH=\"" + cwidth[counter] + "\"";
  }
  
  // DT is set for all of this particular separator : reset next separator

  DT = "<TD" + dw + eacolumn + va + ">"; 
  
  return DT;
}

// for tables:
// find the next separator starting at nend(the end of the last .insert)
// if there is either a delim or newline 
// decide which is first
// set n = to that position
// nsep is position of the next separator. changes in here.

long find_n(string &s, int nend, int &nsep)
  // nsep is position of the next separator. changes in here.
{
  long n;
  nsep = s.find(tabledelim,nend);
  long n2 = s.find('\n',nend);
  long m = s.length() - 1;
  if (nsep >= 0 && n2 >= 0) {
    if (nsep <= n2) n = nsep;
    else n = n2;
  } else {
    if (nsep >= 0) n = nsep;
    else{
      if (n2 < m) n = n2;
      else n = string::npos;
    }
  }
      
  return n;
}
