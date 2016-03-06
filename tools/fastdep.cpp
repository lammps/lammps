// Fast dependency generator for LAMMPS
// (c) 2016 Axel Kohlmeyer <akohlmey@gmail.com>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
// * Neither the name of the <organization> nor the
//   names of its contributors may be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//   ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
//   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
//   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <deque>
#include <set>
#include <map>
#include <iostream>

using namespace std;
const string version("1.1");

deque<string> paths; // directory search path for include files
deque<string> src;   // list of object file sources
deque<string> todo;  // queue of pending files to be processed
set<string>   incl;  // list of unique include file names
map<string, set<string> > deps; // map of direct dependencies of a file

// add to list of paths after removing trailing path separators
void add_path(char *path)
{
    int last = strlen(path) - 1;
    while ((path[last] == '/') || (path[last] == '\\')) 
        --last;
    path[++last] = '\0';
    paths.push_back(string(path));
}

// combine provided file name with one of the include file paths.
// return "", if there is no match with any of them.
string make_path(const char *file)
{
    FILE *fp;
    string full;
    deque<string>::const_iterator i;
    
    for (i=paths.begin(); i != paths.end(); ++i) {
        // build path
        full = *i;        
        full += '/';
        full += file;

        // verify existence by trying to open the file.
        // return path, if successful.
        fp = fopen(full.c_str(),"r");
        if (fp != NULL) {
            fclose(fp);
            return full;
        }
    }
        
    full = string("");
    return full;
}

// scan all provided source files for #include "..." statements
void find_includes(bool append)
{
    FILE *fp;
    const char *file;
    char *buffer,*ptr,*end;

    buffer = new char[4096];
    while (todo.size() > 0) {
        file = todo.front().c_str();
        fp = fopen(file,"r");
        if (fp == NULL) {
            perror("Cannot read source");
            fprintf(stderr,"For file: %s\n",file);
            exit(EXIT_FAILURE);
        }
        
        // read file line by line and look for #include "..."
        while (!feof(fp) && !ferror(fp)) {
            if (fgets(buffer,4096,fp) == NULL) continue;
            ptr = buffer;
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (*ptr != '#') continue;
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (*++ptr != 'i') continue;
            if (*++ptr != 'n') continue;
            if (*++ptr != 'c') continue;
            if (*++ptr != 'l') continue;
            if (*++ptr != 'u') continue;
            if (*++ptr != 'd') continue;
            if (*++ptr != 'e') continue;
            ++ptr;
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (*ptr != '"') continue;
            ++ptr;
            end = ptr;
            while (*end != '"') {
                if (*end == '\0') {
                    fprintf(stderr,"Unmatched '\"': %s\n",buffer);
                    exit(EXIT_FAILURE);
                }
                ++end;
            }
            *end = '\0';

            // record include file with path in set
            string full = make_path(ptr);
            // skip, if not found or readable.
            if (full == "") continue;

            // if this is a yet unknown include, add to the
            // todo list, if append is enabled
            if (append && (incl.count(full) == 0)) {
                todo.push_back(full);
            }

            incl.insert(full);
            deps[todo.front()].insert(full);
        }
        fclose(fp);
        todo.pop_front();
    }
    delete[] buffer;
}

// recurse through dependencies
void add_depend(const set<string> &mydeps)
{
    set<string>::const_iterator d;
    for (d = mydeps.begin(); d != mydeps.end(); ++d) {
        if (incl.count(*d) == 0) {
            incl.insert(*d);
            add_depend(deps[*d]);
        }
    }
}

void do_depend()
{
    string target;
    size_t pos;
    
    deque<string>::const_iterator s;
    for (s = src.begin(); s != src.end(); ++s) {
        target = *s;
        pos = target.rfind("/");
        if (pos != string::npos) target.erase(0,pos+1);
        pos = target.rfind(".");
        if (pos != string::npos) {
            target.replace(target.begin()+pos,target.end(),".o");
        }
        
        cout << target << " : " << *s;
        incl.clear();
        add_depend(deps[*s]);

        set<string>::const_iterator d=incl.begin();
        for (; d != incl.end(); ++d) {
            cout << " " << *d;
        }
        cout << endl;
    }
    cout << endl;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        cout << "FastDep v" << version << " for LAMMPS\nUsage: "
             << argv[0] << " [-I <path> ...] -- <src1> [<src2> ...]\n";
        return 1;
    }
    
    cout << "# fastdep v" << version << endl;
    paths.push_back(string("."));
    paths.push_back(string(".."));
    
    while (++argv, --argc > 0) {

        if (strncmp(*argv, "-I", 2) == 0) {

            if ((*argv)[2] != '\0') {
                add_path(*argv+2);
            } else {
                ++argv;
                --argc;

                if (argc > 0) {
                    if (strcmp(*argv,"--") == 0) {
                        cout << "Error: -I flag without path\n";
                        return 1;
                    } else {
                        add_path(*argv);
                    }
                } else {
                    cout << "Error: -I flag without path\n";
                    return 1;
                }
            }
        } else if (strcmp(*argv,"--") == 0) {
            break;
        } // ignore all unrecognized arguments before '--'.
    }

    while (++argv, --argc > 0) {
        src.push_back(*argv);
    }

    deque<string>::const_iterator p=paths.begin();
    cout << "# Search path: " << *p;
    for (++p; p != paths.end(); ++p) {
        cout << ":" << *p;
    }
    cout << endl;

    // process the main source files first to extract non-system includes
    todo = src;
    find_includes(false);

    // now fill the todo list with currently known included files
    // and start the search again.
    set<string>::const_iterator i;
    todo.clear();               // should not be needed.
    for (i = incl.begin(); i != incl.end(); ++i) {
        todo.push_back(*i);
    }
    find_includes(true);
    cout << "# " << src.size() << " sources\n";
    cout << "# " << incl.size() << " includes\n";
    cout << "# " << deps.size() << " files with dependencies\n";

    do_depend();
}

//
// Local Variables:
// compile-command: "g++ -o fastdep.exe -Wall -g -O fastdep.cpp"
// c-basic-offset: 4
// End:
//
