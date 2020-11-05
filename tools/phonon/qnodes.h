#ifndef QNODES_H
#define QNODES_H

#include <vector>
#include <string>

class QNodes {
public:
   QNodes();
   ~QNodes();
   
   std::vector<double> nodes;
   std::vector<std::string> ndstr;
   std::vector<double *> qs, qe;
   std::vector<int> nqbin;
};
#endif
