#include "qnodes.h"

/* ----------------------------------------------------------------------------
 * Class QNodes stores the high symmetry k-path nodes for a given lattice.
 * The constructor and the deconstructor simply empties the data.
 * ---------------------------------------------------------------------------- */
QNodes::QNodes()
{
   nodes.clear();
   ndstr.clear();
   qs.clear();
   qe.clear();
   nqbin.clear();
}

/* ----------------------------------------------------------------------------
 * The constructor and the deconstructor simply empties the data.
 * ---------------------------------------------------------------------------- */
QNodes::~QNodes()
{
   nodes.clear();
   ndstr.clear();
   qs.clear();
   qe.clear();
   nqbin.clear();
}
