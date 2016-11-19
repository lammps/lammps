#pragma once
#ifndef __SYMBOLIC_TASK_HPP__
#define __SYMBOLIC_TASK_HPP__

/// \file symbolic_task.hpp
/// \brief Provides tasking interface with graphviz output.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 
  
  using namespace std;

  /// \brief Graphviz color mapping for the generated tasks.
  static map<string,string> g_graphviz_color = {
    { "chol/scalar", "indianred2"},
    { "chol/trsm",   "orange2"   },
    { "chol/gemm",   "lightblue2"} };

  class SymbolicTaskQueue;

  class SymbolicTask {
  private:
    string _name;
    set<SymbolicTask*> _dep_tasks;

  public:
    // at this moment, make the queue global
    // but this should be local and work with 
    // multiple queues with separate thread teams
    typedef SymbolicTaskQueue queue;

    SymbolicTask() 
      : _name("no-name") 
    { }
    
    SymbolicTask(const SymbolicTask &b) 
      : _name(b._name)
    { }
    
    SymbolicTask(const string name) 
      : _name(name) 
    { }

    int addDependence(SymbolicTask *b) {
      if (b != NULL) 
        _dep_tasks.insert(b);
      return 0;
    }

    int clearDependence() {
      _dep_tasks.clear();
      return 0;
    }

    ostream& showMe(ostream &os) const {
      os << "    uid = " << this << " , name = " << _name << ", # of deps = " << _dep_tasks.size()  << endl;
      if (_dep_tasks.size()) {
        for (auto it=_dep_tasks.begin();it!=_dep_tasks.end();++it) 
          os << "          " << (*it) << " , name = " << (*it)->_name << endl;
      }
      return os;
    }    

    ostream& graphviz(ostream &os) const {
      os << (long)(this) 
         << " [label=\"" << _name ;
      auto it = g_graphviz_color.find(_name);
      if (it != g_graphviz_color.end())
        os << "\" ,style=filled,color=\"" << it->second << "\" "; 
      os << "];";
      for (auto it=_dep_tasks.begin();it!=_dep_tasks.end();++it) 
        os << (long)(*it) << " -> " << (long)this << ";";
      return (os << endl);
    }

  };

  static vector<SymbolicTask*> g_queue;

  class SymbolicTaskQueue {
  public:
    static SymbolicTask* push(SymbolicTask *task) {
      g_queue.push_back(task);
      return g_queue.back();
    }

    static int clear() {
      for (auto it=g_queue.begin();it!=g_queue.end();++it)
        delete (*it);
      g_queue.clear();
      return 0;
    }

    static ostream& showMe(ostream &os) {
      if (g_queue.size()) {
        os << " -- Symbolic Task Queue -- " << endl;
        for (auto it=g_queue.begin();it!=g_queue.end();++it)
          (*it)->showMe(os);
      } else {
        os << " -- Symbolic Task Queue is empty -- " << endl;
      }
      return os;
    }

    static ostream& graphviz(ostream &os, 
                             const double width = 7.5,
                             const double length = 10.0) {
      os << "digraph TaskGraph {" << endl;
      os << "size=\"" << width << "," << length << "\";" << endl;
      for (auto it=g_queue.begin();it!=g_queue.end();++it) 
        (*it)->graphviz(os);
      os << "}" << endl;
      return (os << endl);
    }
  };
  
}
#endif
