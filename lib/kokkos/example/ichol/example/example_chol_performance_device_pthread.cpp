#include <Kokkos_Core.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "example_chol_performance_device.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  string file_input = "test.mtx";                                                                             
  int nthreads = 1;                                                                                           
  int max_task_dependence = 3;                                                                                
  int max_concurrency = 1024;                                                                                 
  int team_size = 1;                                                                                          
  int fill_level = 0;
  int treecut = 0;
  int prunecut = 0;
  int seed = 0;
  int league_size = 1;                                                                                        
  bool verbose = false;                                                                                       
  for (int i=0;i<argc;++i) {                                                                                  
    if ((strcmp(argv[i],"--file-input")         ==0)) { file_input          = argv[++i];       continue;}     
    if ((strcmp(argv[i],"--nthreads")           ==0)) { nthreads            = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--max-task-dependence")==0)) { max_task_dependence = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--max-concurrency")    ==0)) { max_concurrency     = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--team-size")          ==0)) { team_size           = atoi(argv[++i]); continue;}     

    if ((strcmp(argv[i],"--fill-level")         ==0)) { fill_level          = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--league-size")        ==0)) { league_size         = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--treecut")            ==0)) { treecut             = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--prunecut")           ==0)) { prunecut            = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--seed")               ==0)) { seed                = atoi(argv[++i]); continue;}     
    if ((strcmp(argv[i],"--enable-verbose")     ==0)) { verbose             = true;            continue;}     
  }                                                                                                           

  int r_val = 0;
  {
    exec_space::initialize(nthreads);
    exec_space::print_configuration(cout, true);

    r_val = exampleCholPerformanceDevice
      <value_type,ordinal_type,size_type,exec_space>
      (file_input,
       treecut,
       prunecut,
       seed,
       nthreads,
       max_task_dependence, max_concurrency, team_size,
       fill_level, league_size,
       (nthreads != 1), // skip_serial
       verbose);

    exec_space::finalize();
  }

  return r_val;
}
