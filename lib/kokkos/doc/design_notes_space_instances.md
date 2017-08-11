# Design Notes for Execution and Memory Space Instances

## Objective

 * Enable Kokkos interoperability with coarse-grain tasking models
 
## Requirements

 * Backwards compatable with existing Kokkos API
 * Support existing Host execution spaces (Serial, Threads, OpenMP, maybe Qthreads)
 * Support DARMA threading model (may require a new Host execution space)
 * Support Uintah threading model, i.e. indepentant worker threadpools working of of shared task queues
 
  
## Execution Space

  * Parallel work is *dispatched* on an execution space instance
  
  * Execution space instances are conceptually disjoint/independant from each other 
  

## Host Execution Space Instances

  *  A host-side *control* thread dispatches work to an instance

  * `main` is the initial control thread

  *  A host execution space instance is an organized thread pool

  *  All instances are disjoint, i.e. hardware resources are not shared between instances

  *  Exactly one control thread is associated with
     an instance and only that control thread may
     dispatch work to to that instance

  *  The control thread is a member of the instance

  *  The pool of threads associated with an instances is not mutatable during that instance existance

  *  The pool of threads associated with an instance may be masked

    -  Allows work to be dispatched to a subset of the pool

    -  Example: only one hyperthread per core of the instance

    -  A mask can be applied during the policy creation of a parallel algorithm
 
    -  Masking is portable by defining it as ceiling of fraction between [0.0, 1.0] 
       of the available resources

```
class ExecutionSpace {
public:
  using execution_space = ExecutionSpace;
  using memory_space = ...;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  using array_layout = ...;
  using size_type = ...;
  using scratch_memory_space = ...;
  
  
  class Instance
  {
    int thread_pool_size( int depth = 0 );
    ...
  };
  
  class InstanceRequest
  {
  public:
    using Control = std::function< void( Instance * )>;
    
    InstanceRequest( Control control
                   , unsigned thread_count
                   , unsigned use_numa_count = 0
                   , unsigned use_cores_per_numa = 0
                   );    
  
  };
  
  static bool in_parallel();
  
  static bool sleep();
  static bool wake();
  
  static void fence();
  
  static void print_configuration( std::ostream &, const bool detailed = false );
  
  static void initialize( unsigned thread_count = 0
                        , unsigned use_numa_count = 0
                        , unsigned use_cores_per_numa = 0
                        );
  
  // Partition the current instance into the requested instances
  // and run the given functions on the cooresponding instances
  // will block until all the partitioned instances complete and 
  // the original instance will be restored 
  //
  // Requires that the space has already been initialized
  // Requires that the request can be statisfied by the current instance
  //   i.e. the sum of number of requested threads must be less than the 
  //   max_hardware_threads
  //
  // Each control functor will accept a handle to its new default instance
  // Each instance must be independant of all other instances 
  //   i.e. no assumption on scheduling between instances
  // The user is responible for checking the return code for errors
  static int run_instances( std::vector< InstanceRequest> const& requests );
  
  static void finalize();

  static int is_initialized();
  
  static int concurrency();
  
  static int thread_pool_size( int depth = 0 );
  
  static int thread_pool_rank();
  
  static int max_hardware_threads();
  
  static int hardware_thread_id();
                        
 };

```
 



