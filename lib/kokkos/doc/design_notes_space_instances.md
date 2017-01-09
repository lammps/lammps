# Design Notes for Execution and Memory Space Instances


## Execution Spaces

  *  Work is *dispatched* to an execution space instance



## Host Associated Execution Space Instances

Vocabulary and examples assuming C++11 Threads Support Library

  *  A host-side *control* thread dispatches work to an instance

  * `this_thread` is the control thread

  * `main` is the initial control thread

  *  An execution space instance is a pool of threads

  *  All instances are disjoint thread pools

  *  Exactly one control thread is associated with
     an instance and only that control thread may
     dispatch work to to that instance

  *  A control thread may be a member of an instance,
     if so then it is also the control thread associated
     with that instance

  *  The pool of threads associated with an instances is not mutatable

  *  The pool of threads associated with an instance may be masked

    -  Allows work to be dispatched to a subset of the pool

    -  Example: only one hyperthread per core of the instance

    -  When a mask is applied to an instance that mask
       remains until cleared or another mask is applied

    -  Masking is portable by defining it as using a fraction
       of the available resources (threads)

  *  Instances are shared (referenced counted) objects,
     just like `Kokkos::View`

```
struct StdThread {
  void mask( float fraction );
  void unmask() { mask( 1.0 ); }
};
```



### Requesting an Execution Space Instance

  *  `Space::request(` *who* `,` *what* `,` *control-opt* `)`

  *  *who* is an identifier for subsquent queries regarding
    who requested each instance

  *  *what* is the number of threads and how they should be placed

    -  Placement within locality-topology hierarchy; e.g., HWLOC

    -  Compact within a level of hierarchy, or striped across that level;
       e.g., socket or NUMA region

    -  Granularity of request is core

  *  *control-opt*  optionally specifies whether the instance
     has a new control thread

    -  *control-opt* includes a control function / closure

    -  The new control thread is a member of the instance

    -  The control function is called by the new control thread
       and is passed a `const` instance

    -  The instance is **not** returned to the creating control thread

  *  `std::thread` that is not a member of an instance is
     *hard blocked* on a `std::mutex`

    -  One global mutex or one mutex per thread?

  *  `std::thread` that is a member of an instance is
     *spinning* waiting for work, or are working

```
struct StdThread {

  struct Resource ;

  static StdThread request(); // default

  static StdThread request( const std::string & , const Resource & );

  // If the instance can be reserved then
  // allocate a copy of ControlClosure and invoke
  //   ControlClosure::operator()( const StdThread intance ) const
  template< class ControlClosure >
  static bool request( const std::string & , const Resource &
                     , const ControlClosure & );
};
```

### Relinquishing an Execution Space Instance

  *  De-referencing the last reference-counted instance
     relinquishes the pool of threads

  *  If a control thread was created for the instance then
     it is relinquished when that control thread returns
     from the control function

    -  Requires the reference count to be zero, an error if not

  *  No *forced* relinquish



## CUDA Associated Execution Space Instances

  *  Only a signle CUDA architecture

  *  An instance is a device + stream

  *  A stream is exclusive to an instance

  *  Only a host-side control thread can dispatch work to an instance

  *  Finite number of streams per device

  *  ISSUE:  How to use CUDA `const` memory with multiple streams?

  *  Masking can be mapped to restricting the number of CUDA blocks
     to the fraction of available resources; e.g., maximum resident blocks


### Requesting an Execution Space Instance

  *  `Space::request(` *who* `,` *what* `)`

  *  *who* is an identifier for subsquent queries regarding
    who requested each instance

  *  *what* is which device, the stream is a requested/relinquished resource


```
struct Cuda {

  struct Resource ;

  static Cuda request();

  static Cuda request( const std::string & , const Resource & );
};
```


