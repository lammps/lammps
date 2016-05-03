#pragma once
#ifndef __TASK_FACTORY_HPP__
#define __TASK_FACTORY_HPP__

/// \file task_factory.hpp
/// \brief A wrapper for task policy and future with a provided space type.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  using namespace std;

  /// \class TaskFactory
  /// \brief Minimal interface to Kokkos tasking.
  ///
  /// TaskFactory is attached to blocks as a template argument in order to 
  /// create and manage tasking future objects. Note that policy (shared 
  /// pointer to the task generator) is not a member object in this class.
  /// This class includes minimum interface for tasking with type decralation 
  /// of the task policy and template alias of future so that future objects 
  /// generated in this class will match to their policy and its execution space. 
  ///
  template<typename PolicyType,        
           typename FutureType>
  class TaskFactory {
  private:
    static constexpr int _max_task_dependence = 10 ;

  public:
    typedef PolicyType policy_type;
    typedef FutureType future_type;
    
    template<typename TaskFunctorType>
    static KOKKOS_INLINE_FUNCTION
    future_type create(policy_type &policy, const TaskFunctorType &func) {

      future_type f ;
      // while ( f.is_null() ) {
        f = policy.task_create_team(func, _max_task_dependence);
      // }
      if ( f.is_null() ) Kokkos::abort("task_create_team FAILED, out of memory");
      return f ;
    }
    
    static KOKKOS_INLINE_FUNCTION
    void spawn(policy_type &policy, const future_type &obj, bool priority = false ) {
      policy.spawn(obj,priority);
    }
    
    static KOKKOS_INLINE_FUNCTION
    void addDependence(policy_type &policy, 
                       const future_type &after, const future_type &before) {
      policy.add_dependence(after, before);
    }

    template<typename TaskFunctorType>
    static  KOKKOS_INLINE_FUNCTION
    void addDependence(policy_type &policy, 
                       TaskFunctorType *after, const future_type &before) {
      policy.add_dependence(after, before);
    }

    template<typename TaskFunctorType>
    static  KOKKOS_INLINE_FUNCTION
    void clearDependence(policy_type &policy, TaskFunctorType *func) {
      policy.clear_dependence(func);
    }

    template<typename TaskFunctorType>
    static KOKKOS_INLINE_FUNCTION
    void respawn(policy_type &policy, TaskFunctorType *func) {
      policy.respawn(func);
    }
  };
}

#endif
