
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_fetch_add(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());
    #pragma omp atomic capture                                                    
    { tmp = *dest;  *dest += value; }                                             
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_fetch_sub(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { tmp = *dest;  *dest -= value; }                                             
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_fetch_and(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    #pragma omp atomic capture                                                    
    { tmp = *dest;  *dest &= value; }                                             
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_fetch_or(   
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { tmp = *dest;  *dest |= value; }                                             
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_fetch_xor(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { tmp = *dest;  *dest ^= value; }                                             
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_add_fetch(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { *dest += value; tmp = *dest; }                                              
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_sub_fetch(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { *dest -= value; tmp = *dest; }                                              
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_and_fetch(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { *dest &= value; tmp = *dest; }                                              
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_or_fetch(   
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { *dest |= value; tmp = *dest; }                                              
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }                                                                               
  template <typename T>                                                           
  std::enable_if_t<Impl::is_openmp_atomic_type_v<T>,T> host_atomic_xor_fetch(  
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       
    T tmp;                                                                        
    Impl::openmp_maybe_call_pre_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());    
    #pragma omp atomic capture                                                    
    { *dest ^= value; tmp = *dest; }                                              
    Impl::openmp_maybe_call_post_capture_flush(MEMORY_ORDER(), MEMORY_SCOPE());   
    return tmp;                                                                   
  }
