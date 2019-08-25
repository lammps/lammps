/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace Test {
struct SomeTag {};

template< class ExecutionSpace >
class TestRangePolicyConstruction {
public:
  TestRangePolicyConstruction() {
    test_compile_time_parameters();
    test_runtime_parameters();
  }

private:
  void test_compile_time_parameters() {
    {
      Kokkos::Impl::expand_variadic();
      Kokkos::Impl::expand_variadic( 1, 2, 3 );
    }

    {
      typedef Kokkos::RangePolicy<> policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Static>    >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< ExecutionSpace > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Static>    >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::IndexType<long>, ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::Schedule<Kokkos::Dynamic>, ExecutionSpace, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< SomeTag, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, ExecutionSpace > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::IndexType<long>, Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::RangePolicy< SomeTag, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }
  }
  void test_runtime_parameters() {
    {
      typedef Kokkos::RangePolicy<> policy_t;
      policy_t p(5,15);
      ASSERT_TRUE( (p.begin() == 5) );
      ASSERT_TRUE( (p.end() == 15) );
    }
    {
      typedef Kokkos::RangePolicy<> policy_t;
      policy_t p(Kokkos::DefaultExecutionSpace(),5,15);
      ASSERT_TRUE( (p.begin() == 5) );
      ASSERT_TRUE( (p.end() == 15) );
    }
    {
      typedef Kokkos::RangePolicy<> policy_t;
      policy_t p(5,15,Kokkos::ChunkSize(10));
      ASSERT_TRUE( (p.begin() == 5) );
      ASSERT_TRUE( (p.end() == 15) );
      ASSERT_TRUE( (p.chunk_size() == 10) );
    }
    {
      typedef Kokkos::RangePolicy<> policy_t;
      policy_t p(Kokkos::DefaultExecutionSpace(),5,15,Kokkos::ChunkSize(10));
      ASSERT_TRUE( (p.begin() == 5) );
      ASSERT_TRUE( (p.end() == 15) );
      ASSERT_TRUE( (p.chunk_size() == 10) );
    }
  }
};

template< class ExecutionSpace >
class TestTeamPolicyConstruction {
public:
  TestTeamPolicyConstruction() {
    test_compile_time_parameters();
    test_run_time_parameters();
  }

private:
  void test_compile_time_parameters() {
    {
      typedef Kokkos::TeamPolicy<> policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Static>    >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< ExecutionSpace > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Static>    >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::IndexType<long>, ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::Schedule<Kokkos::Dynamic>, ExecutionSpace, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< SomeTag, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, ExecutionSpace > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, ExecutionSpace                      >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace        >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      typename execution_space::size_type >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::IndexType<long>, Kokkos::Schedule<Kokkos::Dynamic> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        void                                >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, SomeTag > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }

    {
      typedef Kokkos::TeamPolicy< SomeTag, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > policy_t;
      typedef typename policy_t::execution_space  execution_space;
      typedef typename policy_t::index_type       index_type;
      typedef typename policy_t::schedule_type    schedule_type;
      typedef typename policy_t::work_tag         work_tag;

      ASSERT_TRUE( ( std::is_same< execution_space, Kokkos::DefaultExecutionSpace       >::value ) );
      ASSERT_TRUE( ( std::is_same< index_type,      long                                >::value ) );
      ASSERT_TRUE( ( std::is_same< schedule_type,   Kokkos::Schedule<Kokkos::Dynamic>   >::value ) );
      ASSERT_TRUE( ( std::is_same< work_tag,        SomeTag                             >::value ) );
    }
  }


  template< class policy_t >
  void test_run_time_parameters_type() {
    int league_size = 131;
    int team_size = 4 < policy_t::execution_space::concurrency() ? 4 : policy_t::execution_space::concurrency();
#ifdef KOKKOS_ENABLE_HPX
    team_size = 1;
#endif
    int chunk_size = 4;
    int per_team_scratch = 1024;
    int per_thread_scratch = 16;
    int scratch_size = per_team_scratch + per_thread_scratch * team_size;
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    int vector_length = 4;
#endif

    policy_t p1( league_size, team_size );
    ASSERT_EQ  ( p1.league_size(),     league_size                    );
    ASSERT_EQ  ( p1.team_size(),       team_size                      );
    ASSERT_TRUE( p1.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p1.scratch_size( 0 ), 0                              );

    policy_t p2 = p1.set_chunk_size( chunk_size );
    ASSERT_EQ  ( p1.league_size(),     league_size                    );
    ASSERT_EQ  ( p1.team_size(),       team_size                      );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_TRUE( p1.chunk_size()  > 0                                 );
#else
    ASSERT_EQ  ( p1.chunk_size(),      chunk_size                     );
#endif
    ASSERT_EQ  ( p1.scratch_size( 0 ), 0                              );

    ASSERT_EQ  ( p2.league_size(),     league_size                    );
    ASSERT_EQ  ( p2.team_size(),       team_size                      );
    ASSERT_EQ  ( p2.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p2.scratch_size( 0 ), 0                              );

    policy_t p3 = p2.set_scratch_size( 0, Kokkos::PerTeam( per_team_scratch ) );
    ASSERT_EQ  ( p2.league_size(),     league_size                    );
    ASSERT_EQ  ( p2.team_size(),       team_size                      );
    ASSERT_EQ  ( p2.chunk_size(),      chunk_size                     );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_EQ  ( p2.scratch_size( 0 ), 0                              );
#else
    ASSERT_EQ  ( p2.scratch_size( 0 ), per_team_scratch               );
#endif
    ASSERT_EQ  ( p3.league_size(),     league_size                    );
    ASSERT_EQ  ( p3.team_size(),       team_size                      );
    ASSERT_EQ  ( p3.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p3.scratch_size( 0 ), per_team_scratch               );

    policy_t p4 = p2.set_scratch_size( 0, Kokkos::PerThread( per_thread_scratch ) );
    ASSERT_EQ  ( p2.league_size(),     league_size                    );
    ASSERT_EQ  ( p2.team_size(),       team_size                      );
    ASSERT_EQ  ( p2.chunk_size(),      chunk_size                     );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_EQ  ( p2.scratch_size( 0 ), 0                              );
#else
    ASSERT_EQ  ( p2.scratch_size( 0 ), scratch_size );
#endif
    ASSERT_EQ  ( p4.league_size(),     league_size                    );
    ASSERT_EQ  ( p4.team_size(),       team_size                      );
    ASSERT_EQ  ( p4.chunk_size(),      chunk_size                     );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_EQ  ( p4.scratch_size( 0 ), per_thread_scratch * team_size );
#else
    ASSERT_EQ  ( p4.scratch_size( 0 ), scratch_size );
#endif

    policy_t p5 = p2.set_scratch_size( 0, Kokkos::PerThread( per_thread_scratch ), Kokkos::PerTeam( per_team_scratch ) );
    ASSERT_EQ  ( p2.league_size(),     league_size                    );
    ASSERT_EQ  ( p2.team_size(),       team_size                      );
    ASSERT_EQ  ( p2.chunk_size(),      chunk_size                     );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_EQ  ( p2.scratch_size( 0 ), 0                              );
#else
    ASSERT_EQ  ( p2.scratch_size( 0 ), scratch_size                   );
#endif
    ASSERT_EQ  ( p5.league_size(),     league_size                    );
    ASSERT_EQ  ( p5.team_size(),       team_size                      );
    ASSERT_EQ  ( p5.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p5.scratch_size( 0 ), scratch_size                   );

    policy_t p6 = p2.set_scratch_size( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) );
    ASSERT_EQ  ( p2.league_size(),     league_size                    );
    ASSERT_EQ  ( p2.team_size(),       team_size                      );
    ASSERT_EQ  ( p2.chunk_size(),      chunk_size                     );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_EQ  ( p2.scratch_size( 0 ), 0                              );
#else
    ASSERT_EQ  ( p2.scratch_size( 0 ), scratch_size                   );
#endif
    ASSERT_EQ  ( p6.league_size(),     league_size                    );
    ASSERT_EQ  ( p6.team_size(),       team_size                      );
    ASSERT_EQ  ( p6.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p6.scratch_size( 0 ), scratch_size                   );

    policy_t p7 = p3.set_scratch_size( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) );
    ASSERT_EQ  ( p3.league_size(),     league_size                    );
    ASSERT_EQ  ( p3.team_size(),       team_size                      );
    ASSERT_EQ  ( p3.chunk_size(),      chunk_size                     );
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    ASSERT_EQ  ( p3.scratch_size( 0 ), per_team_scratch               );
#else
    ASSERT_EQ  ( p3.scratch_size( 0 ), scratch_size                   );
#endif
    ASSERT_EQ  ( p7.league_size(),     league_size                    );
    ASSERT_EQ  ( p7.team_size(),       team_size                      );
    ASSERT_EQ  ( p7.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p7.scratch_size( 0 ), scratch_size                   );

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    policy_t p8(league_size, team_size, Kokkos::ChunkSize(chunk_size) );
    ASSERT_EQ  ( p8.league_size(),     league_size                    );
    ASSERT_EQ  ( p8.team_size(),       team_size                      );
    ASSERT_EQ  ( p8.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p8.scratch_size( 0 ), 0                              );

    policy_t p10( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p10.league_size(),     league_size                    );
    ASSERT_EQ  ( p10.team_size(),       team_size                      );
    ASSERT_TRUE( p10.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p10.scratch_size( 0 ), per_team_scratch               );

    policy_t p11( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p11.league_size(),     league_size                    );
    ASSERT_EQ  ( p11.team_size(),       team_size                      );
    ASSERT_TRUE( p11.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p11.scratch_size( 0 ), per_thread_scratch * team_size );

    policy_t p12( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ), Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p12.league_size(),     league_size                    );
    ASSERT_EQ  ( p12.team_size(),       team_size                      );
    ASSERT_TRUE( p12.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p12.scratch_size( 0 ), scratch_size                   );

    policy_t p13( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p13.league_size(),     league_size                    );
    ASSERT_EQ  ( p13.team_size(),       team_size                      );
    ASSERT_TRUE( p13.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p13.scratch_size( 0 ), scratch_size                   );

    policy_t p14( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p14.league_size(),     league_size                    );
    ASSERT_EQ  ( p14.team_size(),       team_size                      );
    ASSERT_TRUE( p14.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p14.scratch_size( 0 ), scratch_size                   );

    policy_t p15( league_size, team_size, Kokkos::ChunkSize(chunk_size), Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p15.league_size(),     league_size                    );
    ASSERT_EQ  ( p15.team_size(),       team_size                      );
    ASSERT_TRUE( p15.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p15.scratch_size( 0 ), per_team_scratch               );

    policy_t p16( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ) ), Kokkos::ChunkSize(chunk_size) );
    ASSERT_EQ  ( p16.league_size(),     league_size                    );
    ASSERT_EQ  ( p16.team_size(),       team_size                      );
    ASSERT_EQ  ( p16.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p16.scratch_size( 0 ), per_thread_scratch * team_size );

    policy_t p17( league_size, team_size, Kokkos::ChunkSize(chunk_size), Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ), Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p17.league_size(),     league_size                    );
    ASSERT_EQ  ( p17.team_size(),       team_size                      );
    ASSERT_EQ  ( p17.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p17.scratch_size( 0 ), scratch_size                   );

    policy_t p18( league_size, team_size, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ), Kokkos::ChunkSize(chunk_size) );
    ASSERT_EQ  ( p18.league_size(),     league_size                    );
    ASSERT_EQ  ( p18.team_size(),       team_size                      );
    ASSERT_EQ  ( p18.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p18.scratch_size( 0 ), scratch_size                   );

    policy_t p19( league_size, team_size, Kokkos::ChunkSize(chunk_size), Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p19.league_size(),     league_size                    );
    ASSERT_EQ  ( p19.team_size(),       team_size                      );
    ASSERT_EQ  ( p19.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p19.scratch_size( 0 ), scratch_size                   );

    policy_t p20( league_size, team_size, vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p20.league_size(),     league_size                    );
    ASSERT_EQ  ( p20.team_size(),       team_size                      );
    ASSERT_TRUE( p20.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p20.scratch_size( 0 ), per_team_scratch               );

    policy_t p21( league_size, team_size, vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p21.league_size(),     league_size                    );
    ASSERT_EQ  ( p21.team_size(),       team_size                      );
    ASSERT_TRUE( p21.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p21.scratch_size( 0 ), per_thread_scratch * team_size );

    policy_t p22( league_size, team_size, vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ), Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p22.league_size(),     league_size                    );
    ASSERT_EQ  ( p22.team_size(),       team_size                      );
    ASSERT_TRUE( p22.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p22.scratch_size( 0 ), scratch_size                   );

    policy_t p23( league_size, team_size, (size_t) vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p23.league_size(),     league_size                    );
    ASSERT_EQ  ( p23.team_size(),       team_size                      );
    ASSERT_TRUE( p23.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p23.scratch_size( 0 ), scratch_size                   );

    policy_t p24( league_size, team_size, (size_t) vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p24.league_size(),     league_size                    );
    ASSERT_EQ  ( p24.team_size(),       team_size                      );
    ASSERT_TRUE( p24.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p24.scratch_size( 0 ), scratch_size                   );

    policy_t p25( league_size, team_size, vector_length, Kokkos::ChunkSize(chunk_size), Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p25.league_size(),     league_size                    );
    ASSERT_EQ  ( p25.team_size(),       team_size                      );
    ASSERT_TRUE( p25.chunk_size()  > 0                                 );
    ASSERT_EQ  ( p25.scratch_size( 0 ), per_team_scratch               );

    policy_t p26( league_size, team_size, vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ) ), Kokkos::ChunkSize(chunk_size) );
    ASSERT_EQ  ( p26.league_size(),     league_size                    );
    ASSERT_EQ  ( p26.team_size(),       team_size                      );
    ASSERT_EQ  ( p26.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p26.scratch_size( 0 ), per_thread_scratch * team_size );

    policy_t p27( league_size, team_size, vector_length, Kokkos::ChunkSize(chunk_size), Kokkos::ScratchRequest( 0, Kokkos::PerThread( per_thread_scratch ), Kokkos::PerTeam( per_team_scratch ) ) );
    ASSERT_EQ  ( p27.league_size(),     league_size                    );
    ASSERT_EQ  ( p27.team_size(),       team_size                      );
    ASSERT_EQ  ( p27.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p27.scratch_size( 0 ), scratch_size                   );

    policy_t p28( league_size, team_size, (size_t) vector_length, Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ), Kokkos::ChunkSize(chunk_size) );
    ASSERT_EQ  ( p28.league_size(),     league_size                    );
    ASSERT_EQ  ( p28.team_size(),       team_size                      );
    ASSERT_EQ  ( p28.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p28.scratch_size( 0 ), scratch_size                   );

    policy_t p29( league_size, team_size, (size_t) vector_length, Kokkos::ChunkSize(chunk_size), Kokkos::ScratchRequest( 0, Kokkos::PerTeam( per_team_scratch ), Kokkos::PerThread( per_thread_scratch ) ) );
    ASSERT_EQ  ( p29.league_size(),     league_size                    );
    ASSERT_EQ  ( p29.team_size(),       team_size                      );
    ASSERT_EQ  ( p29.chunk_size(),      chunk_size                     );
    ASSERT_EQ  ( p29.scratch_size( 0 ), scratch_size                   );
#endif

  }

  void test_run_time_parameters() {
    test_run_time_parameters_type< Kokkos::TeamPolicy<ExecutionSpace> >();
    test_run_time_parameters_type< Kokkos::TeamPolicy<ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long> > >();
    test_run_time_parameters_type< Kokkos::TeamPolicy<Kokkos::IndexType<long>, ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic> > >();
    test_run_time_parameters_type< Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<long>, ExecutionSpace, SomeTag > >();
  }
};

TEST_F( TEST_CATEGORY, policy_construction )
{
   TestRangePolicyConstruction< TEST_EXECSPACE >();
   TestTeamPolicyConstruction< TEST_EXECSPACE >();
}

}
