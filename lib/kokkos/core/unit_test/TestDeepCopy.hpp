#include<Kokkos_Core.hpp>

namespace Test {

namespace Impl {
template<class MemorySpaceA, class MemorySpaceB>
struct TestDeepCopy {

  typedef Kokkos::View<double*, Kokkos::LayoutRight, MemorySpaceA> a_base_t;
  typedef Kokkos::View<double*, Kokkos::LayoutRight, MemorySpaceB> b_base_t;
  typedef Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_t;
  typedef Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceB> b_char_t;

  typedef Kokkos::RangePolicy<typename MemorySpaceA::execution_space> policyA_t;
  typedef Kokkos::RangePolicy<typename MemorySpaceB::execution_space> policyB_t;

  static void reset_a_copy_and_b(Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_copy, Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceB> b_char) {
    const int N = b_char.extent(0);
    Kokkos::parallel_for("TestDeepCopy: FillA_copy",policyA_t(0,N), KOKKOS_LAMBDA (const int& i) {
      a_char_copy(i) = char(0);
    });
    Kokkos::parallel_for("TestDeepCopy: FillB",policyB_t(0,N), KOKKOS_LAMBDA (const int& i) {
      b_char(i) = char(0);
    });
  }

  static int compare_equal(Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_copy, Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char) {
    const int N = a_char.extent(0);
    int errors;
    Kokkos::parallel_reduce("TestDeepCopy: FillA_copy",policyA_t(0,N), KOKKOS_LAMBDA (const int& i, int& lsum) {
      if(a_char_copy(i) != a_char(i)) lsum++;
    },errors);
    return errors;
  }

  static void run_test(int num_bytes) {
    a_base_t a_base("test_space_to_space",(num_bytes+128)/8);
    a_base_t a_base_copy("test_space_to_space",(num_bytes+128)/8);
    Kokkos::View<double*, Kokkos::LayoutRight, MemorySpaceB> b_base("test_space_to_space",(num_bytes+128)/8);
    
    Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char((char*) a_base.data(),a_base.extent(0)*8);
    Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_copy((char*) a_base_copy.data(),a_base.extent(0)*8);
    Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceB> b_char((char*) b_base.data(),b_base.extent(0)*8);

    Kokkos::parallel_for("TestDeepCopy: FillA",policyA_t(0,a_char.extent(0)), KOKKOS_LAMBDA (const int& i) {
      a_char(i) = static_cast<char>(i%97)+1;
    });

    reset_a_copy_and_b(a_char_copy, b_char);

    {
      int check = compare_equal(a_char_copy,a_char);
      ASSERT_EQ( check, a_char.extent(0) );
    }

    // (a.data()%8, (a.data()+a.extent(0))%8, b.data()%8, (b.data()+b.extent(0))%8
    // (0,0,0,0) 
    {
      int a_begin = 0;
      int a_end   = 0;
      int b_begin = 0;
      int b_end   = 0;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }

    {
      int a_begin = 0;
      int a_end   = 5;
      int b_begin = 0;
      int b_end   = 5;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }
    
    {
      int a_begin = 3;
      int a_end   = 0;
      int b_begin = 3;
      int b_end   = 0;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }

    {
      int a_begin = 3;
      int a_end   = 6;
      int b_begin = 3;
      int b_end   = 6;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }

    {
      int a_begin = 5;
      int a_end   = 4;
      int b_begin = 3;
      int b_end   = 6;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }

    {
      int a_begin = 0;
      int a_end   = 8;
      int b_begin = 2;
      int b_end   = 6;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }
     
    {
      int a_begin = 2;
      int a_end   = 6;
      int b_begin = 0;
      int b_end   = 8;
      auto a = Kokkos::subview(a_char,std::pair<int,int>(a_begin,a_char.extent(0)-a_end));
      auto b = Kokkos::subview(b_char,std::pair<int,int>(b_begin,b_char.extent(0)-b_end));
      auto a_copy = Kokkos::subview(a_char_copy,std::pair<int,int>(a_begin,a_char_copy.extent(0)-a_end));
      Kokkos::deep_copy(b,a);
      Kokkos::deep_copy(a_copy,b);
      int check = compare_equal(a_copy,a);
      ASSERT_EQ( check, 0 );
    }

  }
};
}

TEST_F( TEST_CATEGORY, deep_copy_alignment )
{
  { Impl::TestDeepCopy< TEST_EXECSPACE::memory_space , TEST_EXECSPACE::memory_space >::run_test( 100000 ); }
  { Impl::TestDeepCopy< Kokkos::HostSpace , TEST_EXECSPACE::memory_space >::run_test( 100000 ); }
  { Impl::TestDeepCopy< TEST_EXECSPACE::memory_space , Kokkos::HostSpace >::run_test( 100000 ); }
}

}
