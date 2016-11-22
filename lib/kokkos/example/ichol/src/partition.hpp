
#ifndef __PARTITION_HPP__
#define __PARTITION_HPP__

/// \file partition.hpp
/// \brief Matrix partitioning utilities.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  using namespace std;

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Part_2x2(const MatView A, MatView &ATL, MatView &ATR, 
           /**************/ MatView &ABL, MatView &ABR,
           const typename MatView::ordinal_type bm, 
           const typename MatView::ordinal_type bn,
           const int quadrant) {
    typename MatView::ordinal_type bmm, bnn;

    switch (quadrant) {
    case Partition::TopLeft:
      bmm = min(bm, A.NumRows());
      bnn = min(bn, A.NumCols());                
      
      ATL.setView(A.BaseObject(),
                  A.OffsetRows(), bmm,
                  A.OffsetCols(), bnn);
      break;
    case Partition::TopRight:
    case Partition::BottomLeft:
      Kokkos::abort("Tacho::Part_2x2 Not yet implemented");
      break;
    case Partition::BottomRight:
      bmm = A.NumRows() - min(bm, A.NumRows());
      bnn = A.NumCols() - min(bn, A.NumCols());                
      
      ATL.setView(A.BaseObject(),
                  A.OffsetRows(), bmm,
                  A.OffsetCols(), bnn);
      break;
    default:
      Kokkos::abort("Tacho::Part_2x2 Invalid Input");
      break;
    }
    
    ATR.setView(A.BaseObject(),
                A.OffsetRows(),                 ATL.NumRows(),
                A.OffsetCols() + ATL.NumCols(), A.NumCols() - ATL.NumCols());
    
    ABL.setView(A.BaseObject(),
                A.OffsetRows() + ATL.NumRows(), A.NumRows() - ATL.NumRows(),
                A.OffsetCols(),                 ATL.NumCols());
    
    ABR.setView(A.BaseObject(),
                A.OffsetRows() + ATL.NumRows(), A.NumRows() - ATL.NumRows(),
                A.OffsetCols() + ATL.NumCols(), A.NumCols() - ATL.NumCols());
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Part_1x2(const MatView A, MatView &AL, MatView &AR, 
           const typename MatView::ordinal_type bn,
           const int side) {
    typename MatView::ordinal_type bmm, bnn;

    switch (side) {
    case Partition::Left:
      bmm = A.NumRows();
      bnn = min(bn, A.NumCols());
      
      AL.setView(A.BaseObject(),
                 A.OffsetRows(), bmm,
                 A.OffsetCols(), bnn);
      break;
    case Partition::Right:
      bmm = A.NumRows();
      bnn = A.NumCols() - min(bn, A.NumCols());

      AL.setView(A.BaseObject(),
                 A.OffsetRows(), bmm,
                 A.OffsetCols(), bnn);
      break;
    default:
      Kokkos::abort("Tacho::Part_1x2 Invalid Input");
      break;
    }

    AR.setView(A.BaseObject(),
               A.OffsetRows(),                A.NumRows(),
               A.OffsetCols() + AL.NumCols(), A.NumCols() - AL.NumCols());
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Part_2x1(const MatView A, MatView &AT, 
           /*************/  MatView &AB, 
           const typename MatView::ordinal_type bm,
           const int side) {
    typename MatView::ordinal_type bmm, bnn;
    
    switch (side) {
    case Partition::Top:
      bmm = min(bm, A.NumRows());
      bnn = A.NumCols();
      
      AT.setView(A.BaseObject(),
                 A.OffsetRows(), bmm,
                 A.OffsetCols(), bnn);
      break;
    case Partition::Bottom:
      bmm = A.NumRows() - min(bm, A.NumRows());
      bnn = A.NumCols();

      AT.setView(A.BaseObject(),
                 A.OffsetRows(), bmm,
                 A.OffsetCols(), bnn);
      break;
    default:
      Kokkos::abort("Tacho::Part_2x1 Invalid Input");
      break;
    }
    
    AB.setView(A.BaseObject(),
               A.OffsetRows() + AT.NumRows(), A.NumRows() - AT.NumRows(),
               A.OffsetCols(),                A.NumCols());
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Part_2x2_to_3x3(const MatView ATL, const MatView ATR, MatView &A00, MatView &A01, MatView &A02,
                  /***********************************/ MatView &A10, MatView &A11, MatView &A12,
                  const MatView ABL, const MatView ABR, MatView &A20, MatView &A21, MatView &A22,
                  const typename MatView::ordinal_type bm, 
                  const typename MatView::ordinal_type bn,
                  const int quadrant) {
    switch (quadrant) {
    case Partition::TopLeft:
      Part_2x2(ATL, A00, A01,
               /**/ A10, A11, 
               bm, bn, Partition::BottomRight);

      Part_2x1(ATR, A02, 
               /**/ A12,
               bm, Partition::Bottom);

      Part_1x2(ABL, A20, A21,
               bn, Partition::Right);

      A22.setView(ABR.BaseObject(),
                  ABR.OffsetRows(), ABR.NumRows(),
                  ABR.OffsetCols(), ABR.NumCols());
      break;
    case Partition::TopRight:
    case Partition::BottomLeft:
      Kokkos::abort("Tacho::Part_???");
      break;
    case Partition::BottomRight:
      A00.setView(ATL.BaseObject(),
                  ATL.OffsetRows(), ATL.NumRows(),
                  ATL.OffsetCols(), ATL.NumCols());

      Part_1x2(ATR, A01, A02,
               bn, Partition::Left);

      Part_2x1(ABL, A10, 
               /**/ A20,
               bm, Partition::Top);

      Part_2x2(ABR, A11, A12,
               /**/ A21, A22, 
               bm, bn, Partition::TopLeft);
      break;
    default:
      Kokkos::abort("Tacho::Part_???");
      break;
    }
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Part_2x1_to_3x1(const MatView AT, MatView &A0, 
                  /***************/ MatView &A1, 
                  const MatView AB, MatView &A2, 
                  const typename MatView::ordinal_type bm, 
                  const int side) {
    switch (side) {
    case Partition::Top:
      Part_2x1(AT,  A0, 
               /**/ A1,
               bm, Partition::Bottom);

      A2.setView(AB.BaseObject(),
                 AB.OffsetRows(), AB.NumRows(),
                 AB.OffsetCols(), AB.NumCols());
      break;
    case Partition::Bottom:
      A0.setView(AT.BaseObject(),
                 AT.OffsetRows(), AT.NumRows(),
                 AT.OffsetCols(), AT.NumCols());

      Part_2x1(AB,  A1, 
               /**/ A2,
               bm, Partition::Top);
      break;
    default:
      Kokkos::abort("Tacho::Part_???");
      break;
    }
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Part_1x2_to_1x3(const MatView AL, const MatView AR, 
                  MatView &A0, MatView &A1, MatView &A2,
                  const typename MatView::ordinal_type bn, 
                  const int side) {
    switch (side) {
    case Partition::Left:
      Part_1x2(AL,  A0, A1,
               bn, Partition::Right);

      A2.setView(AR.BaseObaject(),
                 AR.OffsetRows(), AR.NumRows(),
                 AR.OffsetCols(), AR.NumCols());
      break;
    case Partition::Right:
      A0.setView(AL.BaseObject(),
                 AL.OffsetRows(), AL.NumRows(),
                 AL.OffsetCols(), AL.NumCols());

      Part_1x2(AR,  A1, A2,
               bn, Partition::Left);
      break;
    default:
      Kokkos::abort("Tacho::Part_???");
      break;
    }
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Merge_2x2(const MatView ATL, const MatView ATR, 
            const MatView ABL, const MatView ABR, MatView &A) {
    A.setView(ATL.BaseObject(),
              ATL.OffsetRows(), ATL.NumRows() + ABR.NumRows(), 
              ATL.OffsetCols(), ATL.NumCols() + ABR.NumCols());
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Merge_1x2(const MatView AL, const MatView AR, MatView &A) {
    A.setView(AL.BaseObject(),
              AL.OffsetRows(), AL.NumRows(),
              AL.OffsetCols(), AL.NumCols() + AR.NumCols());
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Merge_2x1(const MatView AT, 
            const MatView AB, MatView &A) {
    A.setView(AT.BaseObject(),
              AT.OffsetRows(), AT.NumRows() + AB.NumRows(),
              AT.OffsetCols(), AT.NumCols());
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Merge_3x3_to_2x2(const MatView A00, const MatView A01, const MatView A02, MatView &ATL, MatView &ATR, 
                   const MatView A10, const MatView A11, const MatView A12,
                   const MatView A20, const MatView A21, const MatView A22, MatView &ABL, MatView &ABR,
                   const int quadrant) {
    switch (quadrant) {
    case Partition::TopLeft:
      Merge_2x2(A00, A01, 
                A10, A11, ATL);
      
      Merge_2x1(A02, 
                A12, ATR);

      Merge_1x2(A20, A21, ABL);
      
      ABR.setView(A22.BaseObject(),
                  A22.OffsetRows(), A22.NumRows(),
                  A22.OffsetCols(), A22.NumCols());
      break;
    case Partition::TopRight:
    case Partition::BottomLeft:
      Kokkos::abort("Tacho::Part_???");
      break;
    case Partition::BottomRight:
      ATL.setView(A00.BaseObject(),
                  A00.OffsetRows(), A00.NumRows(),
                  A00.OffsetCols(), A00.NumCols());

      Merge_1x2(A01, A02, ATR);

      Merge_2x1(A10, 
                A20, ABL);

      Merge_2x2(A11, A12, 
                A21, A22, ABR);
      break;
    default:
      Kokkos::abort("Tacho::Part_???");
      break;
    }
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Merge_3x1_to_2x1(const MatView A0, MatView &AT, 
                   const MatView A1, 
                   const MatView A2, MatView &AB, 
                   const int side) {
    switch (side) {
    case Partition::Top:
      Merge_2x1(A0, 
                A1, AT);

      AB.setView(A2.BaseObject(),
                 A2.OffsetRows(), A2.NumRows(),
                 A2.OffsetCols(), A2.NumCols());
      break;
    case Partition::Bottom:
      AT.setView(A0.BaseObject(),
                 A0.OffsetRows(), A0.NumRows(),
                 A0.OffsetCols(), A0.NumCols());

      Merge_2x1(A1, 
                A2, AB);
      break;
    default:
      Kokkos::abort("Tacho::Part_???");
      break;
    }
  }

  template<typename MatView>
  KOKKOS_INLINE_FUNCTION 
  void 
  Merge_1x3_to_1x2(const MatView A0, const MatView A1, const MatView A2, 
                   MatView &AL, MatView &AR, 
                   const int side) {
    switch (side) {
    case Partition::Left:
      Merge_1x2(A0, A1, AL);

      AR.setView(A2.BaseObject(),
                 A2.OffsetRows(), A2.NumRows(),
                 A2.OffsetCols(), A2.NumCols());
      break;
    case Partition::Right:
      AL.setView(A0.BaseObject(),
                 A0.OffsetRows(), A0.NumRows(),
                 A0.OffsetCols(), A0.NumCols());

      Merge_1x2(A1, A2, AR);
      break;
    default:
      Kokkos::abort("Tacho::Part_???");
      break;
    }
  }


}

#endif
