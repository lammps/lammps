#pragma once
#ifndef __SYMBOLIC_FACTOR_HELPER_HPP__
#define __SYMBOLIC_FACTOR_HELPER_HPP__

/// \file symbolic_factor_helper.hpp
/// \brief The class compute a nonzero pattern with a given level of fills
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"

namespace Tacho {

  using namespace std;

  template<class CrsMatrixType>
  class SymbolicFactorHelper : public Disp {
  public:
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef typename CrsMatrixType::size_type    size_type;

    typedef typename Kokkos::HostSpace::execution_space  host_exec_space ;

    typedef typename CrsMatrixType::ordinal_type_array ordinal_type_array;
    typedef typename CrsMatrixType::size_type_array    size_type_array;
    typedef typename CrsMatrixType::value_type_array   value_type_array;

  private:
    string _label;                   // name of this class

    // matrix index base
    CrsMatrixType _A;                // input matrix
    ordinal_type _m, _n;             // matrix dimension

    struct crs_graph {
      size_type_array _ap;           // row ptr array
      ordinal_type_array _aj;        // col index array
      size_type _nnz;                // # of nonzeros
    };
    typedef struct crs_graph crs_graph_type;
    crs_graph_type _in, _out;

    typedef Kokkos::View<ordinal_type**, Kokkos::LayoutLeft, host_exec_space> league_specific_ordinal_type_array;
    typedef typename league_specific_ordinal_type_array::value_type* league_specific_ordinal_type_array_ptr;

    int _lsize;
    league_specific_ordinal_type_array _queue, _visited, _distance;

    void createInternalWorkSpace() {
      _queue    = league_specific_ordinal_type_array(_label+"::QueueArray",    _m, _lsize);
      _visited  = league_specific_ordinal_type_array(_label+"::VisitedArray",  _m, _lsize);
      _distance = league_specific_ordinal_type_array(_label+"::DistanceArray", _m, _lsize);
    }

    void freeInternalWorkSpace() {
      _queue    = league_specific_ordinal_type_array();
      _visited  = league_specific_ordinal_type_array();
      _distance = league_specific_ordinal_type_array();
    }

  public:

    void setLabel(string label) { _label = label; }
    string Label() const { return _label; }

    SymbolicFactorHelper(const CrsMatrixType &A,
                         const int lsize = (host_exec_space::thread_pool_size(0)/
                                            host_exec_space::thread_pool_size(2)))  {

      _label = "SymbolicFactorHelper::" ;

      // matrix index base and the number of rows
      _A = A;

      _m = _A.NumRows();
      _n = _A.NumCols();

      // allocate memory for input crs matrix
      _in._nnz = _A.NumNonZeros();
      _in._ap  = size_type_array(_label+"::Input::RowPtrArray", _m+1);
      _in._aj  = ordinal_type_array(_label+"::Input::ColIndexArray", _in._nnz);

      // adjust graph structure; A is assumed to have a graph without its diagonal
      A.convertGraph(_in._ap, _in._aj);
      _in._nnz = _in._ap[_m];

      // league size
      _lsize = lsize;

      // create workspace per league
      createInternalWorkSpace();
    }
    virtual~SymbolicFactorHelper() {
      freeInternalWorkSpace();
    }

    class Queue {
    private:
      league_specific_ordinal_type_array_ptr _q;
      ordinal_type _begin, _end;

    public:
      Queue(league_specific_ordinal_type_array_ptr q)
        : _q(q),_begin(0),_end(0) { }

      ordinal_type size() const { return _end - _begin; }
      bool empty() const { return !size(); }

      void push(const ordinal_type val) { _q[_end++] = val; }
      ordinal_type pop() { return _q[_begin++]; }
      ordinal_type end() { return _end; }
      void reset() { _begin = 0; _end = 0; }
    };

    class FunctorComputeNonZeroPatternInRow {
    public:
      typedef Kokkos::TeamPolicy<host_exec_space> policy_type;

    private:
      ordinal_type _level, _m;
      crs_graph_type _graph;

      league_specific_ordinal_type_array _queue;
      league_specific_ordinal_type_array _visited;
      league_specific_ordinal_type_array _distance;

      size_type_array _ap;
      ordinal_type_array _aj;

      ordinal_type _phase;

    public:
      FunctorComputeNonZeroPatternInRow(const ordinal_type level,
                                        const ordinal_type m,
                                        const crs_graph_type &graph,
                                        league_specific_ordinal_type_array &queue,
                                        league_specific_ordinal_type_array &visited,
                                        league_specific_ordinal_type_array &distance,
                                        size_type_array &ap,
                                        ordinal_type_array &aj)
        : _level(level), _m(m), _graph(graph),
          _queue(queue), _visited(visited), _distance(distance),
          _ap(ap), _aj(aj), _phase(0)
      { }

      void setPhaseCountNumNonZeros() { _phase = 0; }
      void setPhaseComputeColIndex()  { _phase = 1; }

      inline
      void operator()(const typename policy_type::member_type &member) const {
        const int lrank = member.league_rank();
        const int lsize = member.league_size();

        league_specific_ordinal_type_array_ptr queue    = &_queue(0, lrank);
        league_specific_ordinal_type_array_ptr distance = &_distance(0, lrank);
        league_specific_ordinal_type_array_ptr visited  = &_visited(0, lrank);

        for (ordinal_type i=0;i<_m;++i)
          visited[i] = 0;

        // shuffle rows to get better load balance;
        // for instance, if ND is applied, more fills are generated in the last seperator.
        for (ordinal_type i=lrank;i<_m;i+=lsize) {

          size_type cnt = 0;

          // account for the diagonal
          switch (_phase) {
          case 0:
            cnt = 1;
            break;
          case 1:
            cnt = _ap[i];
            _aj[cnt++] = i;
            break;
          }

          {
            Queue q(queue); // fixed size queue

            // initialize work space
            q.push(i);
            distance[i] = 0;

            const ordinal_type id = (i+1);
            visited[i] = id;

            // breath first search for i
            while (!q.empty()) {
              const ordinal_type h = q.pop();
              // loop over j adjancy
              const ordinal_type jbegin = _graph._ap[h], jend = _graph._ap[h+1];
              for (ordinal_type j=jbegin;j<jend;++j) {
                const ordinal_type t = _graph._aj[j];
                if (visited[t] != id) {
                  visited[t] = id;

                  if (t < i && (_level < 0 || distance[h] < _level)) {
                    q.push(t);
                    distance[t] = distance[h] + 1;
                  }
                  if (t > i) {
                    switch (_phase) {
                    case 0:
                      ++cnt;
                      break;
                    case 1:
                      _aj[cnt++] = t;
                      break;
                    }
                  }
                }
              }
            }

            // clear work space
            for (ordinal_type j=0;j<q.end();++j) {
              const ordinal_type jj = queue[j];
              distance[jj] = 0;
            }
            q.reset();
          }
          switch (_phase) {
          case 0:
            _ap[i+1] = cnt;
            break;
          case 1:
            sort(_aj.data() + _ap[i] , _aj.data() + _ap[i+1]);
            break;
          }
        }
      }
    };

    class FunctorCountOffsetsInRow {
    public:
      typedef Kokkos::RangePolicy<host_exec_space> policy_type;
      typedef size_type value_type;

    private:
      size_type_array _off_in_rows;

    public:
      FunctorCountOffsetsInRow(size_type_array &off_in_rows)
        : _off_in_rows(off_in_rows)
      { }

      KOKKOS_INLINE_FUNCTION
      void init(value_type &update) const {
        update = 0;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const typename policy_type::member_type &i, value_type &update, const bool final) const {
        update += _off_in_rows(i);
        if (final)
          _off_in_rows(i) = update;
      }

      KOKKOS_INLINE_FUNCTION
      void join(volatile value_type &update,
                volatile const value_type &input) const {
        update += input;
      }
    };

    int createNonZeroPattern(const ordinal_type level,
                             const int uplo,
                             CrsMatrixType &F) {
      // all output array should be local and rcp in Kokkos::View manage memory (de)allocation
      size_type_array ap = size_type_array(_label+"::Output::RowPtrArray", _m+1);

      // later determined
      ordinal_type_array aj;
      value_type_array ax;
      size_type nnz  = 0;

      {
        FunctorComputeNonZeroPatternInRow functor(level, _m, _in,
                                                  _queue,
                                                  _visited,
                                                  _distance,
                                                  ap,
                                                  aj);

        functor.setPhaseCountNumNonZeros();
        Kokkos::parallel_for(typename FunctorComputeNonZeroPatternInRow::policy_type(_lsize, 1), functor);
      }
      {
        FunctorCountOffsetsInRow functor(ap);
        Kokkos::parallel_scan(typename FunctorCountOffsetsInRow::policy_type(0, _m+1), functor);
      }

      nnz  = ap[_m];
      aj = ordinal_type_array(_label+"::Output::ColIndexArray", nnz);
      ax = value_type_array(_label+"::Output::ValueArray", nnz);

      {
        FunctorComputeNonZeroPatternInRow functor(level, _m, _in,
                                                  _queue,
                                                  _visited,
                                                  _distance,
                                                  ap,
                                                  aj);

        functor.setPhaseComputeColIndex();
        Kokkos::parallel_for(typename FunctorComputeNonZeroPatternInRow::policy_type(_lsize, 1), functor);
      }

      {
        F = CrsMatrixType("dummy", _m, _n, nnz, ap, aj, ax);
        F.add(_A);
      }

      // record the symbolic factors
      _out._nnz = nnz;
      _out._ap = ap;
      _out._aj = aj;

      return 0;
    }

    int createNonZeroPattern(const int uplo,
                             CrsMatrixType &F) {
      return createNonZeroPattern(-1, uplo, F);
    }

    ostream& showMe(ostream &os) const {
      streamsize prec = os.precision();
      os.precision(15);
      os << scientific;

      const int w = 6;

      os << " -- Matrix Dimension -- " << endl
         << "    # of Rows  = " << _m << endl
         << "    # of Cols  = " << _n << endl;

      os << endl;

      os << " -- Input Graph Without Diagonals -- " << endl
         << "    # of NonZeros  = " << _in._nnz << endl ;

      os << " -- Input Graph :: RowPtr -- " << endl;
      {
        const ordinal_type n0 = _in._ap.dimension_0();
        for (ordinal_type i=0;i<n0;++i)
          os << setw(w) << i
             << setw(w) << _in._ap[i]
             << endl;
      }

      os << endl;

      os << " -- Output Graph With Diagonals-- " << endl
         << "    # of NonZeros  = " << _out._nnz << endl ;

      os << " -- Output Graph :: RowPtr -- " << endl;
      {
        const ordinal_type n0 = _out._ap.dimension_0();
        for (ordinal_type i=0;i<n0;++i)
          os << setw(w) << i
             << setw(w) << _out._ap[i]
             << endl;
      }

      os.unsetf(ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif



