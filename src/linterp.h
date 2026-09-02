//
// Copyright (c) 2012 Ronaldo Carpio
//                                     
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and   
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.          
// It is provided "as is" without express or implied warranty.
//               

/*
This is a C++ header-only library for N-dimensional linear interpolation on a rectangular grid. Implements two methods:
* Multilinear: Interpolate using the N-dimensional hypercube containing the point. Interpolation step is O(2^N) 
* Simplicial: Interpolate using the N-dimensional simplex containing the point. Interpolation step is O(N log N), but less accurate.
Updated for C++17 - uses standard library containers and modern C++ features.
*/

#ifndef _linterp_h
#define _linterp_h

#include <cassert>
#include <cmath>
#include <cfloat>
#include <string>
#include <vector>
#include <array>
#include <functional>
#include <memory>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iterator>

using std::vector;
using std::array;
using namespace std::literals;

using uint = unsigned int;
using iVec = vector<int>;
using dVec = vector<double>;

// TODO:
//  - specify behavior past grid boundaries.
//    1) clamp
//    2) return a pre-determined value (e.g. NaN)

struct EmptyClass {};

template <int N, class T, bool CopyData = true, bool Continuous = true, 
          class ArrayRefCountT = EmptyClass, class GridRefCountT = EmptyClass>
class NDInterpolator {
public:
  using value_type = T;
  using array_ref_count_type = ArrayRefCountT;
  using grid_ref_count_type = GridRefCountT;
  
  static constexpr int m_N = N;
  static constexpr bool m_bCopyData = CopyData;
  static constexpr bool m_bContinuous = Continuous;
  
  using grid_type = vector<T>;
  
private:
  // Multi-dimensional array implementation using vector
  template<typename ValueType, std::size_t Dims>
  class MultiArray {
  private:
    vector<ValueType> data_;
    array<std::size_t, Dims> sizes_;
    array<std::size_t, Dims> strides_;
    
    void calculate_strides() {
      if constexpr (Dims > 0) {
        strides_[Dims-1] = 1;
        if constexpr (Dims > 1) {
          for (std::size_t i = Dims-1; i > 0; --i) {
            strides_[i-1] = strides_[i] * sizes_[i];
          }
        }
      }
    }
    
  public:
    MultiArray() = default;
    
    template<typename Iterator>
    MultiArray(Iterator data_begin, Iterator data_end, const array<std::size_t, Dims>& sizes)
      : data_(data_begin, data_end), sizes_(sizes) {
      calculate_strides();
    }
    
    template<typename... Indices>
    const ValueType& operator()(Indices... indices) const {
      static_assert(sizeof...(indices) == Dims, "Wrong number of indices");
      array<std::size_t, Dims> idx_array{static_cast<std::size_t>(indices)...};
      std::size_t flat_index = 0;
      for (std::size_t i = 0; i < Dims; ++i) {
        flat_index += idx_array[i] * strides_[i];
      }
      return data_[flat_index];
    }
    
    template<typename IndexArray>
    const ValueType& operator()(const IndexArray& indices) const {
      std::size_t flat_index = 0;
      for (std::size_t i = 0; i < Dims; ++i) {
        flat_index += indices[i] * strides_[i];
      }
      return data_[flat_index];
    }
  };
  
public:
  using array_type = MultiArray<T, N>;
  using array_type_ptr = std::unique_ptr<array_type>;
  
  array_type_ptr m_pF;
  ArrayRefCountT m_ref_F;                    // reference count for m_pF
  vector<T> m_F_copy;                        // if CopyData == true, this holds the copy of F
     
  vector<grid_type> m_grid_list;    
  vector<GridRefCountT> m_grid_ref_list;     // reference counts for grids  
  vector<vector<T>> m_grid_copy_list;        // if CopyData == true, this holds the copies of the grids
  
  // constructors assume that [f_begin, f_end) is a contiguous array in C-order  
  // non ref-counted constructor.
  template <class IterT1, class IterT2, class IterT3>  
  NDInterpolator(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end) {
    init(grids_begin, grids_len_begin, f_begin, f_end);  
  }
  
  // ref-counted constructor
  template <class IterT1, class IterT2, class IterT3, class RefCountIterT>
  NDInterpolator(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end, 
                 ArrayRefCountT &refF, RefCountIterT grid_refs_begin) {
    init_refcount(grids_begin, grids_len_begin, f_begin, f_end, refF, grid_refs_begin);
  }    
  
  template <class IterT1, class IterT2, class IterT3>                                                          
  void init(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end) {    
    set_grids(grids_begin, grids_len_begin, m_bCopyData);
    set_f_array(f_begin, f_end, m_bCopyData);
  }  
  
  template <class IterT1, class IterT2, class IterT3, class RefCountIterT>
  void init_refcount(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end, 
                     ArrayRefCountT &refF, RefCountIterT grid_refs_begin) {        
    set_grids(grids_begin, grids_len_begin, m_bCopyData);
    set_grids_refcount(grid_refs_begin, grid_refs_begin + N);
    set_f_array(f_begin, f_end, m_bCopyData);
    set_f_refcount(refF);
  }    

  template <class IterT1, class IterT2>  
  void set_grids(IterT1 grids_begin, IterT2 grids_len_begin, bool bCopy) {
    m_grid_list.clear();
    m_grid_ref_list.clear();
    m_grid_copy_list.clear();
    
    for (int i = 0; i < N; i++) {
      int gridLength = grids_len_begin[i];
      if (!bCopy) {        
        const T* grid_ptr = &(*grids_begin[i]);
        m_grid_list.emplace_back(grid_ptr, grid_ptr + gridLength);
      } else {
        m_grid_copy_list.emplace_back(grids_begin[i], grids_begin[i] + grids_len_begin[i]);
        const auto& grid_copy = m_grid_copy_list.back();
        m_grid_list.emplace_back(grid_copy.begin(), grid_copy.end());
      }
    }
  }    
  
  template <class RefCountIterT>  
  void set_grids_refcount(RefCountIterT refs_begin, RefCountIterT refs_end) {
    assert(refs_end - refs_begin == N);    
    m_grid_ref_list.assign(refs_begin, refs_begin + N);
  }    
  
  // assumes that [f_begin, f_end) is a contiguous array in C-order  
  template <class IterT>  
  void set_f_array(IterT f_begin, IterT f_end, bool bCopy) {
    std::size_t nGridPoints = 1;
    array<std::size_t, N> sizes;
    
    for (std::size_t i = 0; i < m_grid_list.size(); i++) {
      sizes[i] = m_grid_list[i].size();
      nGridPoints *= sizes[i];
    }

    std::size_t f_len = std::distance(f_begin, f_end);
    if ((m_bContinuous && f_len != nGridPoints) || (!m_bContinuous && f_len != 2 * nGridPoints)) {
      throw std::invalid_argument("f has wrong size");
    }
    
    if constexpr (!m_bContinuous) {
      for (auto& size : sizes) {
        size *= 2;
      }
    }

    m_F_copy.clear();
    if (!bCopy) {
      m_pF = std::make_unique<array_type>(f_begin, f_end, sizes);
    } else {
      m_F_copy.assign(f_begin, f_end);
      m_pF = std::make_unique<array_type>(m_F_copy.begin(), m_F_copy.end(), sizes);
    }
  }  
  
  void set_f_refcount(ArrayRefCountT &refF) {    
    m_ref_F = refF;
  }
  
  // -1 is before the first grid point
  // N-1 (where grid.size() == N) is after the last grid point
  [[nodiscard]] int find_cell(int dim, T x) const {  
    const auto& grid = m_grid_list[dim];
    if (x < grid.front()) return -1;
    else if (x >= grid.back()) return static_cast<int>(grid.size()) - 1;
    else {
      auto i_upper = std::upper_bound(grid.begin(), grid.end(), x);
      return static_cast<int>(std::distance(grid.begin(), i_upper)) - 1;
    }    
  }
  
  // return the value of f at the given cell and vertex
  [[nodiscard]] T get_f_val(const array<int, N>& cell_index, const array<int, N>& v_index) const {
    array<int, N> f_index;
    
    if constexpr (m_bContinuous) {          
      for (int i = 0; i < N; i++) {
        if (cell_index[i] < 0) {
          f_index[i] = 0;              
        } else if (cell_index[i] >= static_cast<int>(m_grid_list[i].size()) - 1) {
          f_index[i] = static_cast<int>(m_grid_list[i].size()) - 1;              
        } else {
          f_index[i] = cell_index[i] + v_index[i];              
        }
      }
    } else {
      for (int i = 0; i < N; i++) {
        if (cell_index[i] < 0) {
          f_index[i] = 0;
        } else if (cell_index[i] >= static_cast<int>(m_grid_list[i].size()) - 1) {
          f_index[i] = (2 * static_cast<int>(m_grid_list[i].size())) - 1;
        } else {
          f_index[i] = 1 + (2 * cell_index[i]) + v_index[i];
        }
      }
    }
    return (*m_pF)(f_index);
  }
  
  [[nodiscard]] T get_f_val(const array<int, N>& cell_index, int v) const {
    array<int, N> v_index;
    for (int dim = 0; dim < N; dim++) {
      v_index[dim] = (v >> (N - dim - 1)) & 1;    // test if the i-th bit is set
    }
    return get_f_val(cell_index, v_index);
  }    
};

template <int N, class T, bool CopyData = true, bool Continuous = true, 
          class ArrayRefCountT = EmptyClass, class GridRefCountT = EmptyClass>
class InterpSimplex : public NDInterpolator<N, T, CopyData, Continuous, ArrayRefCountT, GridRefCountT> {
public:
  using super = NDInterpolator<N, T, CopyData, Continuous, ArrayRefCountT, GridRefCountT>;
  
  template <class IterT1, class IterT2, class IterT3>  
  InterpSimplex(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end)
    : super(grids_begin, grids_len_begin, f_begin, f_end)
  {}
  
  template <class IterT1, class IterT2, class IterT3, class RefCountIterT>  
  InterpSimplex(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end, 
                ArrayRefCountT &refF, RefCountIterT ref_begins)
    : super(grids_begin, grids_len_begin, f_begin, f_end, refF, ref_begins)
  {}

  template <class IterT>
  [[nodiscard]] T interp(IterT x_begin) const {
    array<T, 1> result;
    array<array<T, 1>, N> coord_iter;
    for (int i = 0; i < N; i++) {
      coord_iter[i][0] = x_begin[i];
    }
    interp_vec(1, coord_iter.begin(), coord_iter.end(), result.begin());
    return result[0];
  }
  
  template <class IterT1, class IterT2>
  void interp_vec(int n, IterT1 coord_iter_begin, IterT1 coord_iter_end, IterT2 i_result) const {
    assert(N == coord_iter_end - coord_iter_begin);
    
    array<int, N> cell_index, v_index;
    array<std::pair<T, int>, N> xipair;    
    
    for (int i = 0; i < n; i++) {            // for each point
      for (int dim = 0; dim < N; dim++) {
        const auto& grid = super::m_grid_list[dim];
        int c = this->find_cell(dim, coord_iter_begin[dim][i]);
        
        T y;
        if (c == -1) {                     // before first grid point
          y = 1.0;
        } else if (c == static_cast<int>(grid.size()) - 1) {    // after last grid point
          y = 0.0;
        } else {
          y = (coord_iter_begin[dim][i] - grid[c]) / (grid[c + 1] - grid[c]);
          y = std::clamp(y, T{0.0}, T{1.0});
        }
        xipair[dim] = {y, dim};        
        cell_index[dim] = c;
      }        
      
      // sort xi's and get the permutation    
      std::sort(xipair.begin(), xipair.end(), 
                [](const auto& a, const auto& b) { return a.first < b.first; });
      
      // walk the vertices of the simplex determined by the permutation  
      std::fill(v_index.begin(), v_index.end(), 1);
      
      T v0 = this->get_f_val(cell_index, v_index);
      T y = v0;
      
      for (int j = 0; j < N; j++) {
        const auto& [xi_val, xi_dim] = xipair[j];  // C++17 structured binding
        v_index[xi_dim]--;        
        T v1 = this->get_f_val(cell_index, v_index);
        y += (1.0 - xi_val) * (v1 - v0);        // interpolate
        v0 = v1;
      }
      *i_result++ = y;
    }    
  }  
};

template <int N, class T, bool CopyData = true, bool Continuous = true, 
          class ArrayRefCountT = EmptyClass, class GridRefCountT = EmptyClass>
class InterpMultilinear : public NDInterpolator<N, T, CopyData, Continuous, ArrayRefCountT, GridRefCountT> {
public:
  using super = NDInterpolator<N, T, CopyData, Continuous, ArrayRefCountT, GridRefCountT>;
  
  template <class IterT1, class IterT2, class IterT3>  
  InterpMultilinear(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end)
    : super(grids_begin, grids_len_begin, f_begin, f_end)
  {}
  
  template <class IterT1, class IterT2, class IterT3, class RefCountIterT>  
  InterpMultilinear(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end, 
                    ArrayRefCountT &refF, RefCountIterT ref_begins)
    : super(grids_begin, grids_len_begin, f_begin, f_end, refF, ref_begins)
  {}

  template <class IterT1, class IterT2>
  [[nodiscard]] static T linterp_nd_unitcube(IterT1 f_begin, IterT1 f_end, IterT2 xi_begin, IterT2 xi_end) {
    int n = static_cast<int>(std::distance(xi_begin, xi_end));
    int f_len = static_cast<int>(std::distance(f_begin, f_end));
    assert((1 << n) == f_len);
    
    if (n == 1) {
      return f_begin[0] + (*xi_begin) * (f_begin[1] - f_begin[0]);
    }
    
    T sub_lower = linterp_nd_unitcube(f_begin, f_begin + (f_len / 2), xi_begin + 1, xi_end);
    T sub_upper = linterp_nd_unitcube(f_begin + (f_len / 2), f_end, xi_begin + 1, xi_end);
    
    return sub_lower + (*xi_begin) * (sub_upper - sub_lower);
  }

  template <class IterT>
  [[nodiscard]] T interp(IterT x_begin) const {
    array<T, 1> result;
    array<array<T, 1>, N> coord_iter;
    for (int i = 0; i < N; i++) {
      coord_iter[i][0] = x_begin[i];
    }
    interp_vec(1, coord_iter.begin(), coord_iter.end(), result.begin());
    return result[0];
  }
  
  template <class IterT1, class IterT2>
  void interp_vec(int n, IterT1 coord_iter_begin, IterT1 coord_iter_end, IterT2 i_result) const {
    assert(N == coord_iter_end - coord_iter_begin);
    
    array<int, N> index;
    vector<T> f(1 << N);
    array<T, N> x;
    
    for (int i = 0; i < n; i++) {                                // loop over each point
      for (int dim = 0; dim < N; dim++) {                        // loop over each dimension
        const auto& grid = super::m_grid_list[dim];        
        int c = this->find_cell(dim, coord_iter_begin[dim][i]);
        
        T y;
        if (c == -1) {                     // before first grid point
          y = 1.0;
        } else if (c == static_cast<int>(grid.size()) - 1) {    // after last grid point
          y = 0.0;
        } else {
          y = (coord_iter_begin[dim][i] - grid[c]) / (grid[c + 1] - grid[c]);
          y = std::clamp(y, T{0.0}, T{1.0});
        }
        index[dim] = c;
        x[dim] = y;
      }
      
      // copy f values at vertices
      for (int v = 0; v < (1 << N); v++) {                    // loop over each vertex of hypercube
        f[v] = this->get_f_val(index, v);
      }
      *i_result++ = linterp_nd_unitcube(f.begin(), f.end(), x.begin(), x.end());
    }
  }
};    

// Type aliases for common use cases
using NDInterpolator_1_S = InterpSimplex<1, double>;
using NDInterpolator_2_S = InterpSimplex<2, double>;
using NDInterpolator_3_S = InterpSimplex<3, double>;
using NDInterpolator_4_S = InterpSimplex<4, double>;
using NDInterpolator_5_S = InterpSimplex<5, double>;

using NDInterpolator_1_ML = InterpMultilinear<1, double>;
using NDInterpolator_2_ML = InterpMultilinear<2, double>;
using NDInterpolator_3_ML = InterpMultilinear<3, double>;
using NDInterpolator_4_ML = InterpMultilinear<4, double>;
using NDInterpolator_5_ML = InterpMultilinear<5, double>;

// C interface
extern "C" {
  void linterp_simplex_1(double **grids_begin, int *grid_len_begin, double *pF, int xi_len, double **xi_begin, double *pResult);
  void linterp_simplex_2(double **grids_begin, int *grid_len_begin, double *pF, int xi_len, double **xi_begin, double *pResult);
  void linterp_simplex_3(double **grids_begin, int *grid_len_begin, double *pF, int xi_len, double **xi_begin, double *pResult);  
}

inline void linterp_simplex_1(double **grids_begin, int *grid_len_begin, double *pF, int xi_len, double **xi_begin, double *pResult) {
  constexpr int N = 1;
  std::size_t total_size = std::accumulate(grid_len_begin, grid_len_begin + N, std::size_t{1}, std::multiplies<>{});
  InterpSimplex<N, double, false> interp_obj(grids_begin, grid_len_begin, pF, pF + total_size);
  interp_obj.interp_vec(xi_len, xi_begin, xi_begin + N, pResult);
}

inline void linterp_simplex_2(double **grids_begin, int *grid_len_begin, double *pF, int xi_len, double **xi_begin, double *pResult) {
  constexpr int N = 2;
  std::size_t total_size = std::accumulate(grid_len_begin, grid_len_begin + N, std::size_t{1}, std::multiplies<>{});
  InterpSimplex<N, double, false> interp_obj(grids_begin, grid_len_begin, pF, pF + total_size);
  interp_obj.interp_vec(xi_len, xi_begin, xi_begin + N, pResult);
}

inline void linterp_simplex_3(double **grids_begin, int *grid_len_begin, double *pF, int xi_len, double **xi_begin, double *pResult) {
  constexpr int N = 3;
  std::size_t total_size = std::accumulate(grid_len_begin, grid_len_begin + N, std::size_t{1}, std::multiplies<>{});
  InterpSimplex<N, double, false> interp_obj(grids_begin, grid_len_begin, pF, pF + total_size);
  interp_obj.interp_vec(xi_len, xi_begin, xi_begin + N, pResult);
}

#endif //_linterp_h