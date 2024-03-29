// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "adjacency_matrix.h"

#include "verbose.h"

#include <vector>

template <typename DerivedF, typename T>
IGL_INLINE void igl::adjacency_matrix(
  const Eigen::MatrixBase<DerivedF> & F, 
  Eigen::SparseMatrix<T>& A)
{
  using namespace std;
  using namespace Eigen;
  typedef typename DerivedF::Scalar Index;

  typedef Triplet<T> IJV;
  vector<IJV > ijv;
  ijv.reserve(F.size()*2);
  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      // Get indices of edge: s --> d
      Index s = F(i,j);
      Index d = F(i,(j+1)%F.cols());
      ijv.push_back(IJV(s,d,1));
      ijv.push_back(IJV(d,s,1));
    }
  }

  const Index n = F.maxCoeff()+1;
  A.resize(n,n);
  switch(F.cols())
  {
    case 3:
      A.reserve(6*(F.maxCoeff()+1));
      break;
    case 4:
      A.reserve(26*(F.maxCoeff()+1));
      break;
  }
  A.setFromTriplets(ijv.begin(),ijv.end());

  // Force all non-zeros to be one

  // Iterate over outside
  for(int k=0; k<A.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (A,k); it; ++it)
    {
      assert(it.value() != 0);
      A.coeffRef(it.row(),it.col()) = 1;
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::adjacency_matrix<Eigen::Matrix<int, -1, -1, 0, -1, -1>, bool>(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<bool, 0, int>&);
template void igl::adjacency_matrix<Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int>&);
template void igl::adjacency_matrix<Eigen::Matrix<int, -1, -1, 0, -1, -1>, int>(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<int, 0, int>&);
#endif
