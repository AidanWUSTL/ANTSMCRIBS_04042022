/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * ATTENTION: This source file has been automatically generated using the code
 *            generator mirtkForEachVoxelFunction.py! This generator is
 *            invoked during CMake configuration of the build system when this
 *            source file is missing from the project.
 *
 *            DO NOT modify this file manually. Instead, modify the code
 *            generator, remove any existing mirtkForEach*VoxelFunction.h
 *            header file from the include/ directory and then re-run CMake.
 *            This will invoke the code generator to re-generate the source files.
 */

#ifndef MIRTK_ForEachSeptenaryVoxelFunction_H
#define MIRTK_ForEachSeptenaryVoxelFunction_H

#include "mirtk/Stream.h"
#include "mirtk/VoxelFunction.h"


namespace mirtk {


inline void _foreachseptenaryvoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
               " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}


// =============================================================================
// 7 const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 7 const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
  const GenericImage<T5> &im5;
  const GenericImage<T6> &im6;
  const GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_Const(const GenericImage<T1> &im1,
                                  const GenericImage<T2> &im2,
                                  const GenericImage<T3> &im3,
                                  const GenericImage<T4> &im4,
                                  const GenericImage<T5> &im5,
                                  const GenericImage<T6> &im6,
                                  const GenericImage<T7> &im7,
                                  VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_Const(const SeptenaryForEachVoxelBody_Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_Const(SeptenaryForEachVoxelBody_Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 7 const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
  const GenericImage<T5> &im5;
  const GenericImage<T6> &im6;
  const GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_Const(const GenericImage<T1> &im1,
                                    const GenericImage<T2> &im2,
                                    const GenericImage<T3> &im3,
                                    const GenericImage<T4> &im4,
                                    const GenericImage<T5> &im5,
                                    const GenericImage<T6> &im6,
                                    const GenericImage<T7> &im7,
                                    VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_Const(const SeptenaryForEachVoxelIfBody_Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_Const(SeptenaryForEachVoxelIfBody_Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, const GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, const GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 6 const, 1 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 6 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_6Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
  const GenericImage<T5> &im5;
  const GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_6Const(const GenericImage<T1> &im1,
                                   const GenericImage<T2> &im2,
                                   const GenericImage<T3> &im3,
                                   const GenericImage<T4> &im4,
                                   const GenericImage<T5> &im5,
                                   const GenericImage<T6> &im6,
                                         GenericImage<T7> &im7,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_6Const(const SeptenaryForEachVoxelBody_6Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_6Const(SeptenaryForEachVoxelBody_6Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 6 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_6Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
  const GenericImage<T5> &im5;
  const GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_6Const(const GenericImage<T1> &im1,
                                     const GenericImage<T2> &im2,
                                     const GenericImage<T3> &im3,
                                     const GenericImage<T4> &im4,
                                     const GenericImage<T5> &im5,
                                     const GenericImage<T6> &im6,
                                           GenericImage<T7> &im7,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_6Const(const SeptenaryForEachVoxelIfBody_6Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_6Const(SeptenaryForEachVoxelIfBody_6Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, const GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, const GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 5 const, 2 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 5 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_5Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
  const GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_5Const(const GenericImage<T1> &im1,
                                   const GenericImage<T2> &im2,
                                   const GenericImage<T3> &im3,
                                   const GenericImage<T4> &im4,
                                   const GenericImage<T5> &im5,
                                         GenericImage<T6> &im6,
                                         GenericImage<T7> &im7,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_5Const(const SeptenaryForEachVoxelBody_5Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_5Const(SeptenaryForEachVoxelBody_5Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 5 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_5Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
  const GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_5Const(const GenericImage<T1> &im1,
                                     const GenericImage<T2> &im2,
                                     const GenericImage<T3> &im3,
                                     const GenericImage<T4> &im4,
                                     const GenericImage<T5> &im5,
                                           GenericImage<T6> &im6,
                                           GenericImage<T7> &im7,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_5Const(const SeptenaryForEachVoxelIfBody_5Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_5Const(SeptenaryForEachVoxelIfBody_5Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, const GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, const GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 4 const, 3 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 4 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_4Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_4Const(const GenericImage<T1> &im1,
                                   const GenericImage<T2> &im2,
                                   const GenericImage<T3> &im3,
                                   const GenericImage<T4> &im4,
                                         GenericImage<T5> &im5,
                                         GenericImage<T6> &im6,
                                         GenericImage<T7> &im7,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_4Const(const SeptenaryForEachVoxelBody_4Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_4Const(SeptenaryForEachVoxelBody_4Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 4 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_4Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_4Const(const GenericImage<T1> &im1,
                                     const GenericImage<T2> &im2,
                                     const GenericImage<T3> &im3,
                                     const GenericImage<T4> &im4,
                                           GenericImage<T5> &im5,
                                           GenericImage<T6> &im6,
                                           GenericImage<T7> &im7,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_4Const(const SeptenaryForEachVoxelIfBody_4Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_4Const(SeptenaryForEachVoxelIfBody_4Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 3 const, 4 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 3 const, 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_3Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
        GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_3Const(const GenericImage<T1> &im1,
                                   const GenericImage<T2> &im2,
                                   const GenericImage<T3> &im3,
                                         GenericImage<T4> &im4,
                                         GenericImage<T5> &im5,
                                         GenericImage<T6> &im6,
                                         GenericImage<T7> &im7,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_3Const(const SeptenaryForEachVoxelBody_3Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_3Const(SeptenaryForEachVoxelBody_3Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 3 const, 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_3Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
        GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_3Const(const GenericImage<T1> &im1,
                                     const GenericImage<T2> &im2,
                                     const GenericImage<T3> &im3,
                                           GenericImage<T4> &im4,
                                           GenericImage<T5> &im5,
                                           GenericImage<T6> &im6,
                                           GenericImage<T7> &im7,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_3Const(const SeptenaryForEachVoxelIfBody_3Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_3Const(SeptenaryForEachVoxelIfBody_3Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 2 const, 5 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 2 const, 5 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_2Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_2Const(const GenericImage<T1> &im1,
                                   const GenericImage<T2> &im2,
                                         GenericImage<T3> &im3,
                                         GenericImage<T4> &im4,
                                         GenericImage<T5> &im5,
                                         GenericImage<T6> &im6,
                                         GenericImage<T7> &im7,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_2Const(const SeptenaryForEachVoxelBody_2Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_2Const(SeptenaryForEachVoxelBody_2Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 2 const, 5 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_2Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_2Const(const GenericImage<T1> &im1,
                                     const GenericImage<T2> &im2,
                                           GenericImage<T3> &im3,
                                           GenericImage<T4> &im4,
                                           GenericImage<T5> &im5,
                                           GenericImage<T6> &im6,
                                           GenericImage<T7> &im7,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_2Const(const SeptenaryForEachVoxelIfBody_2Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_2Const(SeptenaryForEachVoxelIfBody_2Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 1 const, 6 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 const, 6 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody_1Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
        GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody_1Const(const GenericImage<T1> &im1,
                                         GenericImage<T2> &im2,
                                         GenericImage<T3> &im3,
                                         GenericImage<T4> &im4,
                                         GenericImage<T5> &im5,
                                         GenericImage<T6> &im6,
                                         GenericImage<T7> &im7,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody_1Const(const SeptenaryForEachVoxelBody_1Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody_1Const(SeptenaryForEachVoxelBody_1Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 const, 6 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody_1Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
        GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;
        GenericImage<T5> &im5;
        GenericImage<T6> &im6;
        GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody_1Const(const GenericImage<T1> &im1,
                                           GenericImage<T2> &im2,
                                           GenericImage<T3> &im3,
                                           GenericImage<T4> &im4,
                                           GenericImage<T5> &im5,
                                           GenericImage<T6> &im6,
                                           GenericImage<T7> &im7,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody_1Const(const SeptenaryForEachVoxelIfBody_1Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody_1Const(SeptenaryForEachVoxelIfBody_1Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 7 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 7 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct SeptenaryForEachVoxelBody : public ForEachVoxelBody<VoxelFunc>
{
  GenericImage<T1> &im1;
  GenericImage<T2> &im2;
  GenericImage<T3> &im3;
  GenericImage<T4> &im4;
  GenericImage<T5> &im5;
  GenericImage<T6> &im6;
  GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelBody(GenericImage<T1> &im1,
                            GenericImage<T2> &im2,
                            GenericImage<T3> &im3,
                            GenericImage<T4> &im4,
                            GenericImage<T5> &im5,
                            GenericImage<T6> &im6,
                            GenericImage<T7> &im7,
                            VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelBody(const SeptenaryForEachVoxelBody &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelBody(SeptenaryForEachVoxelBody &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<SeptenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 7 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct SeptenaryForEachVoxelIfBody : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  GenericImage<T1> &im1;
  GenericImage<T2> &im2;
  GenericImage<T3> &im3;
  GenericImage<T4> &im4;
  GenericImage<T5> &im5;
  GenericImage<T6> &im6;
  GenericImage<T7> &im7;

  /// Constructor
  SeptenaryForEachVoxelIfBody(GenericImage<T1> &im1,
                              GenericImage<T2> &im2,
                              GenericImage<T3> &im3,
                              GenericImage<T4> &im4,
                              GenericImage<T5> &im5,
                              GenericImage<T6> &im6,
                              GenericImage<T7> &im7,
                              VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  SeptenaryForEachVoxelIfBody(const SeptenaryForEachVoxelIfBody &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  SeptenaryForEachVoxelIfBody(SeptenaryForEachVoxelIfBody &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<SeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<SeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  SeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, GenericImage<T5> *im5, GenericImage<T6> *im6, GenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  SeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, GenericImage<T5> &im5, GenericImage<T6> &im6, GenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _foreachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}


} // namespace mirtk

#endif
