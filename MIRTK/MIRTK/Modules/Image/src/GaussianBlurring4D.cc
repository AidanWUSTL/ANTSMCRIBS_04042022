/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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
 */

#include "mirtk/GaussianBlurring4D.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring4D<VoxelType>::GaussianBlurring4D(double sigma)
:
  GaussianBlurring<VoxelType>(sigma, sigma, sigma, sigma)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring4D<VoxelType>
::GaussianBlurring4D(double xsigma, double ysigma, double zsigma, double tsigma)
:
  GaussianBlurring<VoxelType>(xsigma, ysigma, zsigma, tsigma)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GaussianBlurring4D<VoxelType>::~GaussianBlurring4D()
{
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GaussianBlurring4D<unsigned char>;
template class GaussianBlurring4D<short>;
template class GaussianBlurring4D<unsigned short>;
template class GaussianBlurring4D<float>;
template class GaussianBlurring4D<double>;


} // namespace mirtk
