/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#include "mirtk/Erosion.h"
#include "mirtk/BinaryVoxelFunction.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
Erosion<VoxelType>::Erosion()
:
  _Connectivity(ConnectivityType::CONNECTIVITY_26)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
Erosion<VoxelType>::~Erosion()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void Erosion<VoxelType>::Initialize()
{
  ImageToImage<VoxelType>::Initialize();
  _Offsets.Initialize(this->Input(), _Connectivity);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void Erosion<VoxelType>::Run()
{
  this->Initialize();

  const GenericImage<VoxelType> * const input  = this->Input();
  GenericImage<VoxelType>       * const output = this->Output();
  const ImageAttributes                &attr   = input->Attributes();

  ParallelForEachVoxel(BinaryVoxelFunction::Erode(_Offsets), attr, input, output);

  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class Erosion<BytePixel>;
template class Erosion<GreyPixel>;
template class Erosion<RealPixel>;


} // namespace mirtk
