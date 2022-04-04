/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_IntrinsicLeastEdgeLengthDistortionSurfaceMapper_H
#define MIRTK_IntrinsicLeastEdgeLengthDistortionSurfaceMapper_H

#include "mirtk/NearOptimalIntrinsicSurfaceMapper.h"


namespace mirtk {


/**
 * Compute a near-optimal intrinsic surface map which minimizes edge-length distortion
 *
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 */
class IntrinsicLeastEdgeLengthDistortionSurfaceMapper : public NearOptimalIntrinsicSurfaceMapper
{
  mirtkObjectMacro(IntrinsicLeastEdgeLengthDistortionSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const IntrinsicLeastEdgeLengthDistortionSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  IntrinsicLeastEdgeLengthDistortionSurfaceMapper();

  /// Copy constructor
  IntrinsicLeastEdgeLengthDistortionSurfaceMapper(
    const IntrinsicLeastEdgeLengthDistortionSurfaceMapper &
  );

  /// Assignment operator
  IntrinsicLeastEdgeLengthDistortionSurfaceMapper &operator =(
    const IntrinsicLeastEdgeLengthDistortionSurfaceMapper &
  );

  /// Destructor
  virtual ~IntrinsicLeastEdgeLengthDistortionSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Compute affine combination weight that minimizes a given distortion measure
  ///
  /// \param[in] u0 Discrete authalic  surface map values, i.e., lambda=0.
  /// \param[in] u1 Discrete conformal surface map values, i.e., lambda=1.
  ///
  /// \returns Affine combination weight \f$\lambda\f$, where final surface map values
  ///          are computed as \f$\lambda u1 + (1 - \lambda) * u0\f$.
  virtual double ComputeLambda(vtkDataArray *u0, vtkDataArray *u1) const;

};


} // namespace mirtk

#endif // MIRTK_IntrinsicLeastEdgeLengthDistortionSurfaceMapper_H
