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
 */

#ifndef MIRTK_InflationStoppingCriterion_H
#define MIRTK_InflationStoppingCriterion_H

#include "mirtk/StoppingCriterion.h"


namespace mirtk {


/**
 * Stopping criterion for cortical surface inflation
 */
class InflationStoppingCriterion : public StoppingCriterion
{
  mirtkObjectMacro(InflationStoppingCriterion);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Convexity/concavity threshold
  mirtkPublicAttributeMacro(double, Threshold);

  /// Normalized average distance of vertex neighbors to the tangent plane
  mirtkReadOnlyAttributeMacro(double, AverageDistance);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const InflationStoppingCriterion &other);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  InflationStoppingCriterion(const ObjectiveFunction * = NULL);

  /// Copy constructor
  InflationStoppingCriterion(const InflationStoppingCriterion &);

  /// Assignment operator
  InflationStoppingCriterion &operator =(const InflationStoppingCriterion &);

  /// Create new copy of this instance
  virtual StoppingCriterion *New() const;

  /// Destructor
  virtual ~InflationStoppingCriterion();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Test surface inflation stopping criterion
  ///
  /// Fischl et al., Cortical Surface-Based Analysis II: Inflation, Flattening,
  /// and a Surface-Based Coordinate System. NeuroImage, 9(2), 195???207 (1999).
  ///
  /// \param[in] iter  Current number of iterations (unused).
  /// \param[in] value Objective function value at currentiteration (unused).
  /// \param[in] delta Last change of objective function parameters (unused).
  ///
  /// \returns Whether stopping criterion is fulfilled.
  virtual bool Fulfilled(int iter, double value, const double *delta);

  // ---------------------------------------------------------------------------
  // Logging

  /// Print current stopping criterion status / value
  ///
  /// This function must be called after Fulfilled, which should update any
  /// cached values that are needed by this function to avoid a costly
  /// reevaluation of the stopping criterion.
  ///
  /// \note The printed string is expected to be considerably short and must
  ///       not end with a newline or space characters.
  virtual void Print(ostream &) const;

};


} // namespace mirtk

#endif // MIRTK_InflationStoppingCriterion_H
