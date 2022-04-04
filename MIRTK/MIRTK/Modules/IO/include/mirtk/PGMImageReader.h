/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_PGMImageReader_H
#define MIRTK_PGMImageReader_H

#include "mirtk/ImageReader.h"


namespace mirtk {


/**
 * Class for reading images in PGM file format.
 *
 * This is a class which reads images in PGM file format and converts them
 * into images. The PGM (portable graymap) file format is a file format for
 * 2D images and is defined in pgm(1). At the moment only images in PGM raw
 * file format are supported.
 */
class PGMImageReader : public ImageReader
{
  mirtkObjectMacro(PGMImageReader);

protected:

  /// Read header of PGM file
  virtual void ReadHeader();

public:

  /// Return whether file has correct header
  static bool CheckHeader(const char *);

  /// Check if this reader can read a given image file
  virtual bool CanRead(const char *) const;

};


} // namespace mirtk

#endif // MIRTK_PGMImageReader_H
