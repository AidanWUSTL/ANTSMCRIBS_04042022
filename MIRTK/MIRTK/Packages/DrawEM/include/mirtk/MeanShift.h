/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
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

#ifndef _MIRTKMEANSHIFT_H

#define _MIRTKMEANSHIFT_H
// queue::push/pop

#include "mirtk/GenericImage.h"
#include "mirtk/Dilation.h"
#include "mirtk/Erosion.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Queue.h"
#include "mirtk/GaussianBlurring.h"

namespace mirtk {

class MeanShift
{
  typedef Vector3D<int> Voxel;

	int _nBins;
	int _padding;
	GreyImage _image;
	Queue<Voxel> _q;
	GreyImage _map, _orig_image;
	GreyImage *_brain;
	GreyImage *_output;
	GreyPixel _imin, _imax;
	double _limit1, _limit2, _limit, _treshold;
	double _bin_width;
	Array<double> _density;
	int _clusterSize;
public:
	double _bg,_wm,_gm,_split1,_split2;


public:

	MeanShift(GreyImage& image, int padding = -1, int nBins = 256);
	~MeanShift();
	void SetOutput( GreyImage *_output);
	int ValueToBin(double value);
	double BinToValue(int bin);
	void AddPoint(int x, int y, int z);
	void AddPoint(int x, int y, int z, int label);
	double msh(double y, double h);
	double findMax(double tr1, double tr2);
	double findMin(double tr1, double tr2);
	double findGMvar();
	double split(double pos1, double pos2, double bw, double h1, double h2);
	double GenerateDensity(double cut_off=0.02);
	void Grow(int x, int y, int z, int label);
	int Lcc(int label, bool add_second = false);
	int LccS(int label, double treshold = 0.5);
	void RemoveBackground();
	void RegionGrowing();
	void FindWMGMmeans();
	void Write(char *output_name);
	void WriteMap(char *output_name);
	void SetTreshold();
	void SetTreshold(double treshold);

	RealImage ReturnMask();
};

}

#endif

