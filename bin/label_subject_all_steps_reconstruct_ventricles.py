#!/usr/bin/env python3

import nibabel
import numpy
import os
import sys
import skimage.morphology 
import skimage.measure 

if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <initial ventricle labels> <mask> <output>")
    quit()

initVentLabelsFileName = sys.argv[1]
maskFileName = sys.argv[2]

if not os.path.isfile(initVentLabelsFileName):
    print("Init File not found")
    quit()

if not os.path.isfile(maskFileName):
    print("Mask file not found")
    quit()

initVentLabelsNII = nibabel.load(initVentLabelsFileName)
initVentLabelsIMG = (initVentLabelsNII.get_fdata() > 0)

maskNII = nibabel.load(maskFileName)
maskIMG = (maskNII.get_fdata() > 0)

initVentLabelsIMG = numpy.logical_and(initVentLabelsIMG, maskIMG)

reconIMG = skimage.morphology.reconstruction(initVentLabelsIMG, maskIMG, method='dilation', selem=None, offset=None)

reconIMG = skimage.morphology.binary_opening(reconIMG)

# retain the two largest components
all_labels = skimage.measure.label(reconIMG)

for z in range(1, numpy.max(all_labels) + 1):
    if not numpy.any(numpy.logical_and(initVentLabelsIMG, all_labels == z)):
        reconIMG[all_labels == z] = 0
#reconIMG = numpy.reshape(numpy.in1d(all_labels.ravel(), numpy.where(N > 100)[0]), reconIMG.shape)
reconIMG = skimage.morphology.binary_dilation(reconIMG)
reconNII = nibabel.Nifti1Image(numpy.uint8(reconIMG), initVentLabelsNII.affine)
nibabel.save(reconNII, sys.argv[3])
