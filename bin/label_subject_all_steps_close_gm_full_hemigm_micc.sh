#!/bin/bash
#SBATCH --job-name=skullstrip
#SBATCH -n1
#SBATCH -N1
#SBATCH --time=24:00:00
#SBATCH --mem=6GB
#SBATCH --cpus-per-task=2

if [ -z "${1}" ]
then
	echo "$0 <subjid>"
	exit
fi

# register the template to the target non-linearly
#       SBATCH --constraint=SouthLevel1

export USET1=NO
export FORCEOVERWRITE=NO
export NOSKULLSTRIP=NO
export LARGEVENTRICLES=NO
export DOPARALLEL=NO

while (( "$#" ))
do
	case $1 in
		-T1)
			export USET1=YES
			;;
		-force)
			export FORCEOVERWRITE=YES
			;;
		-noskullstrip)
			export NOSKULLSTRIP=YES
			;;
		-largeventricles)
			export LARGEVENTRICLES=YES
			;;
		-parallel)
			export DOPARALLEL=YES
			;;
		-h)
			echo "$0 [options] <subject id>"
			echo
			echo "Options:"
			echo -e "\t-force: Force overwrite"
			echo -e "\t-largeventricles: Large ventricles (unused at the moment)"
			echo -e "\t-noskullstrip: Input T2 is already skull stripped (unusued)"
			echo -e "\t-T1: Use T1 in deformable (unusued)"
			echo -e "\t-parallel: Perform registrations to training data in parallel"
			;;
		*)
			export SUBJID=$1
		;;
	esac
	shift;
done

if [ -z "$SUBJID" ]
then

fi
H=`hostname`
export TISSUESEGDIR=TissueSegMCRIBS
mkdir -p ${TISSUESEGDIR}/${SUBJID}

export TEMPLATEDIR=/home/addo/MCRIownCloud/deve2-chris.adamson/neonatal/OrigImagesLabelledLaPrem/ANTST1T2TemplateGMAIMIHighGMWeightDemons
#export TEMPLATEDIR=/group/deve2/data/addo/neonatal/OrigImagesLabelledLaPrem/ANTST1T2TemplateGMAIMIHighGMWeightDemons
export T2TEMPLATE=$TEMPLATEDIR/Finaltemplate0.nii.gz
export T2TEMPLATELAPLACIAN=$TEMPLATEDIR/Finaltemplate0Laplacian.nii.gz
export T1TEMPLATE=$TEMPLATEDIR/Finaltemplate1.nii.gz
export T2TEMPLATEBRAIN=$TEMPLATEDIR/Finaltemplate0Brain.nii.gz
export T1TEMPLATEBRAIN=$TEMPLATEDIR/Finaltemplate1Brain.nii.gz

#TEMPLATEDIR=$HOME/MCRIownCloud/deve2-chris.adamson/neonatal/BrainMaskAtlas/
FSDIR=../LaPrem/freesurfer

export T2TARGET=../RawT2RadiologicalIsotropicCropped/${SUBJID}.nii.gz
export T1TARGET=../RawT1RadiologicalIsotropicCropped/${SUBJID}.nii.gz

#if [ "$USET1" == "YES" ]
#then
#	fslcpgeom $T2TARGET $T1TARGET
#fi

export OUTPUTPREFIX=${TISSUESEGDIR}/${SUBJID}/${SUBJID}

#rm -fr ${TISSUESEGDIR}/${SUBJID}
mkdir -p ${TISSUESEGDIR}/${SUBJID}

NUMPROC=`nproc`
export NUMTHREADS=`expr $NUMPROC / 10`
if [ "$NUMTHREADS" == "0" ]
then
	NUMTHREADS=1
fi
#echo $NUMTHREADS
if [ "$H" == "beast" ]
then
	NUMTHREADS=4
fi

#rm -f ${OUTPUTPREFIX}BrainExtraction*
#if [ ! -f "${OUTPUTPREFIX}_brain_mask.nii.gz" -o "$FORCEOVERWRITE" == "YES" ]
PADAMOUNT=5
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=`nproc`
#export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
export OMP_NUM_THREADS=`nproc`


#DOPARALLEL=YES
#export OMP_NUM_THREADS=1
#if [ -f "${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz" ]
#then
#	echo "$SUBJID found"
#	exit
#else
#	echo "$SUBJID not found"
#fi
if [ ! -f "${OUTPUTPREFIX}_t2w_init_restore.nii.gz" -o "$FORCEOVERWRITE" == "YES" ]
then
	ImageMath 3 ${OUTPUTPREFIX}_t2w_rescaled.nii.gz RescaleImage $T2TARGET 10 100
	N4BiasFieldCorrection --verbose --input-image ${OUTPUTPREFIX}_t2w_rescaled.nii.gz --image-dimensionality 3 -s 4 --output ${OUTPUTPREFIX}_t2w_init_restore.nii.gz
	#p ${OUTPUTPREFIX}_t2w_rescaled.nii.gz ${OUTPUTPREFIX}_t2w_init_restore.nii.gz
	#NINETYFIVET2=`fslstats  ${OUTPUTPREFIX}_t2w_init_restore.nii.gz -P 95`
	#mageMath 3 ${OUTPUTPREFIX}_t2w_init_restore.nii.gz WindowImage ${OUTPUTPREFIX}_t2w_init_restore.nii.gz 0 $NINETYFIVET2 0 1000
	RescaleNinetyFivePercentile ${OUTPUTPREFIX}_t2w_init_restore.nii.gz ${OUTPUTPREFIX}_t2w_init_restore.nii.gz
	ImageMath 3 ${OUTPUTPREFIX}_t2w_init_restore.nii.gz PadImage ${OUTPUTPREFIX}_t2w_init_restore.nii.gz $PADAMOUNT
	rm -f ${OUTPUTPREFIX}_t2w_rescaled.nii.gz
fi
if [ ! -f "${OUTPUTPREFIX}InitialAffine.mat" ]
then
	ResampleImageBySpacing 3 $T2TEMPLATE ${OUTPUTPREFIX}T2TemplateDown.nii.gz 3 3 3 1 0
	ResampleImageBySpacing 3 ${OUTPUTPREFIX}_t2w_init_restore.nii.gz ${OUTPUTPREFIX}_t2w_init_restore_down.nii.gz 3 3 3 1 0
	antsAI -d 3 -v 1 \
		-m Mattes[${OUTPUTPREFIX}T2TemplateDown.nii.gz,${OUTPUTPREFIX}_t2w_init_restore_down.nii.gz,32,Random,0.2 ] \
		-t Affine[ 0.1 ] \
		-s [ 20,0.12 ] \
		-g [ 40,0x40x40 ] \
		-p 0 \
		-c 10 \
		-o ${OUTPUTPREFIX}InitialAffine.mat
		#-m Mattes[RawT1RadiologicalIsotropicCropped/${i}.nii.gz,$T1TEMPLATE,32,Random,0.05 ] \
	rm -f ${OUTPUTPREFIX}T2TemplateDown.nii.gz ${OUTPUTPREFIX}_t2w_init_restore_down.nii.gz
fi

#if [ ! -f "${OUTPUTPREFIX}_majority_dkt_init_reg.nii.gz" ]
#then
#	antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateDKTMajority.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz \
#	       --transform [${OUTPUTPREFIX}InitialAffine.mat,1] \
#	       --interpolation GenericLabel \
#	       --output-data-type short \
#	       --output ${OUTPUTPREFIX}_majority_dkt_init_reg.nii.gz
#fi
#if [ ! -f "${OUTPUTPREFIX}_template_t2w_restore_initaffine_reg.nii.gz" ]
#then
#	antsApplyTransforms -v -d 3 --reference-image $TEMPLATEDIR/FinaltemplateDKTMajority.nii.gz --input ${OUTPUTPREFIX}_t2w_restore.nii.gz \
#	       --transform [${OUTPUTPREFIX}InitialAffine.mat,0] \
#	       --interpolation Linear \
#	       --output-data-type short \
#	       --output ${OUTPUTPREFIX}_template_t2w_restore_initaffine_reg.nii.gz
#fi
#

if [ ! -f "${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat" ]
then
	#-x $TEMPLATEDIR/FinaltemplateBrainProbRegistrationMask.nii.gz \
	antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1  \
		--initial-moving-transform [${OUTPUTPREFIX}InitialAffine.mat,0] \
		--transform Rigid[ 0.1 ] --metric MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,32 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --masks $TEMPLATEDIR/FinaltemplateBrainProbRegistrationMask.nii.gz \
		--transform Affine[ 0.1 ] --metric MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,32 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --masks $TEMPLATEDIR/FinaltemplateBrainProbRegistrationMask.nii.gz \
		--output ${OUTPUTPREFIX}_skullstrip_affine
#--output [${OUTPUTPREFIX}_affine,${OUTPUTPREFIX}_affineWarped.nii.gz]
		#-m MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,32 ] \
		#--transform Rigid[ 0.1 ] --metric MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_restore.nii.gz,1,32 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --masks $TEMPLATEDIR/FinaltemplateBrainMask.nii.gz \
		#--transform Affine[ 0.1 ] --metric MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_restore.nii.gz,1,32 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --masks $TEMPLATEDIR/FinaltemplateBrainMask.nii.gz \
fi

#fslmaths ${OUTPUTPREFIX}_skullstrip_reg_affineWarped.nii.gz -bin -mul $TEMPLATEDIR/FinaltemplateBrainProbRegistrationMask.nii.gz ${OUTPUTPREFIX}_skullstrip_reg_NonLinearMask.nii.gz -odt char
antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/Finaltemplate0.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_init_restore.nii.gz \
	--transform [${OUTPUTPREFIX}InitialAffine.mat,1] \
		--output-data-type float \
	--output ${OUTPUTPREFIX}_template_initaffine_reg.nii.gz
antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/Finaltemplate0.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_init_restore.nii.gz \
	--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
		--output-data-type float \
	--output ${OUTPUTPREFIX}_template_affine_reg.nii.gz

#./GaussianLaplacian ${OUTPUTPREFIX}_t2w_restore.nii.gz ${OUTPUTPREFIX}_t2w_restore_laplacian2.nii.gz

if [ ! -f "${OUTPUTPREFIX}_skullstrip_reg1Warp.nii.gz" ]
then
#	rm -f Stage*.nii.gz
#antsRegistration -v -d 3 -u 1 --verbose 1 --float 1 --collapse-output-transforms 1 --write-interval-volumes 10 \
	antsRegistration -v -d 3 -u 1 --verbose 1 --float 1 --collapse-output-transforms 1 \
		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat \
		--transform BSplineSyN[0.13,26,0,3] \
		--metric MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,32 ] \
		--metric CC[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,2 ] \
		--convergence [ 200x200x200x0x0,1e-6,40 ] --shrink-factors 8x6x4x2x1 --smoothing-sigmas 6x4x2x2x1vox \
		--transform SyN[0.13,3,0] \
		--metric MI[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,32 ] \
		--metric CC[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,2 ] \
		--convergence [ 100x0,1e-6,40 ] --shrink-factors 2x1 --smoothing-sigmas 2x1vox \
		--output ${OUTPUTPREFIX}_skullstrip_reg
		
#--metric CC[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_init_restore.nii.gz,1,4 ] \
#--output [${OUTPUTPREFIX}_ants_reg,${OUTPUTPREFIX}_ants_regWarped.nii.gz]
#fslmerge -a ${OUTPUTPREFIX}_ants_reg_all_iters Stage*.nii.gz
#	rm -f Stage*.nii.gz
		#--metric CC[ $TEMPLATEDIR/Finaltemplate0Laplacian2.nii.gz,${OUTPUTPREFIX}_t2w_restore_laplacian2.nii.gz,0.8,4 ] \
	#antsRegistration -v -d 3 -u 1 --verbose 1 --float 1 --collapse-output-transforms 1 \
	#	--initial-moving-transform ${OUTPUTPREFIX}_affine0GenericAffine.mat \
	#	--transform BSplineSyN[0.1,26,0,3] \
	#	--metric CC[ $T2TEMPLATE,${OUTPUTPREFIX}_t2w_restore.nii.gz,1,4 ] \
	#	--convergence [ 200x100x0,1e-6,40 ] --shrink-factors 4x2x1 --smoothing-sigmas 2x2x1vox \
	#	--output [${OUTPUTPREFIX}_ants_reg,${OUTPUTPREFIX}_ants_regWarped.nii.gz]
fi

antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/Finaltemplate0.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz \
	--transform [${OUTPUTPREFIX}_skullstrip_reg0GenericAffine.mat,1] \
	--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
		--output-data-type float \
	--output ${OUTPUTPREFIX}_template_skullstrip_reg.nii.gz

#if [ ! -f "${OUTPUTPREFIX}_brain_mask.nii.gz" ]
#then
	antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateBrainMask.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_init_restore.nii.gz \
		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
		--output-data-type float \
		--output ${OUTPUTPREFIX}_brainmask_skullstrip_reg.nii.gz

	./MRIBinarize --i ${OUTPUTPREFIX}_brainmask_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_brain_mask.nii.gz --min 0.25 --noverbose 
	rm -f ${OUTPUTPREFIX}_brainmask_skullstrip_reg.nii.gz
#fi
#mri_mask ${OUTPUTPREFIX}_t2w_init_restore.nii.gz ${OUTPUTPREFIX}_brain_mask.nii.gz ${OUTPUTPREFIX}_t2w_init_restore_brain.nii.gz 
fslcpgeom ${OUTPUTPREFIX}_t2w_init_restore.nii.gz ${OUTPUTPREFIX}_brain_mask.nii.gz

if [ ! -f "${OUTPUTPREFIX}_t2w_restore.nii.gz" -o "$FORCEOVERWRITE" == "YES" ]
then
	ImageMath 3 ${OUTPUTPREFIX}_t2w_rescaled.nii.gz RescaleImage $T2TARGET 10 100
	ImageMath 3 ${OUTPUTPREFIX}_brain_mask_unpadded.nii.gz PadImage ${OUTPUTPREFIX}_brain_mask.nii.gz -$PADAMOUNT
	fslcpgeom ${OUTPUTPREFIX}_t2w_rescaled.nii.gz ${OUTPUTPREFIX}_brain_mask_unpadded.nii.gz
	N4BiasFieldCorrection --verbose --input-image ${OUTPUTPREFIX}_t2w_rescaled.nii.gz --image-dimensionality 3 -s 2 --output ${OUTPUTPREFIX}_t2w_init_restore2.nii.gz -x ${OUTPUTPREFIX}_brain_mask_unpadded.nii.gz
	#p ${OUTPUTPREFIX}_t2w_rescaled.nii.gz ${OUTPUTPREFIX}_t2w_init_restore.nii.gz
	#NINETYFIVET2=`fslstats  ${OUTPUTPREFIX}_t2w_init_restore.nii.gz -P 95`
	#mageMath 3 ${OUTPUTPREFIX}_t2w_init_restore.nii.gz WindowImage ${OUTPUTPREFIX}_t2w_init_restore.nii.gz 0 $NINETYFIVET2 0 1000
	RescaleNinetyFivePercentile ${OUTPUTPREFIX}_t2w_init_restore2.nii.gz ${OUTPUTPREFIX}_t2w_init_restore2.nii.gz
	ImageMath 3 ${OUTPUTPREFIX}_t2w_restore.nii.gz PadImage ${OUTPUTPREFIX}_t2w_init_restore2.nii.gz $PADAMOUNT
	rm -f ${OUTPUTPREFIX}_t2w_rescaled.nii.gz ${OUTPUTPREFIX}_t2w_init_restore2.nii.gz
fi
fslmaths ${OUTPUTPREFIX}_t2w_restore.nii.gz -mas ${OUTPUTPREFIX}_brain_mask.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz 
fslcpgeom ${OUTPUTPREFIX}_t2w_restore ${OUTPUTPREFIX}_t2w_restore_brain
if [ ! -f "${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz" ]
then
	DenoiseImage -d 3 -i ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz -x ${OUTPUTPREFIX}_brain_mask.nii.gz -o ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -v 1
fi
if [ "$USET1" == "YES" ]
then
	if [ ! -f "${OUTPUTPREFIX}_t1w_restore.nii.gz" -o "$FORCEOVERWRITE" == "YES" ]
	then
		ImageMath 3 ${OUTPUTPREFIX}_t1w_rescaled.nii.gz RescaleImage $T1TARGET 10 100
		N4BiasFieldCorrection --verbose --input-image ${OUTPUTPREFIX}_t1w_rescaled.nii.gz --image-dimensionality 3 -s 4 --output ${OUTPUTPREFIX}_t1w_init_restore.nii.gz
		RescaleNinetyFivePercentile ${OUTPUTPREFIX}_t1w_init_restore.nii.gz ${OUTPUTPREFIX}_t1w_init_restore.nii.gz
		ImageMath 3 ${OUTPUTPREFIX}_t1w_restore.nii.gz PadImage ${OUTPUTPREFIX}_t1w_init_restore.nii.gz $PADAMOUNT
		rm -f ${OUTPUTPREFIX}_t1w_rescaled.nii.gz ${OUTPUTPREFIX}_t1w_init_restore.nii.gz
	fi
	fslmaths ${OUTPUTPREFIX}_t1w_restore -mas ${OUTPUTPREFIX}_brain_mask ${OUTPUTPREFIX}_t1w_restore_brain
fi

if [ ! -f "${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz" ]
then
	./RobustOtsu2 ${OUTPUTPREFIX}_t2w_restore.nii.gz ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz
	ImageMath 3 ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz GetLargestComponent ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz
	./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz --match 1 --dilate 3 --noverbose 
	fslmaths ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz ${OUTPUTPREFIX}_t2w_restore_regmask.nii.gz -odt char
	fslcpgeom ${OUTPUTPREFIX}_t2w_restore ${OUTPUTPREFIX}_t2w_restore_regmask	
fi

#if [ ! -f "${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz" ]
#then
#	#Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o [${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_posterior%01d.nii.gz] -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 5 ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
#	Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 5 ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
#fi
#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz -odt char
#
#if [ ! -f "${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos4_segmentation.nii.gz" ]
#then
#	#Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o [${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_posterior%01d.nii.gz] -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 5 ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
#	Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos4_segmentation.nii.gz -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 4 ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
#fi
#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos4_segmentation.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos4_segmentation.nii.gz -odt char
#
for j in `seq 1000 1035`
do
        GMMATCH="$GMMATCH $j `expr $j + 1000`"
done

#if [ ! -f "${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz" ]
#then
	antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateDKTWithSkullLabelMajority.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
		--transform [${OUTPUTPREFIX}_skullstrip_reg0GenericAffine.mat,1] \
		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
		--interpolation GenericLabel \
		--output-data-type short \
		--output ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz
#fi
#if [ ! -f "${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg.nii.gz" ]
#then
	antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateDKTWithLatVentRingsMajority.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
		--transform [${OUTPUTPREFIX}_skullstrip_reg0GenericAffine.mat,1] \
		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
		--interpolation GenericLabel \
		--output-data-type short \
		--output ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg.nii.gz
#fi

#if [ ! -f "${OUTPUTPREFIX}_template_skullstrip_reg.nii.gz" ]
#then
#	antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/Finaltemplate0.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_reg0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--output-data-type float \
#		--output ${OUTPUTPREFIX}_template_skullstrip_reg.nii.gz
#fi

fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_brain_mask.nii.gz 
if [ ! -f "${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz" ]
then
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz --match 2 41 --noverbose 
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_gm.nii.gz --match $GMMATCH 9 48 51 52 12 13 --noverbose 
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz --match 4 43 24 --noverbose 

	CSFMEAN=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -k ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz -m`
	GMMEAN=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -k ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_gm.nii.gz -m`
	WMMEAN=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -k ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz -m`

	rm -f ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_gm.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz
	#Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o [${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_posterior%01d.nii.gz] -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 5 ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
	Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 3,${GMMEAN}x${WMMEAN}x${CSFMEAN} ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
fi
fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz -odt char
#			antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/Finaltemplate0.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_reg_cc0GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg_cc1InverseWarp.nii.gz \
#				--output ${OUTPUTPREFIX}_template_skullstrip_reg_cc.nii.gz
#			
#			antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateDKTWithSkullLabelMajority.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_reg_cc0GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg_cc1InverseWarp.nii.gz \
#				--interpolation GenericLabel \
#				--output-data-type short \
#				--output ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cc.nii.gz
#			antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateBrainProbRegistrationMaskNoLatVent.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_reg_cc0GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg_cc1InverseWarp.nii.gz \
#				--interpolation GenericLabel \
#				--output-data-type short \
#				--output ${OUTPUTPREFIX}_FinaltemplateBrainProbRegistrationMaskNoLatVent_skullstrip_reg_cc.nii.gz
#			antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/Finaltemplate0.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_reg_cc20GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg_cc21InverseWarp.nii.gz \
#				--output ${OUTPUTPREFIX}_template_skullstrip_reg_cc2.nii.gz
#			
#			antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateDKTWithSkullLabelMajority.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_reg_cc20GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg_cc21InverseWarp.nii.gz \
#				--interpolation GenericLabel \
#				--output-data-type short \
#				--output ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cc2.nii.gz
#			antsApplyTransforms -v -d 3 --input $TEMPLATEDIR/FinaltemplateBrainProbRegistrationMaskNoLatVent.nii.gz --reference-image ${OUTPUTPREFIX}_t2w_restore_brain.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_reg_cc20GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg_cc21InverseWarp.nii.gz \
#				--interpolation GenericLabel \
#				--output-data-type short \
#				--output ${OUTPUTPREFIX}_FinaltemplateBrainProbRegistrationMaskNoLatVent_skullstrip_reg_cc2.nii.gz
#
#
#		exit		
#fi	
	#for i in `seq -w 1 10`
	#do
	#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_t2_dkt_dark_wm.nii.gz \
	#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
	#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
	#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
	#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
	#		--interpolation GenericLabel \
	#		--output-data-type short \
	#		--output $TISSUESEGDIR/$SUBJID/P${i}_t2_dkt_dark_wm_to_$SUBJID.nii.gz
	#done
	#if [ ! -f "${OUTPUTPREFIX}_t2w_dark_wm_prob_skullstrip_reg.nii.gz" ]
	#then
		antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/FinaltemplateT2DarkWMProb.nii.gz \
			--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
			--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
			--interpolation Linear \
			--output-data-type float \
			--output ${OUTPUTPREFIX}_t2w_dark_wm_prob_skullstrip_reg.nii.gz
	#fi
#		ATLAST2IMAGESNAMES=
#		ATLAST2IMAGES=
#		ATLASDKTIMAGES=
#		ATLASDKTIMAGESNAMES=	
#		for i in `seq -w 1 10`
#		do
#			antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_t2.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#				--interpolation Linear \
#				--output-data-type float \
#				--output $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_${SUBJID}_skullstrip_reg.nii.gz
#			antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_dkt_with_skull_label.nii.gz \
#				--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#				--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#				--interpolation GenericLabel \
#				--output-data-type short \
#				--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_${SUBJID}_skullstrip_reg.nii.gz
#			ATLAST2IMAGESNAMES="$ATLAST2IMAGESNAMES $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_${SUBJID}_skullstrip_reg.nii.gz"
#			ATLAST2IMAGES="$ATLAST2IMAGES -g $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_${SUBJID}_skullstrip_reg.nii.gz"
#			ATLASDKTIMAGES="$ATLASDKTIMAGES -l $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_${SUBJID}_skullstrip_reg.nii.gz"
#			ATLASDKTIMAGESNAMES="$ATLASDKTIMAGESNAMES $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_${SUBJID}_skullstrip_reg.nii.gz"
#			fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_${SUBJID}_skullstrip_reg
#			fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_${SUBJID}_skullstrip_reg
#
#		
#		done
#		#if [ ! -f "${OUTPUTPREFIX}_labelfusionimage_skullstrip_reg_dkt_antsinit.nii.gz" ]
#		#then
#			antsJointFusion $ATLAST2IMAGES $ATLASDKTIMAGES --target-image ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -x ${OUTPUTPREFIX}_brain_mask.nii.gz --output ${OUTPUTPREFIX}_labelfusionimage_skullstrip_reg_dkt_antsinit.nii.gz -v --patch-metric PC -s 3 -p 3
#		#fi
#		fslmaths ${OUTPUTPREFIX}_labelfusionimage_skullstrip_reg_dkt_antsinit.nii.gz  ${OUTPUTPREFIX}_labelfusionimage_skullstrip_reg_dkt_antsinit.nii.gz  -odt short
#		exit
#	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_not_3rd_ventricle.nii.gz --match 14 24 170 --inv --erode 1 --noverbose 
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_lateral_ventricles.nii.gz --match 4 43 31 63 --noverbose 
	#./CSFFromAtropos5 ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf.nii.gz
	#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_not_3rd_ventricle.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf -odt char
	
	#./label_subject_all_steps_reconstruct_ventricles.py ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_lateral_ventricles.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf.nii.gz ${OUTPUTPREFIX}_segmentation_latvent.nii.gz
	#./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_segmentation_latvent.nii.gz --match 4 43 31 63 --noverbose 
	#ARGEVENTRICLEPART="-m Demons[$TEMPLATEDIR/FinaltemplateDKTMajorityLateralVentricles.nii.gz,${OUTPUTPREFIX}_lateral_ventricles_mask.nii.gz,2]"
	
	#MAXINTENSITY=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn -R | awk '{ print $2 }'`
	#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation -thr 2 -bin -mul -1 -add 1 -mul 1000 -add ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad.nii.gz
	#fslmaths ${OUTPUTPREFIX}_brain_mask -mul -1 -add 1 -mul 1000 -add ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad.nii.gz

	#./GaussianLaplacian -s 1 ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian1.nii.gz
	#./GaussianLaplacian -s 3 ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian3.nii.gz
#if [ ! -f "${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2.nii.gz" ]
#then
#	./GaussianLaplacian -s 2 ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2.nii.gz
#	fi
#	./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2_pos.nii.gz --min 0 --noverbose 
#	./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2_neg.nii.gz --max 0 --noverbose 
	

#	./MRIBinarize --i ${OUTPUTPREFIX}_brain_mask.nii.gz --o ${OUTPUTPREFIX}_brain_mask_eroded.nii.gz --match 1 --erode 10 --noverbose
#	fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2_neg.nii.gz -mas ${OUTPUTPREFIX}_brain_mask_eroded.nii.gz ${OUTPUTPREFIX}_brain_mask_eroded_to_remove.nii.gz
	
	#./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings.nii.gz --match 997 998 --dilate 15 --noverbose
	#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian2_neg.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings.nii.gz ${OUTPUTPREFIX}_rings_dilated_to_remove.nii.gz
	#./GaussianLaplacian -s 1 ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_laplacian1.nii.gz
	#${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz
	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian1.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian1_neg.nii.gz --max 0
	#ri_binarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian1.nii.gz --o ${OUTPUTPREFIX}_segmentation_gm_init.nii.gz --min 0 --mask ${OUTPUTPREFIX}_brain_mask.nii.gz --noverbose 
	./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz --o ${OUTPUTPREFIX}_segmentation_gm_init.nii.gz --match 1 --noverbose 

	#ThresholdImage 3 ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_kmeans4.nii.gz Kmeans 4 ${OUTPUTPREFIX}_brain_mask.nii.gz
	#ThresholdImage 3 ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_kmeans5.nii.gz Kmeans 5 ${OUTPUTPREFIX}_brain_mask.nii.gz

	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_kmeans5.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_kmeans5_dark_gm_remove.nii.gz --match 1 --noverbose
	# remove lateral ventricles
	#ImageMath 3 ${OUTPUTPREFIX}_segmentation_latvent_dilated.nii.gz MD ${OUTPUTPREFIX}_segmentation_latvent.nii.gz 3

	#fslmaths ${OUTPUTPREFIX}_segmentation_latvent_dilated.nii.gz -mul -1 -add 1 -mul ${OUTPUTPREFIX}_segmentation_gm_init.nii.gz ${OUTPUTPREFIX}_segmentation_gm.nii.gz -odt char
	fslmaths ${OUTPUTPREFIX}_segmentation_gm_init.nii.gz ${OUTPUTPREFIX}_segmentation_gm.nii.gz -odt char
	#rm -f ${OUTPUTPREFIX}_segmentation_latvent_dilated.nii.gz ${OUTPUTPREFIX}_segmentation_gm_init.nii.gz
	#fslmaths ${OUTPUTPREFIX}_inter.nii.gz -mul -1 -add 1 -mul ${OUTPUTPREFIX}_segmentation_gm_init.nii.gz ${OUTPUTPREFIX}_segmentation_gm.nii.gz -odt char
	# remove near edge of brain
	#./MRIBinarize --i ${OUTPUTPREFIX}_brain_mask.nii.gz --o ${OUTPUTPREFIX}_brain_mask_at_border.nii.gz --match 1 --inv --dilate 3 --noverbose
	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation1.nii.gz --match 1 --dilate 3 --noverbose
	
	#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_border_dilated.nii.gz -mas ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf_dilated.nii.gz ${OUTPUTPREFIX}_inter.nii.gz
	#rm -f ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf_dilated.nii.gz

	#rm -f ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_border_dilated.nii.gz
	#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_border.nii.gz -mas ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_csf_dilated.nii.gz ${OUTPUTPREFIX}_border_remove_init.nii.gz
	#./MRIBinarize --i ${OUTPUTPREFIX}_border_remove_init.nii.gz --o ${OUTPUTPREFIX}_border_remove_init.nii.gz --match 1 --dilate 1 --noverbose
	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_border.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_border.nii.gz --match 1 --dilate 2 --noverbose
	
	#./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove.nii.gz --match 51 52 12 13 9 48 54 17 53 18 --inv --erode 2 --noverbose
	
	# remove subcortical
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove.nii.gz --match 28 60 170 11 50 51 52 12 13 9 48 54 18 --dilate 7 --erode 3 --noverbose 
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_hippo.nii.gz --match 17 53 --noverbose 
	fslmaths ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove -add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_hippo.nii.gz -bin ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove -odt char
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz --match 91 93 --dilate 3 --erode 3 --noverbose  
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_dilated.nii.gz --match 91 93 --dilate 4 --noverbose  
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_skull_dilated.nii.gz --match 165 --dilate 4 --noverbose  
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz --match 24 --noverbose  
	fslmaths ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_skull_dilated.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_dilated.nii.gz -add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz -bin ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz -odt char
	fslmaths ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_dilated.nii.gz -add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz -bin ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz -odt char

	#rm -f ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_skull_dilated.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_dilated.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz
	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian3.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian3_neg.nii.gz --max 0 --noverbose 
	
	fslmaths ${OUTPUTPREFIX}_t2w_dark_wm_prob_skullstrip_reg.nii.gz -thr 0.3 -bin -mul -1 -add 1 -mas ${OUTPUTPREFIX}_segmentation_gm ${OUTPUTPREFIX}_segmentation_gm -odt char
	#rm -f ${OUTPUTPREFIX}_t2w_dark_wm_prob_skullstrip_reg.nii.gz
	#fslmaths ${OUTPUTPREFIX}_segmentation_gm ${OUTPUTPREFIX}_segmentation_gm_before_last -odt char
	
	./MRIBinarize --i ${OUTPUTPREFIX}_brain_mask.nii.gz --o ${OUTPUTPREFIX}_brain_mask_border.nii.gz --match 0 --dilate 5 --noverbose 

	./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz --o ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation_csf.nii.gz --match 3 --dilate 5 --noverbose 
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_skull_dilated.nii.gz --match 165 --dilate 7 --noverbose  
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_dilated.nii.gz --match 170 --dilate 7 --noverbose  
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz --match 24 --noverbose
	fslmaths ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_skull_dilated.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_dilated.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_to_remove.nii.gz
	fslmaths ${OUTPUTPREFIX}_brain_mask_border.nii.gz -mul ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation_csf.nii.gz -bin ${OUTPUTPREFIX}_brain_border_to_remove.nii.gz -odt char

	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz --match 2 41 --noverbose
	
	./MakeBrightVentricleMask $SUBJID $TISSUESEGDIR
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings.nii.gz --match 997 998 --dilate 3 --noverbose
	./MRIBinarize --i ${OUTPUTPREFIX}_brightmask_kmeans_class3.nii.gz --o ${OUTPUTPREFIX}_brightmask_kmeans_class3_dilated.nii.gz --match 1 --dilate 2 --noverbose
	fslmaths ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz -mas ${OUTPUTPREFIX}_brightmask_kmeans_class3.nii.gz ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings_to_remove.nii.gz
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_to_remove.nii.gz --match 31 63 --dilate 3 
	
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_isthmus_bigdilate.nii.gz --match 1010 2010 --dilate 10
	#$fslmaths ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_isthmus_to_remove.nii.gz -mas ${OUTPUTPREFIX}_brightmask_kmeans_class3_dilated.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_isthmus_to_remove.nii.gz

	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_bigdilate.nii.gz --match 31 63 --dilate 15
	./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_mask.nii.gz --match 2 41 31 63 4 43

	fslmaths ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_isthmus_bigdilate.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_bigdilate.nii.gz -mas ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_mask.nii.gz ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_isthmus_to_remove.nii.gz

	#cp ${OUTPUTPREFIX}_segmentation_gm.nii.gz ${OUTPUTPREFIX}_segmentation_gm1.nii.gz
	
	# remove GM between CSF and brain border
	./MRIBinarize --i ${OUTPUTPREFIX}_brain_mask.nii.gz --o ${OUTPUTPREFIX}_brain_mask_inv_dilated.nii.gz --match 0 --dilate 4

	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation.nii.gz --o  --match 0 --dilate 4

	fslmaths \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove.nii.gz \
		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz \
		-add ${OUTPUTPREFIX}_brain_border_to_remove.nii.gz \
		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_to_remove.nii.gz \
		-add ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings_to_remove.nii.gz \
		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_to_remove.nii.gz \
		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_isthmus_to_remove.nii.gz \
		-bin -mul -1 -add 1 -mul ${OUTPUTPREFIX}_segmentation_gm ${OUTPUTPREFIX}_segmentation_gm -odt char
		#-add ${OUTPUTPREFIX}_brain_mask_eroded_to_remove.nii.gz \
#	fslmaths \
#		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_brain_border_to_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_to_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings_to_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_to_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_brain_mask_eroded_to_remove.nii.gz \
#		-add ${OUTPUTPREFIX}_rings_dilated_to_remove.nii.gz \
#		-bin -mul -1 -add 1 -mul ${OUTPUTPREFIX}_segmentation_gm1 ${OUTPUTPREFIX}_segmentation_gm1 -odt char
#	
	# remove cerebellum and WM
	#ThresholdImage 3 ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_segmentation_gm_otsu2.nii.gz Otsu 2 ${OUTPUTPREFIX}_segmentation_gm.nii.gz
	cp ${OUTPUTPREFIX}_segmentation_gm.nii.gz ${OUTPUTPREFIX}_segmentation_gm_no_dark_wm_remove.nii.gz
	./RemoveDarkWMFromGMSegmentation $SUBJID
	#fslmaths ${OUTPUTPREFIX}_dark_wm_bright_gm_to_remove -mul -1 -add 1 -mas ${OUTPUTPREFIX}_segmentation_gm ${OUTPUTPREFIX}_segmentation_gm -odt char
	./ComponentAreaFilter -a 500 ${OUTPUTPREFIX}_segmentation_gm.nii.gz ${OUTPUTPREFIX}_segmentation_gm.nii.gz
	./ComponentAreaFilter -a 500 ${OUTPUTPREFIX}_segmentation_gm_no_dark_wm_remove.nii.gz ${OUTPUTPREFIX}_segmentation_gm_old_no_dark_wm_remove.nii.gz
	#fslmaths ${OUTPUTPREFIX}_segmentation_gm.nii.gz ${OUTPUTPREFIX}_segmentation_gm_lapinit.nii.gz
	rm -f \
		${OUTPUTPREFIX}_brain_border_to_remove.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos3_priors_segmentation_csf.nii.gz \
		${OUTPUTPREFIX}_brain_mask_border.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_kmeans5_dark_gm_remove.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainsteam_to_remove.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainsteam_dilated.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_skull_dilated.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings_to_remove.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_with_latvent_rings_skullstrip_reg_rings.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_choroid_to_remove.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_dilated.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_brainstem_to_remove.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_dilated.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_cerebellum_remove.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_hippo.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_lateral_ventricles.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_subcort_remove.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_brightpad_laplacian1.nii.gz \
		${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_not_3rd_ventricle.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_45.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_3.nii.gz \
		${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation_border.nii.gz
	
	./SeparateGMSegHemis $SUBJID
	#####
	#fslmaths ${OUTPUTPREFIX}_segmentation_gm.nii.gz ${OUTPUTPREFIX}_segmentation_gm_atroposinit.nii.gz
	#./MRIBinarize --i ${OUTPUTPREFIX}_t2w_restore_brain_dn_atropos5_segmentation.nii.gz --o ${OUTPUTPREFIX}_segmentation_csf.nii.gz --match 4 5 --noverbose 
	#./MRIBinarize --i ${OUTPUTPREFIX}_segmentation_latvent.nii.gz --o ${OUTPUTPREFIX}_segmentation_latvent_dilated_tmp.nii.gz --match 1 --dilate 2 --noverbose
	#fslmaths ${OUTPUTPREFIX}_segmentation_latvent_dilated_tmp.nii.gz -mul -1 -add 1 -mul ${OUTPUTPREFIX}_segmentation_csf.nii.gz ${OUTPUTPREFIX}_segmentation_csf.nii.gz -odt char	
	#rm -f ${OUTPUTPREFIX}_segmentation_latvent_dilated_tmp.nii.gz
	# cleanup
				#	${OUTPUTPREFIX}_t2w_dark_wm_prob_skullstrip_reg.nii.gz \
		#fslmaths ${OUTPUTPREFIX}_brain_mask_at_border.nii.gz -mul -1 -add 1 -mul ${OUTPUTPREFIX}_segmentation_gm ${OUTPUTPREFIX}_segmentation_gm -odt char
#	ls $TEMPLATEDIR/FinaltemplateDKTMajorityLateralVentriclesWithChoroid.nii.gz
#	antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 \
#		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_reg1Warp.nii.gz \
#		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat \
#		--transform BSplineSyN[0.1,26,0,3] \
#		--metric MI[ $T2TEMPLATEBRAIN,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.7,32 ] \
#		--metric Demons[ $TEMPLATEDIR/FinaltemplateDKTMajorityLateralVentriclesWithChoroid.nii.gz,${OUTPUTPREFIX}_segmentation_latvent.nii.gz,0.5 ] \
#		--metric Demons[ $TEMPLATEDIR/FinaltemplateAffineDKTTissueProbGM.nii.gz,${OUTPUTPREFIX}_segmentation_gm.nii.gz,0.5 ] \
#		--convergence [ 200x100x50x0,1e-6,10 ] --shrink-factors 6x4x2x1 --smoothing-sigmas 6x4x2x1vox \
#		--output [${OUTPUTPREFIX}_gmvent_reg,${OUTPUTPREFIX}_gmvent_regWarped.nii.gz]
	#if [ ! -f "${OUTPUTPREFIX}_gmvent_reg2Warp.nii.gz" ]
	#then
	#	antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 \
	#		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_reg1Warp.nii.gz \
	#		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat \
	#		--transform SyN[0.1,1.5,0] \
	#		--metric MI[ $T2TEMPLATEBRAIN,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.1,32 ] \
	#		--metric Demons[ $TEMPLATEDIR/FinaltemplateDKTMajorityLateralVentriclesWithChoroid.nii.gz,${OUTPUTPREFIX}_segmentation_latvent.nii.gz,0.5 ] \
	#		--metric Demons[ $TEMPLATEDIR/FinaltemplateAffineDKTTissueProbGM.nii.gz,${OUTPUTPREFIX}_segmentation_gm.nii.gz,0.5 ] \
	#		--convergence [ 200x100x100x5,1e-6,50 ] --shrink-factors 6x4x2x1 --smoothing-sigmas 6x4x2x1vox --masks [$TEMPLATEDIR/FinaltemplateBrainMask.nii.gz,${OUTPUTPREFIX}_brain_mask.nii.gz] \
	#		--output [${OUTPUTPREFIX}_gmvent_reg_syn,${OUTPUTPREFIX}_gmvent_reg_synWarped.nii.gz]
	#fi
		#fslmerge -a tmp Stage1_level*
	#rm -f Stage*.nii.gz
#fi
###

#T=`mktemp`
#rm -f $T
#for i in `seq -w 1 10`
#do#
#	#echo "./label_subject_all_steps_tissuesegs_to_training_one.sh $i" >> $T
#done

#export NUMTHREADS=`expr $NUMPROC / 5`
if [ "$NUMTHREADS" == "0" ]
then
	NUMTHREADS=1
fi
NUMTHREADS=2
NUMTHREADS=`expr $(nproc) / 5 + 1`
#NUMTHREADS=2

if [ "$DOPARALLEL" == "YES" ]
then
	export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NUMTHREADS
fi
#export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=`nproc`
#parallel -j10 --ungroup < $T
#rm -f $T
#if [ ! -f "${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz" ]
#then
	T=`tmpnam`
	rm -f $T
	for i in `seq -w 1 10`
	do
	#echo "./label_subject_all_steps_tissuesegs_transfornonly_one_highbright.sh $i" >> $T
		#./label_subject_all_steps_tissuesegs_transfornonly_one.sh $i
		#./label_subject_all_steps_tissuesegs_transformonly_syn.sh $i
		#./label_subject_all_steps_tissuesegs_transformonly_bspline.sh $i
		#echo "./label_subject_all_steps_tissuesegs_transformonly_bspline.sh $i" >> $T
	#	#echo "./label_subject_all_steps_tissuesegs_transformonly_syn.sh $i" >> $T
		#echo "./label_subject_all_steps_tissuesegs_transformonly_syn_noring.sh $i" >> $T
		#echo "./label_subject_all_steps_tissuesegs_transformonly_bspline_noring.sh $i" >> $T
		if [ "$DOPARALLEL" == "YES" ]
		then
			echo "./label_subject_all_steps_tissuesegs_transformonly_bspline_noring_gmsep.sh $i" >> $T
		else
			./label_subject_all_steps_tissuesegs_transformonly_bspline_noring_gmsep.sh $i
		fi
#./label_subject_all_steps_tissuesegs_transformonly_bspline_noring.sh $i
		#exit

	done

	if [ "$DOPARALLEL" == "YES" ]
	then
		parallel -j5 --ungroup < $T
		export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=`nproc`
	fi
	#parallel -j1 --ungroup < $T
	rm -f $T
for i in `seq -w 1 10`
do
#	if [ ! -f "${OUTPUTPREFIX}_to_P${i}_reg2InverseWarp.nii.gz" ]
#	then
#		antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 \
#			--initial-moving-transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,1] \
#			--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat \
#			--transform SyN[0.1,0.5,0] \
#			--metric Demons[ $TEMPLATEDIR/P${i}_gm.nii.gz,${OUTPUTPREFIX}_segmentation_gm.nii.gz,0.8 ] \
#			--metric MI[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.1,32 ] \
#			--metric Demons[ $TEMPLATEDIR/P${i}_dkt_latvent.nii.gz,${OUTPUTPREFIX}_segmentation_latvent.nii.gz,0.1 ] \
#			--convergence [ 200x200x200x100x80,1e-6,50 ] --shrink-factors 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
#			--output ${OUTPUTPREFIX}_to_P${i}_reg
#			#--output [${OUTPUTPREFIX}_to_P${i}_reg,${OUTPUTPREFIX}_to_P${i}_regWarped.nii.gz]
#			#--convergence [ 200x100x100x0,1e-6,50 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 6x4x2x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
#			#--metric Demons[ $TEMPLATEDIR/P${i}_csf.nii.gz,$P/${OUTPUTPREFIX}_segmentation_csf.nii.gz,0.3 ] \
#	fi
	#rm -f ${OUTPUTPREFIX}_to_P${i}_reg2Warp.nii.gz
	
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_t2.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz  \
#		--interpolation Linear \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_$SUBJID.nii.gz
#			
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_dkt_with_skull_label.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_$SUBJID.nii.gz
#	
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_dkt_with_latvent_rings_gm.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_latvent_rings_gm_to_$SUBJID.nii.gz
#	
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_vent_and_centre_bright.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_vent_and_centre_bright_to_$SUBJID.nii.gz
#		
	ATLAST2IMAGESNAMES="$ATLAST2IMAGESNAMES $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_$SUBJID.nii.gz"
	ATLAST2IMAGES="$ATLAST2IMAGES -g $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_$SUBJID.nii.gz"
	ATLASDKTIMAGES="$ATLASDKTIMAGES -l $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_$SUBJID.nii.gz"
	ATLASDKTIMAGESNAMES="$ATLASDKTIMAGESNAMES $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_$SUBJID.nii.gz"
	fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_${SUBJID}
	fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_${SUBJID}

#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_t2.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_reg2InverseWarp.nii.gz  \
#		--interpolation Linear \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_t2_to_$SUBJID.nii.gz
#			
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_dkt_with_skull_label.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_reg2InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_with_skull_label_to_$SUBJID.nii.gz
	#ATLAST2IMAGESNAMES="$ATLAST2IMAGESNAMES $TISSUESEGDIR/$SUBJID/P${i}_t2_to_$SUBJID.nii.gz"
	#ATLAST2IMAGES="$ATLAST2IMAGES -g $TISSUESEGDIR/$SUBJID/P${i}_t2_to_$SUBJID.nii.gz"
	#ATLASDKTIMAGES="$ATLASDKTIMAGES -l $TISSUESEGDIR/$SUBJID/P${i}_dkt_with_skull_label_to_$SUBJID.nii.gz"
	#ATLASDKTIMAGESNAMES="$ATLASDKTIMAGESNAMES $TISSUESEGDIR/$SUBJID/P${i}_dkt_with_skull_label_to_$SUBJID.nii.gz"
	#fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain $TISSUESEGDIR/$SUBJID/P${i}_t2_to_${SUBJID}
	#fslcpgeom ${OUTPUTPREFIX}_t2w_restore_brain $TISSUESEGDIR/$SUBJID/P${i}_dkt_with_skull_label_to_${SUBJID}
done
fslmerge -a $TISSUESEGDIR/$SUBJID/all_dkt_antsinit_to_${SUBJID} $ATLASDKTIMAGESNAMES
fslmerge -a $TISSUESEGDIR/$SUBJID/all_t2_antsinit_to_${SUBJID} $ATLAST2IMAGESNAMES
#fslmerge -a $TISSUESEGDIR/$SUBJID/all_dkt_antsinit_with_vent_and_centre_to_${SUBJID} $TISSUESEGDIR/$SUBJID/P[0-9][0-9]_dkt_antsinit_with_vent_and_centre_bright_to_$SUBJID.nii.gz
#fslmerge -a $TISSUESEGDIR/$SUBJID/all_antsinit_with_latvent_rings_gm_to_${SUBJID} $TISSUESEGDIR/$SUBJID/P[0-9][0-9]_dkt_antsinit_with_latvent_rings_gm_to_$SUBJID.nii.gz
#fslmerge -a $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID} $ATLASDKTIMAGESNAMES
#fslmerge -a $TISSUESEGDIR/$SUBJID/all_t2_to_${SUBJID} $ATLAST2IMAGESNAMES
#fslmaths $TISSUESEGDIR/$SUBJID/all_t2_to_${SUBJID} $TISSUESEGDIR/$SUBJID/all_t2_to_${SUBJID} -odt short
#if [ ! -f "${OUTPUTPREFIX}_labelfusionimage_dkt2.nii.gz" ]
#then
#	antsJointFusion $ATLAST2IMAGES $ATLASDKTIMAGES --target-image ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -x ${OUTPUTPREFIX}_brain_mask.nii.gz --output ${OUTPUTPREFIX}_labelfusionimage_dkt2.nii.gz -v --patch-metric PC -s 2 -p 2
#fi
#fslmaths ${OUTPUTPREFIX}_labelfusionimage_newmask.nii.gz ${OUTPUTPREFIX}_labelfusionimage_newmask.nii.gz -odt short
#fslmaths ${OUTPUTPREFIX}_labelfusionimage_dkt2.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt2.nii.gz -odt short

#if [ ! -f "${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz" ]
#then
	antsJointFusion $ATLAST2IMAGES $ATLASDKTIMAGES --target-image ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -x ${OUTPUTPREFIX}_brain_mask.nii.gz --output ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz -v --patch-metric PC -s 2 -p 2
#fi
#fslmaths ${OUTPUTPREFIX}_labelfusionimage_newmask.nii.gz ${OUTPUTPREFIX}_labelfusionimage_newmask.nii.gz -odt short
#fslmaths ${OUTPUTPREFIX}_labelfusionimage_dkt2.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt2.nii.gz -odt short
fslmaths ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz -odt short

mri_binarize --i ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz --o ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_lh_hippo_closed.nii.gz --match 17 --dilate 2 --erode 2 --binval 17 --noverbose
mri_binarize --i ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz --o ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_rh_hippo_closed.nii.gz --match 53 --dilate 2 --erode 2 --binval 53 --noverbose

mri_mask -transfer 17 ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_lh_hippo_closed.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz
mri_mask -transfer 53 ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_rh_hippo_closed.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz
rm -f ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_lh_hippo_closed.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_rh_hippo_closed.nii.gz

mri_binarize --i ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz --o ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_notgm.nii.gz --match $GMMATCH --inv --noverbose

./MRIBinarize --i ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz --o ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_lh_gm_dilated.nii.gz --min 1000 --max 1036 --noverbose --dilate 5
./MRIBinarize --i ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz --o ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_rh_gm_dilated.nii.gz --min 2000 --max 2036 --noverbose --dilate 5

./MRIBinarize --i ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig.nii.gz --o ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_wm_dilated.nii.gz --match 2 41 --noverbose --dilate 1
fslmaths ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_lh_gm_dilated.nii.gz -mas ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_rh_gm_dilated.nii.gz -mas ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_notgm.nii.gz -mas ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_wm_dilated.nii.gz -bin -mul 24 ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_to_csf.nii.gz
mri_mask -transfer 24 ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_orig_to_csf.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz

cp ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit.nii.gz ${OUTPUTPREFIX}_labelfusionimage_dkt_antsinit_edited.nii.gz
#rm -f ${TISSUESEGDIR}/${SUBJID}/*reg6* ${TISSUESEGDIR}/${SUBJID}/P*

# ./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz --match 2 41 --noverbose
#./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_gm.nii.gz --min 1000 --noverbose
#./MRIBinarize --i ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg.nii.gz --o ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz --match 4 43 24 --noverbose
#
#CSFMEAN=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -k ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_csf.nii.gz -m`
#GMMEAN=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -k ${OUTPUTPREFIX}_segmentation_gm.nii.gz -m`
#WMMEAN=`fslstats ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -k ${OUTPUTPREFIX}_majority_dkt_skullstrip_reg_wm.nii.gz -m`
#
##echo $CSFMEAN $GMMEAN $WMMEAN
#
#
##if [ ! -f "${OUTPUTPREFIX}_t2w_restore_brain_dn_majority_gm_segmentation.nii.gz" ]
##then
#    Atropos -a ${OUTPUTPREFIX}_t2w_restore_brain_dn.nii.gz -d 3 -o ${OUTPUTPREFIX}_t2w_restore_brain_dn_majority_gm_segmentation.nii.gz -x ${OUTPUTPREFIX}_brain_mask.nii.gz --winsorize-outliers BoxPlot[0.25,0.75,1.5] -i kmeans[ 4,0x${GMMEAN}x${WMMEAN}x${CSFMEAN} ] -c [ 3,0.0 ] -k Gaussian -m [ 0.1,1x1x1 ] -r 1 --verbose 1
##fi
#
#fslmaths ${OUTPUTPREFIX}_t2w_restore_brain_dn_majority_gm_segmentation.nii.gz ${OUTPUTPREFIX}_t2w_restore_brain_dn_majority_gm_segmentation.nii.gz -odt char
#
#exit
