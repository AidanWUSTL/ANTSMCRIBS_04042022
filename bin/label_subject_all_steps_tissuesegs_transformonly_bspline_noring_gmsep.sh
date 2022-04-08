#!/bin/bash

# AL changes:
# removed pwd var P from path for files, since 
# a change in the previous labelling script uses teh OUTPUTPREFIX var

i=$1

#if [ ! -f "${OUTPUTPREFIX}_to_P${i}_reg2InverseWarp.nii.gz" ]
#then
#	T=`mktemp -d`
#	P=`pwd`
#	cd $T
#
#	#rm -f Stage*.nii.gz
#	#antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 \
#	antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 --write-interval-volumes 10 \
#		--initial-moving-transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,1] \
#		--initial-moving-transform $P/${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat \
#		--transform SyN[0.1,0.5,0] \
#		--metric MI[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,$P/${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.1,32 ] \
#		--metric Demons[ $TEMPLATEDIR/P${i}_dkt_with_latvent_rings_gm.nii.gz,${P}/${OUTPUTPREFIX}_segmentation_gm.nii.gz,0.45 ] \
#		--metric Demons[ $TEMPLATEDIR/P${i}_vent_and_centre_bright.nii.gz,${P}/${OUTPUTPREFIX}_brightmask_kmeans_class3.nii.gz,0.45 ] \
#		--convergence [ 200x200x200x100x50,1e-6,50 ] --shrink-factors 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
#		--output ${P}/${OUTPUTPREFIX}_to_P${i}_reg
#		#--output [${OUTPUTPREFIX}_to_P${i}_reg,${OUTPUTPREFIX}_to_P${i}_regWarped.nii.gz]
#		#--convergence [ 200x100x100x0,1e-6,50 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 6x4x2x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
#		#--metric Demons[ $TEMPLATEDIR/P${i}_csf.nii.gz,$P/${OUTPUTPREFIX}_segmentation_csf.nii.gz,0.3 ] \
#	fslmerge -a ${P}/${OUTPUTPREFIX}_to_P${i}_gm_all_iters Stage*.nii.gz
#	#rm -f Stage*.nii.gz
#	cd $P
#	rm -fr $T
#
P=`pwd`
if [ ! -f "${OUTPUTPREFIX}_to_P${i}_antsinit_reg6InverseWarp.nii.gz" ]
then
	T=`mktemp -d`
	cd $T

	#rm -f Stage*.nii.gz
	#antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 \
#ImageMath 3 ${P}/${OUTPUTPREFIX}_segmentation_gm_dt.nii.gz D ${P}/${OUTPUTPREFIX}_segmentation_gm.nii.gz
	#antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 --write-interval-volumes 10 \

	antsRegistration -v -d 3 -u 1 -w [ 0.025,0.975 ] --verbose 1 --float 1 --collapse-output-transforms 0 \
		--initial-moving-transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,1] \
		--initial-moving-transform $TEMPLATEDIR/FinalP${i}1InverseWarp.nii.gz \
		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_reg1Warp.nii.gz \
		--initial-moving-transform ${OUTPUTPREFIX}_skullstrip_reg0GenericAffine.mat \
		--transform BSplineSyN[0.4,26,0,3] \
		--metric MI[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.05,32 ] \
		--metric CC[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.5,2 ] \
		--metric Demons[ $TEMPLATEDIR/P${i}_dkt_gm_lh.nii.gz,${OUTPUTPREFIX}_segmentation_gm_sep_lh.nii.gz,1 ] \
		--metric Demons[ $TEMPLATEDIR/P${i}_dkt_gm_rh.nii.gz,${OUTPUTPREFIX}_segmentation_gm_sep_rh.nii.gz,1 ] \
		--metric Demons[ $TEMPLATEDIR/P${i}_vent_and_centre_bright.nii.gz,${OUTPUTPREFIX}_brightmask_kmeans_class3.nii.gz,0.1 ] \
		--convergence [ 200x200x50x0x0,1e-6,75 ] --shrink-factors 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
		--transform SyN[0.1,0.5,0] \
		--metric MI[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.01,128 ] \
		--metric CC[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.5,2 ] \
		--metric Demons[ $TEMPLATEDIR/P${i}_dkt_gm_lh.nii.gz,${OUTPUTPREFIX}_segmentation_gm_sep_lh.nii.gz,1 ] \
		--metric Demons[ $TEMPLATEDIR/P${i}_dkt_gm_rh.nii.gz,${OUTPUTPREFIX}_segmentation_gm_sep_rh.nii.gz,1 ] \
		--metric Demons[ $TEMPLATEDIR/P${i}_vent_and_centre_bright.nii.gz,${OUTPUTPREFIX}_brightmask_kmeans_class3.nii.gz,0.1 ] \
		--convergence [ 200x100x50x0,1e-6,75 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 4x2x1x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
		--output ${OUTPUTPREFIX}_to_P${i}_antsinit_reg
		#--convergence [ 200x200x200x100x0,1e-6,50 ] --shrink-factors 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
		#--output [${OUTPUTPREFIX}_to_P${i}_reg,${OUTPUTPREFIX}_to_P${i}_regWarped.nii.gz]
		#--transform SyN[0.1,0.5,0] \
		#--convergence [ 200x100x100x0,1e-6,50 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 6x4x2x1vox --masks $TEMPLATEDIR/P${i}_reg_mask.nii.gz \
		#--metric MI[ $TEMPLATEDIR/P${i}_t2_brain.nii.gz,$P/${OUTPUTPREFIX}_t2w_restore_brain.nii.gz,0.1,32 ] \
		#--metric Demons[ $TEMPLATEDIR/P${i}_dkt_with_latvent_rings_gm.nii.gz,${P}/${OUTPUTPREFIX}_segmentation_gm.nii.gz,0.45 ] \
		#--metric Demons[ $TEMPLATEDIR/P${i}_csf.nii.gz,$P/${OUTPUTPREFIX}_segmentation_csf.nii.gz,0.3 ] \
	#fslmerge -a ${P}/${OUTPUTPREFIX}_to_P${i}_gm_antsinit_all_iters Stage*.nii.gz
	#rm -f Stage*.nii.gz
	cd $P
	rm -fr $T
fi	
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_t2.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz  \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5InverseWarp.nii.gz  \
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
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5InverseWarp.nii.gz  \
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
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5InverseWarp.nii.gz  \
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
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_vent_and_centre_bright_to_$SUBJID.nii.gz
#		
	# composite transform	
	#-------------------------------------------------------------	
	if [ -f "${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz" ]
	then	
		antsApplyTransforms -v -d 3 --reference-image $TEMPLATEDIR/P${i}_t2.nii.gz \
			--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz  \
			--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5InverseWarp.nii.gz  \
			--output [${OUTPUTPREFIX}_to_P${i}_antsinit_reg6InverseWarp.nii.gz,1]
	fi
	rm -f ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4InverseWarp.nii.gz ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5InverseWarp.nii.gz 
	rm -f ${OUTPUTPREFIX}_to_P${i}_antsinit_reg4Warp.nii.gz ${OUTPUTPREFIX}_to_P${i}_antsinit_reg5Warp.nii.gz 
	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_t2.nii.gz \
		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg6InverseWarp.nii.gz  \
		--interpolation Linear \
		--output $TISSUESEGDIR/$SUBJID/P${i}_t2_antsinit_to_$SUBJID.nii.gz
			
	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_dkt_with_skull_label.nii.gz \
		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg6InverseWarp.nii.gz  \
		--interpolation GenericLabel \
		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_skull_label_to_$SUBJID.nii.gz
	
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_dkt_with_latvent_rings_gm.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg6InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_latvent_rings_gm_to_$SUBJID.nii.gz
#	
#	antsApplyTransforms -v -d 3 --reference-image ${OUTPUTPREFIX}_t2w_restore.nii.gz --input $TEMPLATEDIR/P${i}_vent_and_centre_bright.nii.gz \
#		--transform [${OUTPUTPREFIX}_skullstrip_affine0GenericAffine.mat,1] \
#		--transform ${OUTPUTPREFIX}_skullstrip_reg1InverseWarp.nii.gz \
#		--transform $TEMPLATEDIR/FinalP${i}1Warp.nii.gz \
#		--transform [$TEMPLATEDIR/FinalP${i}0GenericAffine.mat,0] \
#		--transform ${OUTPUTPREFIX}_to_P${i}_antsinit_reg6InverseWarp.nii.gz  \
#		--interpolation GenericLabel \
#		--output-data-type short \
#		--output $TISSUESEGDIR/$SUBJID/P${i}_dkt_antsinit_with_vent_and_centre_bright_to_$SUBJID.nii.gz
#		
	
#fi
