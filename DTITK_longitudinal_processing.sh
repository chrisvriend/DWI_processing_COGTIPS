#!/bin/bash
set -euo pipefail

# Written by Anders L. Thorsen & C. Vriend April 2020

# This script will register two longitudinal DTI volumes to a common within-subject template.
# The within-subject template can then be used to register mutliple subjects to each other.
# A key benefit of this approach is that you minimize the number of transformations of the data,
# rather than adding extra interpolations when tranforming between time points (within subject) and between subjects.
# URL of manual is available here (http://dti-tk.sourceforge.net/pmwiki/pmwiki.php?n=Documentation.HomePage)
# and useful online discussion with code is available here (https://groups.google.com/forum/#!topic/dtitk/ciF3INwnGt8)

#########################################
# Setup relevant software and variables
#########################################

# Sets up DTI-TK and FSL
dtitkdir=/path/to/dtitk-2.3.1
export PATH=${dtitkdir}/bin:${PATH}
export PATH=${dtitkdir}/utilities:${PATH}
export PATH=${dtitkdir}/scripts:${PATH}
export DTITK_ROOT=${dtitkdir}
echo "DTI-TK sourced"

# source FSL

# = server and version specific

###

# Sets up variables for folder with tensor images from all subjects and recommended template from DTI-TK
tensor_dir=/path/to/inputdir
QCdir=${tensor_dir}/QC
# location of tract atlas for extracting tract-specific diffusion metrics
atlastractsdir=${headdir}/atlastracts

mkdir -p ${QCdir}


# IXI aging template
# https://www.nitrc.org/frs/shownotes.php?release_id=2273
templatedir=/path/to/ixi_aging_template_v3.0/template
# number of iterations to use for affine registration
Niter=3

# ASSUMED DATA STRUCTURE
# one folder for each subject containing the eddy corrected DWI data (data.nii.gz) and the accompanying bvals and bvecs files + a brain mask (nodif_brain_mask.nii.gz)
# each subject folder can contain more than one session (in this case T0 and T1) for longitudinal processing
# IMPORTANT NOTE: the data should only contain one shell (+ b0 volumes) for the tensor fit. Multishell data should therefore first be split.

cd ${tensor_dir}
#########################################
# dtifit
#########################################


echo "run dtifit"

for subj in `cat subjects.txt`; do

	cd ${thuis}/${subj}
	echo ${subj}
	# assumes two sessions; change if necessary
	for session in T0 T1; do

		if [ ! -f DWI_${subj}_${session}_FA.nii.gz ]; then

			dtifit --data=data_${session} --mask=nodif_brain_mask.nii.gz --bvecs=${session}.bvecs --bvals=${session}.bvals --out=DWI_${subj}_${session}

		fi

	done
done




#########################################
# Make dtifitâ€™s dti_V{123} and dti_L{123} compatible with DTI-TK
#########################################

for subj in `cat subjects.txt`; do

	cd ${thuis}/${subj}

	echo ${subj}
	for session in T0 T1; do


		if [ ! -f DWI_${subj}_${session}_dtitk.nii.gz ]; then

			${dtitkdir}/scripts/fsl_to_dtitk DWI_${subj}_${session}

		fi

	done
done


#########################################
# intra-subject registration
#########################################


for subj in `cat subjects.txt`; do

	cd ${tensor_dir}/${subj}

	#########################################
	# write txt file with *dtitk.nii.gz files
	#########################################

	if [ ! -f ${subj}.txt ]; then
		ls -1 *dtitk.nii.gz > ${subj}.txt
	fi

	#########################################
	# Performs initial construction of the subject-specific template
	#########################################

	if [ ! -f ${subj}_mean_initial.nii.gz ]; then
		echo "running intial template construction"
		dti_template_bootstrap ${templatedir}/ixi_aging_template.nii.gz ${subj}.txt EDS
		mv mean_initial.nii.gz ${subj}_mean_initial.nii.gz
	else
		echo "template bootstrapping has already been run"
	fi

	#########################################
	# Performs first affine (linear) registration of subject images to subject-specific template (named ${subj}_mean_affine3.nii.gz)
	#########################################

	if [ ! -f ${subj}_mean_affine3.nii.gz ]; then
		echo "running affine registration to initial template"
		dti_affine_population ${subj}_mean_initial.nii.gz ${subj}.txt EDS ${Niter}
		mv mean_affine3.nii.gz ${subj}_mean_affine3.nii.gz
	else
		echo "affine registration has already been run"
	fi

	# Creates a binary mask of the affine subject-specific template
	if [ ! -f ${subj}_mask.nii.gz ]; then
		echo "making binary mask for intial template construction"
		TVtool -in ${subj}_mean_affine3.nii.gz -tr
		BinaryThresholdImageFilter ${subj}_mean_affine3_tr.nii.gz ${subj}_mask.nii.gz 0.01 100 1 0
	else
		echo "binary mask has already been made"
	fi

	#########################################
	# Improves the subject-specific template
	# and creates deformation field,
	# stores aligned volumes and puts their filenames in "{subj_aff_diffeo.txt"
	#########################################

	if [ ! -f ${subj}_diffeomorphic.nii.gz ]; then
		echo "making diffeomorphic warps"
		dti_diffeomorphic_population ${subj}_mean_affine3.nii.gz ${subj}_aff.txt ${subj}_mask.nii.gz 0.002
		mv mean_diffeomorphic_initial6.nii.gz ${subj}_diffeomorphic.nii.gz
	else
		echo "diffeomorphic warps already made"
	fi


	##########################################
	# Creates non-linear transform from
	#individual timepoint to subject-specific template
	#########################################
	for session in T0 T1; do
		if [ ! -f DWI_${subj}_${session}_combined.df.nii.gz ]; then
			echo "Making non-inear transform for first timepoint"

			dfRightComposeAffine -aff DWI_${subj}_${session}_dtitk.aff -df DWI_${subj}_${session}_dtitk_aff_diffeo.df.nii.gz -out DWI_${subj}_${session}_combined.df.nii.gz

		else
			echo "Transform to subject-specific template already exists"
		fi

		#########################################
		# Warps individual timepoint to subject-specific template
		#########################################

		if [ ! -f ${subj}_${session}_combined.nii.gz ]; then
			echo "warping ${session} timepoint to subject-specific template"
			deformationSymTensor3DVolume -in DWI_${subj}_${session}_dtitk.nii.gz -trans DWI_${subj}_${session}_combined.df.nii.gz -target ${subj}_mean_initial.nii.gz -out ${subj}_${session}_combined.nii.gz
		else
			echo "Warped timepoint to subject-specific template already exists for ${session}"
		fi


	done
	#########################################
	# Calculates mean  image for both time points in subject-template space
	#########################################

	if [ ! -f ${subj}_mean_intra_template.nii.gz ]; then
		echo "creating mean image of time points in subject-specific template space"

		ls ${subj}_T*_combined.nii.gz > ${subj}_intra_reg_volumes.txt
		TVMean -in ${subj}_intra_reg_volumes.txt -out ${subj}_mean_intra_template.nii.gz
	else
		echo "mean image of  timepoints in subject-specific template space already exists"
	fi

done


#########################################
# QC
#########################################

echo "check registration by comparing the warped images in intra-subject template space (*combined.nii.gz) to the intra-subject template"

cd ${tensor_dir}

for subj in `cat subjects.txt`; do
	echo ${subj}
	cd ${subj}

	for session in T0 T1; do
		overlay 1 0 ${subj}_mean_affine3.nii.gz -a ${subj}_${session}_combined.nii.gz 2 5 ${subj}_intrareg_${session}.nii.gz
		overlay 1 0 ${subj}_mean_affine3.nii.gz -a ${subj}_${session}_combined.nii.gz 2 5 ${subj}_intrareg_${session}.nii.gz

		slicer ${subj}_intrareg_${session}.nii.gz -e 0 -L -t -a ${QCdir}/${subj}_${session}_overlay3D.png
		slicer ${subj}_intrareg_${session}.nii.gz -e 0 -L -t -A 2000 ${QCdir}/${subj}_${session}_overlayslices.png

		rm ${subj}_intrareg_${session}.nii.gz
	done
	cd ${tensor_dir}

done


#########################################
# START INTER SUBJECT REGISTRATION
#########################################
echo "start inter-subject registration pipeline"


mkdir -p ${tensor_dir}/interreg
cd ${tensor_dir}/interreg

# DTITK can't handle having these files in separate folders.
# these files therefore need to be combined into one folder
find .. -name "*_mean_intra_template.nii.gz" -maxdepth 2 > scans.txt

for scan in `cat scans.txt`; do
	base=$(echo ${scan} | sed 's:.*/::')
	# make symbolic link
	ln -s ${scan} ${base}
	# alternatively use move
	#mv ${scan} .

done

# important here to check that this went correctly

#########################################
# Creates intial group template from within-subject templates
#########################################
ls -1 *_mean_intra_template.nii.gz > inter_subjects.txt

if [ ! -f mean_initial.nii.gz ]; then
	echo "Running intial template construction"
	dti_template_bootstrap ${templatedir}/ixi_aging_template.nii.gz inter_subjects.txt EDS
else
	echo "Initial template already exists"
fi

#########################################
# Create affine (linear) warps for each subject to group template
#########################################

if [ ! -f mean_affine3.nii.gz ]; then
	echo "Running affine registration to initial template"
	dti_affine_population mean_initial.nii.gz inter_subjects.txt EDS ${Niter}
else
	echo "Affine registrations already exists"
fi

# make pngs of overlay with slicer for QC
fslroi mean_initial mean_initial_vslicer 0 1
for subj in $(ls *_aff.nii.gz); do
	echo ${subj}
	base=${subj%.nii.gz*}
	fslroi ${subj} aff 0 1
	slicer mean_initial_vslicer aff -a ${base}_overlay.png
	rm aff.nii.gz
done

# use dti_affine_sn for when to register individual subjects to existing template

#########################################
# Create mask of initial template to exclude voxels outside the brain
#########################################

if [ ! -f mask.nii.gz ]; then
	echo "Creating binary mask of initial template"
	TVtool -in mean_affine3.nii.gz -tr
	BinaryThresholdImageFilter mean_affine3_tr.nii.gz mask.nii.gz 0.01 100 1 0
else
	echo "Binary mask already exists"
fi


#########################################
# Improves the subject-specific template and creates deformation fields, stores aligned volumes and puts their filenames in "subj_aff_diffeo.txt"
#########################################

if [ ! -e "mean_diffeomorphic_initial6.nii.gz" ]; then
	echo "making diffeomorphic warps"
	dti_diffeomorphic_population mean_affine3.nii.gz inter_subjects_aff.txt mask.nii.gz 0.002
else
	echo "diffeomorphic warps already made"
fi

# make pngs of overlay with slicer for QC
fslroi mean_diffeomorphic_initial6.nii.gz mean_diffeomorphic_initial6_vslicer 0 1
for subj in $(ls *_aff_diffeo.nii.gz); do
	echo ${subj}
	base=${subj%_aff_diffeo.nii.gz*}_diffeo
	fslroi ${subj} diffeo 0 1
	slicer mean_diffeomorphic_initial6_vslicer diffeo -a ${base}_overlay.png
	rm diffeo.nii.gz
done

# when *intial6.nii.gz available but certain subject-specific warps have not yet been created
# ${dtitkdir}/scripts/dti_diffeomorphic_sn mean_final.nii.gz subjects_aff.txt mean_affine${Niter}_mask.nii.gz 6 0.002

#############################################
# Generate the spatially normalized DTI subject data with the isotropic 1mm3 resolution
#############################################
# Do for subjects with longitudinal data
warp_dir=${tensor_dir}interreg/warps
mkdir -p ${warp_dir}


echo "Warping images from native space to group-template space"
for subj in `cat ${tensor_dir}/subjects.txt`; do


	if [ ! -f ${warp_dir}/${subj}_inter_subject_combined.df.nii.gz ]; then
		# Combine affine and diffeomorphic warp fields from intra-subject space -> inter-subject template space
		dfRightComposeAffine -aff ${subj}_mean_intra_template.aff -df ${subj}_mean_intra_template_aff_diffeo.df.nii.gz -out ${warp_dir}/${subj}_inter_subject_combined.df.nii.gz
	fi

	for session in T0 T1; do
		if [ ! -f ${warp_dir}/${subj}_${session}_native_to_inter_subject_combined.df.nii.gz ]; then
			# Combine warp fields from T0 intra-subject space -> inter-subject template space
			dfComposition -df1 ${tensor_dir}/${subj}/DWI_${subj}_${session}_combined.df.nii.gz -df2 ${warp_dir}/${subj}_inter_subject_combined.df.nii.gz -out ${warp_dir}/${subj}_${session}_native_to_inter_subject_combined.df.nii.gz
		fi
		if [ ! -f ${warp_dir}/${subj}_${session}_2templatespace.dtitk.nii.gz ]; then
			# Warp image for T0 from native -> inter-subjecte space
			deformationSymTensor3DVolume -in ${tensor_dir}/${subj}/DWI_${subj}_${session}_dtitk.nii.gz -trans ${warp_dir}/${subj}_${session}_native_to_inter_subject_combined.df.nii.gz \
			-target ${tensor_dir}/interreg/mean_diffeomorphic_initial6.nii.gz -out ${warp_dir}/${subj}_${session}_2templatespace.dtitk.nii.gz -vsize 1 1 1
		fi


	done
done

# for subject with cross-sectional data
#dti_warp_to_template_group subjects.txt mean_final 1 1 1


#############################################
#Generate the population-specific DTI template with the isotropic 1mm3 spacing
#############################################
#subjs_normalized.txt is an ASCII text file that contains a list of the file names of the normalized high-resolution DTI volumes from the previous step

cd ${warp_dir}

ls -1 *2templatespace.dtitk.nii.gz > subjs_warped.txt

if [ ! -f ${warp_dir}/mean_final_high_res.nii.gz ]; then
	TVMean -in subjs_warped.txt -out mean_final_high_res.nii.gz
fi

#############################################
# Generate the FA map of the high-resolution population-specific DTI template
#############################################
if [ ! -f ${warp_dir}/mean_FA.nii.gz ]; then
	TVtool -in mean_final_high_res.nii.gz -fa

	#############################################
	#Rename the FA map to be consistent with the TBSS pipeline
	#############################################

	mv mean_final_high_res_fa.nii.gz mean_FA.nii.gz

fi

#############################################
#Generate the white matter skeleton from the high-resolution FA map of the DTI template
#############################################

if [ ! -f ${warp_dir}/mean_FA_skeleton.nii.gz ]; then
	tbss_skeleton -i mean_FA -o mean_FA_skeleton
fi

#############################################
# Generate the FA map of the spatially normalized high-resolution DTI data
#############################################

for subj in `cat subjs_warped.txt`; do
	TVtool -in ${subj} -fa
	TVtool -in ${subj} -ad
	TVtool -in ${subj} -rd
	# tr = 3 times md
	TVtool -in ${subj} -tr

done

fslmerge -t all_FA *2templatespace.dtitk_fa.nii.gz
fslmerge -t all_RD *2templatespace.dtitk_rd.nii.gz
fslmerge -t all_TR *2templatespace.dtitk_tr.nii.gz
fslmerge -t all_AD *2templatespace.dtitk_ad.nii.gz

# create 4D image of all MD maps
fslmaths all_TR -div 3 all_MD
# if you want to have subject/session specific maps you will need to split this file

#############################################
#apply fslmaths to all_FA to create a combined binary mask volume called mean_FA_mask
#############################################

fslmaths all_FA -max 0 -Tmin -bin mean_FA_mask -odt char
fslmaths all_AD -max 0 -Tmin -bin mean_AD_mask -odt char
fslmaths all_RD -max 0 -Tmin -bin mean_RD_mask -odt char
fslmaths all_MD -max 0 -Tmin -bin mean_MD_mask -odt char

fslmaths mean_FA_skeleton -mas mean_FA_mask mean_FA_skeleton_mskd

#############################################
#Place the TBSS relevant files into a folder that TBSS expects
#############################################

cd ${tensor_dir}
mkdir -p tbss
cd tbss
mkdir -p stats
cd ..
cp ${warp_dir}/mean_FA.nii.gz ${warp_dir}/mean_FA_skeleton_mskd.nii.gz ${warp_dir}/all_*.nii.gz ${warp_dir}/mean_FA_mask.nii.gz ${tensor_dir}/tbss/stats
cd ${tensor_dir}/tbss/stats
mv mean_FA_skeleton_mskd.nii.gz mean_FA_skeleton.nii.gz
cd ..

#############################################
# final tbss step before stats, including randomize
#############################################

thresh=0.2

# on FA image
tbss_4_prestats ${thresh}

cd stats

# non-FA processing
for diff in MD AD RD; do

	echo "working on ${diff}"

	fslmaths all_${diff} -mas mean_FA_mask all_${diff}

	echo "projecting all_${diff} onto mean FA skeleton"

	tbss_skeleton -i mean_FA -p ${thresh} mean_FA_skeleton_mask_dst ${FSLDIR}/data/standard/LowerCingulum_1mm all_FA all_${diff}_skeletonised -a all_${diff}

done


######################################
# EXTRACT TRACT SPECIFIC STATS       #
######################################

# warp atlas tracts in MNI space to mean_FA

if [ ! -f JHUlabels_thr25_groupspace.nii.gz ] ; then
	echo "running FLIRT on JHU ICBM FA"
	flirt -in ${atlastractsdir}/JHU-ICBM-FA-1mm.nii.gz -ref mean_FA.nii.gz -out JHU2FA -omat JHU2FA.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear
	echo "running FNIRT on JHU ICBM FA"
	fnirt --iout=JHU2FU_head --in=JHU-ICBM-FA-1mm.nii.gz --aff=JHU2FA.mat --cout=JHU2FA_warp --iout=JHU2FA --jout=JHU2FA_jac --ref=mean_FA.nii.gz --refmask=mean_FA_mask.nii.gz --warpres=10,10,10

	applywarp -i JHU-ICBM-FA-1mm.nii.gz -r mean_FA -o JHU2FA -w JHU2FA_warp
	applywarp -i ${atlastractsdir}/JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz -w JHU2FA_warp -r mean_FA -o JHUlabels_thr25_groupspace --interp=nn
fi
if [ ! -f JHUlabels_CC1mm_groupspace.nii.gz ] ; then
	applywarp -i ${atlastractsdir}/JHU-ICBM-labels-1mm.nii.gz -w JHU2FA_warp -r mean_FA -o JHUlabels_CC1mm_groupspace --interp=nn
fi


# extract specific tracts from JHU atlas in group space and multiply with skeleton
fslmaths JHUlabels_CC1mm_groupspace -thr 3 -uthr 3 -mul mean_FA_skeleton_mask.nii.gz -bin CC_genu.nii.gz
fslmaths JHUlabels_CC1mm_groupspace -thr 4 -uthr 4 -mul mean_FA_skeleton_mask.nii.gz -bin CC_body.nii.gz
fslmaths JHUlabels_CC1mm_groupspace -thr 5 -uthr 5 -mul mean_FA_skeleton_mask.nii.gz -bin CC_splenium.nii.gz
fslmaths JHUlabels_thr25_groupspace -thr 13 -uthr 13 -mul mean_FA_skeleton_mask.nii.gz -bin ILF_L.nii.gz
fslmaths JHUlabels_thr25_groupspace -thr 14 -uthr 14 -mul mean_FA_skeleton_mask.nii.gz -bin ILF_R.nii.gz
fslmaths JHUlabels_thr25_groupspace -thr 1 -uthr 1 -mul mean_FA_skeleton_mask.nii.gz -bin ATR_L.nii.gz
fslmaths JHUlabels_thr25_groupspace -thr 2 -uthr 2 -mul mean_FA_skeleton_mask.nii.gz -bin ATR_R.nii.gz

# make bilateral
fslmaths ATR_L.nii.gz -add ATR_R.nii.gz ATR_bil.nii.gz
fslmaths ILF_L.nii.gz -add ILF_R.nii.gz ILF_bil.nii.gz


for diff in FA MD AD RD; do

	echo ${diff}

	for ROI in CC_body CC_genu CC_splenium ATR_bil ILF_bil; do

		echo ${ROI}

		# extract median diffusivity within a ROI
		fslstats -t all_${diff}_skeletonised.nii.gz -k ${ROI} -p 50 > temp.txt
		paste subjects.txt temp.txt > ${diff}_${ROI}.txt
		rm temp.txt

	done
done

echo "DONE - finally!"
