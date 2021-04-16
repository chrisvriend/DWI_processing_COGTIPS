#!/bin/bash

# (C) C. vriend - Amsterdam UMC dpt. Psychiatry / Anatomy and Neurosciences - 2020
# for details on the individual steps see the MRtrix3 documentation: mrtrix.readthedocs.io/en/latest/

#necessary software: mrtrix3  + FSL 6.0.1 + Freesurfer 6 or 7, ANTS

# assumes that DWI data has already been motion and distortion corrected (using eddy) and that freesurfer has been run on the T1 weighted scan.
# in current script there are 3 shells: b1000, b2000, b3000. Change according to your own data

#--------------------------------------------------------------------------------------------------------
# ASSUMED DATA STRUCTURE
#--------------------------------------------------------------------------------------------------------
# headdirectory
#	|
#       |________DWI directory
#	|
#	|		       |____ data, bvals, bvecs
#	|
#	|________T1 directory
#		|
#		|-------
#			    |		Single subject T1 image + skull stripped or robust template of T0 - T1 scans in case of longitudinal data
#			|____
#			    |		Freesurfer
#   		        	|____ FreeSurfer folder w/ DK, BNA, ... atlas
#

#
#---------------------------------------------------------------------------------------------------------


# input variables
headdir=/path/to/input

# folder with atlas specific relabeling text files
labeldir=/path/to/labeltxtfiles
Ncores=30 # number of cores to use


subj=sub-0001
# path to subject specific T1w image
T1_in=${headdir}/T1/${subj}_T1w.nii.gz
# path to subject specific skull stripped T1w image
T1_BET_in=${headdir}/T1/${subj}_T1w_brain.nii.gz
# path to subject specific freesurfer output directory
FSdir=${headdir}/T1/${subj}_FS



cd ${headdir}/DWI
###############################################################################
# CREATE nodif and MASK                                                       #
###############################################################################
if [ ! -f nodif.nii.gz ]; then

# data.nii.gz = preprocessed multishell DWI data file

fslroi ${subj}_data dwinodif 0 2
fslmaths dwinodif -Tmean nodif
bet2 nodif ${subj}_nodif_brain -m -f 0.2
fi

cp nodif_brain_mask.nii.gz ${subj}_data_mask.nii.gz

# this command doesn't work well; holes in the mask that interferes with seeding in the white matter during tckgen
#dwi2mask data.nii.gz data_mask.nii.gz -fslgrad bvecs bvals -force

DWI_nodif=${headdir}/DWI/${subj}_nodif_brain.nii.gz

###############################################################################
# T1 to DWI registration                                                      #
###############################################################################

mkdir -p xfms

flirt -in ${T1_BET_in} -ref ${DWI_nodif} -omat ${headdir}/DWI/xfms/${subj}_str2diff.mat -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -cost corratio
convert_xfm -omat ${headdir}/DWI/xfms/${subj}_diff2str.mat -inverse ${headdir}/DWI/xfms/${subj}_str2diff.mat
epi_reg --epi=${DWI_nodif} --t1=${T1_in} --t1brain=${T1_BET_in} --out=${headdir}/DWI/xfms/${subj}_epireg
convert_xfm -omat ${headdir}/DWI/xfms/${subj}_epireg_inversed.mat -inverse ${headdir}/DWI/xfms/${subj}_epireg.mat

###############################################################################
# Derive tissue-segmented image (generate 5TT data) (MRtrix)                  #
###############################################################################
     # premasked with ANTS (NO SKULL)
     # with nthreads >0 hpc invokes unavailable veryshort.q
     # nocrop otherwise warping to DWI may fail
cd ${headdir}/T1

if [ ! -f ${subj}_5TT.nii.gz ]; then

5ttgen fsl ${T1_BET_in} ${subj}_5TT.nii.gz -premasked -nthreads 16 -nocrop -scratch ${headdir}/temp_5ttgen

rm -rf ${headdir}/temp_5ttgen
else
echo "${subj}_5TT.nii.gz already exists"
fi

# apply DWIwarping to 5TT scan
flirt -in ${subj}_5TT.nii.gz -ref ${DWI_nodif} -applyxfm -init ${headdir}/DWI/xfms/${subj}_epireg_inversed.mat -out ${headdir}/DWI/${subj}_5TT2DWI -interp nearestneighbour

mrconvert ${headdir}/DWI/${subj}_5TT2DWI.nii.gz ${headdir}/DWI/${subj}_5TT2DWI.mif


# to visualize
#-- 5tt2vis t1_mpr_ns_sag_iso_5TT.mif t1_mpr_ns_sag_iso_5TT_vis.mif

###############################################################################
# FREESURFER to DWI registration                                              #
###############################################################################
# mri_vol2vol | Invert the transform. The movvol becomes the geometry template for the output, and the targvol becomes the input that will be resampled.

SUBJECTS_DIR=${headdir}/T1


bbregister --s ${subj} --mov ${DWI_nodif} --init-fsl --reg ${FSdir}/DWI/register.dat --dti

# this assumes that you have already previously warped the brainnetome atlas to freesurfer space
# https://atlas.brainnetome.org/brainnetome.html

# Brainnetome atlas + subcortical structures
mri_vol2vol --mov ${DWI_nodif} --targ ${FSdir}/mri/BNA+aseg.mgz \
--o ${FSdir}/mri/BNA+aseg_wprd.mgz \
--reg ${FSdir}/DWI/register.dat --inv --no-save-reg --interp nearest

# DK - native freesurfer parcellation
mri_vol2vol --mov ${DWI_nodif} --targ ${FSdir}/mri/aparc+aseg.mgz \
--o ${FSdir}/mri/aparc+aseg_wprd.mgz \
--reg ${FSdir}/DWI/register.dat --inv --no-save-reg --interp nearest



###############################################################################
# REGISTER LABELS TO DWI                                                      #
###############################################################################
# https://mrtrix.readthedocs.io/en/latest/quantitative_structural_connectivity/labelconvert_tutorial.html

mrconvert ${FSdir}/mri/aparc+aseg_wprd.mgz ${headdir}/DWI/${subj}_FS_DK_orig.nii.gz
mrconvert ${FSdir}/mri/BNA+aseg_wprd.mgz ${headdir}/DWI/${subj}_FS_BNA_orig.nii.gz

# DK atlas (orig - modified labels)
labelconvert ${headdir}/DWI/${subj}_FS_DK_orig.nii.gz ${labeldir}/DK_labels_orig.txt ${labeldir}/MRTRIX_relabels/DK_labels_modified.txt ${headdir}/DWI/${subj}_FS_DK_relab.nii.gz

# BN Atlas
labelconvert ${headdir}/DWI/${subj}_FS_BNA_orig.nii.gz ${labeldir}/BNA_labels_orig.txt ${labeldir}/MRTRIX_relabels/BNA_labels_modified.txt ${headdir}/DWI/${subj}_FS_BNA_relab.nii.gz


###############################################################################
# DWI PROCESSING                                                   #
###############################################################################

cd ${headdir}/DWI


# BIAS CORRECTION
  # - with fsl or ANTS
if [ ! -f ${subj}_data_biascor.nii.gz ]; then

dwibiascorrect ants ${subj}_data.nii.gz ${subj}_data_biascor.nii.gz -nthreads ${Ncores} -bias bias_est.nii.gz -scratch ${headdir}/tempbiascorrect -force -fslgrad ${subj}_bvecs ${subj}_bvals

fi

# N4BiasFieldCorrection option -b. [initial mesh resolution in mm, spline order]

rm -rf ${headdir}/tempbiascorrect



###############################################################################
# ESTIMATE RESPONSE FUNCTION                                                  #
###############################################################################
# multi shell = dhollander or msmt_5tt

# msm_5tt - change shells if necessary
dwi2response msmt_5tt ${subj}_data_biascor.nii.gz ${subj}_5TT2DWI.mif ${subj}_wm_response.txt ${subj}_gm_response.txt ${subj}_csf_response.txt -shell 0,1000,2000,3000 -nthreads ${Ncores} -fslgrad ${subj}_bvecs ${subj}_bvals -force -tempdir ${headdir}/tempdwi2response

# clean up
rm -rf ${headdir}/tempdwi2response

# for visual inspection using mrview
#shview <Output response text file>

###############################################################################
# DWI SPHERICAL DECONVOLUTION                                                 #
###############################################################################

dwi2fod msmt_csd ${subj}_data_biascor.nii.gz ${subj}_wm_response.txt ${subj}_wmfod.mif ${subj}_gm_response.txt ${subj}_gm.mif ${subj}_csf_response.txt ${subj}_csf.mif \
-mask ${subj}_data_mask.nii.gz -fslgrad ${subj}_bvecs ${subj}_bvals -shell 0,1000,2000,3000 -nthreads ${Ncores} -force

# for visualisatioN | generates a 4D image with 3 volumes, corresponding to the tissue densities of CSF, GM and WM, which will then be displayed in mrview as an RGB image with CSF as red, GM as green and WM as blue (as was presented in the MSMT CSD manuscript).
mrconvert ${subj}_wmfod.mif - -coord 3 0 | mrcat ${subj}_csf.mif ${subj}_gm.mif - ${subj}_tissueRGB.mif -axis 3
mrconvert ${subj}_tissueRGB.mif ${subj}_tissueRGB.nii.gz

#mrview <Input DWI> -odf.load_sh <Output FOD image>

###############################################################################
# DWI STREAMLINE TRACTOGRAPHY                                                 #
###############################################################################

# since the propagation and termination of streamlines is primarily handled by the 5TT image, it is no longer necessary to provide a mask using the -mask option.
#In fact, for whole-brain tractography, it is recommend that you _not_ provide such an image when using ACT: depending on the accuracy of the DWI brain mask,
# its inclusion may only cause erroneous termination of streamlines inside the white matter due to exiting this mask. ( -mask ###)

# this step takes a very long time to run
tckgen wmfod.mif ${headdir}/${subj}_dwi_100M.tck -seed_image ${subj}_data_mask.nii.gz -number 100M -force -maxlength 250 -act ${headdir}/DWI/${subj}_5TT2DWI.mif -nthreads ${Ncores}

tckmap ${headdir}/${subj}_dwi_100M.tck ${headdir}/DWI/DWI_${subj}_100M_dens.nii.gz -template ${headdir}/DWI/${subj}_data_mask.nii.gz -force -nthreads ${Ncores}
# visualization
# mrview <Input DWI> -tractography.load <Output track file>


###############################################################################
# DWI Spherical-deconvolution Informed Filtering of Tractograms (SIFT)        #
###############################################################################
# https://mrtrix.readthedocs.io/en/latest/quantitative_structural_connectivity/sift.html?highlight=tcksift
#the number of streamlines connecting two regions of the brain becomes a proportional estimate of the total cross-sectional area of the white matter fibre pathway connecting those regions; this is inherently a highly biologically relevant measure of ‘structural connectivity’

# EM framework to find an appropriate cross-section multiplier for each streamline

tcksift2 ${headdir}/${subj}_dwi_100M.tck ${headdir}/DWI/${subj}_wmfod.mif ${headdir}/${subj}_sift.txt -act ${headdir}/DWI/${subj}_5TT2DWI.mif -force -nthreads 14

# -act image use an ACT five-tissue-type segmented anatomical image to derive the processing mask
tckmap ${headdir}/${subj}_dwi_100M.tck ${headdir}/DWI/DWI_${subj}_100M_dens_sift.nii.gz -template ${headdir}/DWI/${subj}_data_mask.nii.gz -tck_weights ${headdir}/${subj}_sift.txt -force -nthreads ${Ncores}


###############################################################################
# DWI CONNECTOME GENERATION                                                   #
###############################################################################
# https://mrtrix.readthedocs.io/en/latest/quantitative_structural_connectivity/structural_connectome.html


# --- Brainnetome Atlas

# streamline count - optional scale each contribution to the connectome edge by the inverse of the two node volumes with -scale_invnodevol but this is discouraged by the developers of MRtrix3
# e.g . https://community.mrtrix.org/t/tck2connectome/1133/2
tck2connectome ${headdir}/${subj}_dwi_100M.tck ${headdir}/DWI/${subj}_FS_BNA_relab.nii.gz ${headdir}/DWI/${subj}_connectome.BNA.csv -zero_diagonal -tck_weights ${headdir}/${subj}_sift.txt -nthreads ${Ncores}
# --- DK atlas
tck2connectome ${headdir}/DWI/${subj}_dwi_100M.tck ${headdir}/DWI/${subj}_FS_DK_relab.nii.gz ${headdir}/DWI/${subj}_connectome.DK.csv -zero_diagonal -tck_weights_in ${headdir}/DWI/${subj}_sift.txt -nthreads ${Ncores} -symmetric
