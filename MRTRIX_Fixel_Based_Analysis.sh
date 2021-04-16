#!/bin/bash

# (C) C. vriend - Amsterdam UMC dpt. Psychiatry / Anatomy and Neurosciences - 2020
# documentation:
# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
# tract specific:
#https://community.mrtrix.org/t/how-to-restrict-fixelcfestats-to-subset-of-tracts/845/2

# --------------------------------------------------------------------------------------------
#  SOURCE SOFTWARE
#fsl 6.0.1
#mrtrix3

# -----------------------------------------------------------------------------------------------
# SET UP INPUT VARIABLES
# N simultaneous jobs
cores=30
# N threads to use per job
threads=10

# -----------------------------------------------------------------------------------------------
# SET UP INPUT DIRECTORIES
headdir=/path/to/inputdir
templatedir=${headdir}/group_template
atlastractsdir=${headdir}/atlastracts

# -----------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------
# ASSUMED DATA STRUCTURE
# one folder for each subject containing the eddy corrected DWI data (data.nii.gz) and the accompanying bvals and bvecs files
# script assumes that input has already been denoised and eddy corrected
# in the script below each subject specific folder starts with 735 (e.g. 735_001, 735_002, etc.) Change according to your own sample coding.
# for each allows the commands behind the :  to be run in parallel for -nthreads = XX
# specifying -nthreads behind the :  will use that number of cores per subject.
# BE CAREFUL in setting -nthreads before and after the : to as not to overload your system


echo "convert to mif"
for_each -nthreads ${threads} 735* : mrconvert IN/data.nii.gz IN/data.mif -fslgrad IN/bvecs IN/bvals
#for_each -nthreads 10 735* : mrconvert IN/data_biascor.nii.gz IN/data_biascor.mif -fslgrad IN/bvecs IN/bvals


# bias field correction (if not yet applied)
for_each -nthreads ${threads} 735* : dwibiascorrect ants IN/data.mif IN/data_biascor.mif

# compute tissue response functions
for_each -nthreads ${threads} 735* : dwi2response dhollander IN/data_biascor.mif IN/wm_response.txt IN/gm_response.txt IN/csf_response.txt


mkdir -p ${templatedir}/fod_input
mkdir -p ${templatedir}/mask_input

# group average response functions
if [ ! -f ${templatedir}/group_average_response_gm.txt ]; then
  responsemean 735*/gm_response.txt ${templatedir}/group_average_response_gm.txt
fi
if [ ! -f ${templatedir}/group_average_response_wm.txt ]; then
  responsemean 735*/wm_response.txt ${templatedir}/group_average_response_wm.txt
fi
if [ ! -f ${templatedir}/group_average_response_csf.txt ]; then
  responsemean 735*/csf_response.txt ${templatedir}/group_average_response_csf.txt
fi


# convert to mif format
# upsample to 1.25 mm resolution

echo "upsample to 1.25 mm resolution"
for_each -nthreads ${threads} 735* : mrgrid IN/data_biascor.mif regrid -vox 1.25 IN/data_biascor_upsample.mif -force
echo "convert upsampled DWI from mif to nii.gz"

# run loop for all subject specific folders that start with 735* (=sample code)
# for_each doesn't always properly for this step
for subj in $(ls -d 735*); do
  mrconvert ${subj}/data_biascor_upsample.mif ${subj}/data_biascor_upsample.nii.gz -nthreads 6
done

# for_each -nthreads ${threads} 735* : mrconvert IN/data_biascor_upsample.mif IN/data_biascor_upsample.nii.gz

# create mask using BET because dwimask does not always do a very good job
for subj in $(ls -d 735*); do
  echo ${subj}
  cd ${subj}
  #if [ ! -f data_mask.nii.gz ]; then
  bet2 data_biascor_upsample.nii.gz data -m -f 0.3
  #fi
  cd ..
done


echo
for_each -nthreads ${threads} 735* : mrconvert IN/data_mask.nii.gz IN/data_mask_upsample.mif -force


# fibre orientation distribution estimation
#if [ ! -f wmfod.mif ] | [ ! -f gm.mif ] \ [ ! -f csf.mif ]
echo
echo "estimate fibre orientation distribution (FOD) using average response functions"
echo

# change shell number according to your own dataset
for_each -nthreads ${threads} 735* : dwi2fod msmt_csd IN/data_biascor_upsample.mif \
${templatedir}/group_average_response_wm.txt IN/wmfod.mif \
${templatedir}/group_average_response_gm.txt IN/gm.mif \
${templatedir}/group_average_response_csf.txt IN/csf.mif \
-mask IN/data_mask_upsample.mif \
-shell 0,1000,2000,3000 \
-force

echo
echo "joint bias field correction and intensity normalisation"
echo

for_each -nthreads ${threads} 735* : mtnormalise IN/wmfod.mif IN/wmfod_norm.mif IN/gm.mif IN/gm_norm.mif IN/csf.mif IN/csf_norm.mif -mask IN/data_mask_upsample.mif

cd ${headdir}
fi

echo "create symbolic links"
for_each 735* : ln -sr IN/wmfod_norm.mif ${templatedir}/fod_input/PRE.mif
for_each 735* : ln -sr IN/data_mask_upsample.mif ${templatedir}/mask_input/PRE.mif


echo
echo "Generate a study-specific unbiased FOD template"
echo

if [ ! -f ${templatedir}/wmfod_template.mif ]; then
  population_template ${templatedir}/fod_input -mask_dir ${templatedir}/mask_input ${templatedir}/wmfod_template.mif -voxel_size 1.25 -nthreads ${cores}
fi

echo "Register subject FOD images to the FOD template"
echo


for_each -nthreads ${threads} 735* : mrregister IN/wmfod_norm.mif -mask1 IN/data_mask_upsample.mif ${templatedir}/wmfod_template.mif -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif
#for_each -nthreads 6 735* : mrregister IN/wmfod_norm.mif -mask1 IN/data_mask_upsample.mif ${templatedir}/wmfod_template.mif -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif -nthreads 3
for_each -nthreads ${threads} 735* : mrtransform IN/data_mask_upsample.mif -warp IN/subject2template_warp.mif -interp nearest -datatype bit IN/data_mask_in_template_space.mif

echo
echo "Compute the template mask (intersection of all subject masks in template space)"
mrmath 735*/data_mask_in_template_space.mif min ${templatedir}/template_mask.mif -datatype bit
# It is absolutely crucial to check at this stage that the resulting template mask includes all regions of the brain that are intended to be analysed

echo
echo "Compute a white matter template analysis fixel mask"
fod2fixel -mask ${templatedir}/template_mask.mif -fmls_peak_value 0.06 ${templatedir}/wmfod_template.mif ${templatedir}/fixel_mask -nthreads ${cores}


echo "Warp FOD images to template space"
echo

#if [ ! -f fod_in_template_space_NOT_REORIENTED.mif ]; then
for_each -nthreads ${threads} 735* : mrtransform IN/wmfod_norm.mif -warp IN/subject2template_warp.mif -reorient_fod no IN/fod_in_template_space_NOT_REORIENTED.mif
#fi

echo
echo "Segment FOD images to estimate fixels and their apparent fibre density (FD)"

for_each -nthreads ${threads} 735* : fod2fixel -mask ${templatedir}/template_mask.mif IN/fod_in_template_space_NOT_REORIENTED.mif IN/fixel_in_template_space_NOT_REORIENTED -afd fd.mif

echo
echo "Reorient fixels"
for_each -nthreads ${threads} 735* : fixelreorient IN/fixel_in_template_space_NOT_REORIENTED IN/subject2template_warp.mif IN/fixel_in_template_space

# clean up
rm -r 735*/fixel_in_template_space_NOT_REORIENTED


echo
echo "Assign subject fixels to template fixels"
# matching the fixels of each individual subject to the single common set of template fixels (which then inherently also defines how they match across subjects). This is achieved by, for each fixel in the template fixel mask, identifying the corresponding fixel in the matching voxel of the subject image and assigning the FD value of this corresponding subject fixel to that fixel in template space. If no fixel exists or can be found in a subject that corresponds to a given template fixel then it is assigned a value of zero (as the absence of a subject fixel at this stage is most likely due to a very low, or even zero, FD).

for_each 735* : fixelcorrespondence IN/fixel_in_template_space/fd.mif ${templatedir}/fixel_mask ${templatedir}/fd PRE.mif -nthreads ${cores}

#Note that the output fixel directory ${templatedir}/fd is the same for all subjects. This makes sense, since after this operation, there is only a single remaining set of fixels (i.e. the template fixels), with corresponding FD values as obtained from each subject. This resulting directory ${templatedir}/fd now stores these data as individual fixel data files: one for each subject, and all with respect to a single set of corresponding template fixels. This way of storing the entire population’s FD data is then ready for input to fixelcfestats later on.


# Compute the fibre cross-section (FC) metric

for_each 735* : warp2metric IN/subject2template_warp.mif -fc ${templatedir}/fixel_mask ${templatedir}/fc IN.mif -nthreads ${cores}


mkdir ${templatedir}/logfc
cp ${templatedir}/fc/index.mif ${templatedir}/fc/directions.mif ${templatedir}/logfc
for_each 735* : mrcalc ${templatedir}/fc/IN.mif -log ${templatedir}/logfc/IN.mif


#Compute a combined measure of fibre density and cross-section (FDC)

mkdir ${templatedir}/fdc
cp ${templatedir}/fc/index.mif ${templatedir}/fdc
cp ${templatedir}/fc/directions.mif ${templatedir}/fdc
for_each 735* : mrcalc ${templatedir}/fd/IN.mif ${templatedir}/fc/IN.mif -mult ${templatedir}/fdc/IN.mif

cd ${templatedir}
# Perform whole-brain fibre tractography on the FOD template
# https://github.com/MRtrix3/mrtrix3/issues/1704
tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_dynamic wmfod_template.mif -mask template_mask.mif -select 10M -cutoff 0.06 tracks_dyn_10M.tck -nthreads ${cores}

# Reduce biases in tractogram densities - use tcksift not tcksift2 here
tcksift tracks_dyn_10M.tck wmfod_template.mif tracks_dyn_1M_sift.tck -term_number 1000000

#The appropriate FOD amplitude cutoff for FOD template tractography can vary considerably between different datasets, as well as different versions of MRtrix3 due to historical software bugs. While the value of 0.06 is suggested as a reasonable value for multi-tissue data, it may be beneficial to first generate a smaller number of streamlines (e.g. 100,000) using this value, and visually confirm that the generated streamlines exhibit an appropriate extent of propagation at the ends of white matter pathways, before committing to generation of the dense tractogram.

# Generate fixel-fixel connectivity matrix
fixelconnectivity fixel_mask/ tracks_dyn_1M_sift.tck matrix/ -nthreads ${cores}
#fixelconnectivity fixel_mask/ tracks_dyn_10M.tck -tck_weights_in sift.txt matrix/


# Smooth fixel data using fixel-fixel connectivity
fixelfilter fd smooth fd_smooth -matrix matrix/ -force -nthreads ${cores}
fixelfilter logfc smooth logfc_smooth -matrix matrix/ -force -nthreads ${cores}
fixelfilter fdc smooth fdc_smooth -matrix matrix/ -force -nthreads ${cores}

# Perform statistical analysis of FD, FC, and FDC
fixelcfestats fd_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ stats_fd/
fixelcfestats logfc_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ statslog_fc/
fixelcfestats fdc_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ stats_fdc/

#extracting roi
# https://community.mrtrix.org/t/extracting-rois-for-fixel-based-analysis/1614/2
# restrict stats to subset of tracts
#https://community.mrtrix.org/t/how-to-restrict-fixelcfestats-to-subset-of-tracts/845
# Visualise the results
# To view the results load the population FOD template image in mrview, and overlay the fixel images using the vector plot tool. Note that p-value images are saved as (1 - p-value). Therefore to visualise all results at a threshold of p < 0.05, within the mrview fixel plot tool, apply a lower threshold at a value of 0.95.


######################################
# EXTRACT TRACT SPECIFIC STATS       #
######################################

cd ${headdir}
mkdir xfms
mrconvert template_mask.mif template_mask.nii.gz

# warp atlas tracts in MNI space to tracks in group template
flirt -in ${atlastractsdir}/JHU-ICBM-FA-1mm.nii.gz -ref tracks_10M.nii.gz -out xfms/JHU2tracks \
-omat xfms/JHU2tracks.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear
fnirt --iout=xfms/JHU2tracks_head --in=${atlastractdir}/JHU-ICBM-FA-1mm.nii.gz --aff=xfms/JHU2tracks.mat \
--cout=xfms/JHU2tracks_warp --iout=xfms/JHU2tracks --jout=xfms/JHU2tracks_jac --ref=tracks_10M.nii.gz --refmask=template_mask.nii.gz --warpres=10,10,10

applywarp -i ${atlastractsdir}/JHU-ICBM-FA-1mm.nii.gz -r tracks_10M.nii.gz -o JHU2tracks -w xfms/JHU2tracks_warp

#applywarp -i JHU-ICBM-tracts-maxprob-thr50-1mm.nii.gz -w JHU2FA_warp -r mean_FA -o JHUlabels_thr50_groupspace --interp=nn
applywarp -i ${atlastractsdir}/JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz -w ${atlastractsdir}/JHU2tracks_warp -r tracks_10M.nii.gz -o JHUlabels_thr25_groupspace --interp=nn


# example of tracts to use
fslmaths JHUlabels_thr25_groupspace -thr 1 -uthr 1 -bin ATR_L_vx.nii.gz
fslmaths JHUlabels_thr25_groupspace -thr 2 -uthr 2 -bin ATR_R_vx.nii.gz

fslmaths ATR_L_vx.nii.gz -add ATR_R_vx.nii.gz -bin ATR_bil_vx.nii.gz

mrconvert ATR_bil_vx.nii.gz ATR_bil_vx.mif

# convert voxel to fixel
voxel2fixel ATR_bil_vx.mif ${templatedir}/fixel_mask/ ATRinfixelspace/ ATR_bil_fix.mif


# https://community.mrtrix.org/t/calculating-average-fba-metrics-of-specific-tracts/1805/2


echo "extract fiber-based stats"
rm -f fdc_smooth_stats.txt fd_smooth_stats.txt logfc_smooth_stats.txt
for FBM in fdc_smooth fd_smooth logfc_smooth; do
  echo ${FBM}
  cd ${FBM}

  rm -f ${FBM}temp_median.txt ${FBM}temp_mean.txt subjects.txt
  echo "subj" > subjects.txt
  echo "mean" > ${FBM}temp_mean.txt
  echo "median" > ${FBM}temp_median.txt

  # change 735 according to your own data
  for scan in $(ls 735*.mif); do

    subj=${scan%.mif*}
    echo ${subj} >> subjects.txt

    mrstats ${scan} -mask ${headdir}/ATRinfixelspace/ATR_bil_fix.mif -ignorezero -output mean >> ${FBM}temp_mean.txt
    mrstats ${scan} -mask ${headdir}/ATRinfixelspace/ATR_bil_fix.mif -ignorezero -output median >> ${FBM}temp_median.txt

  done
  paste subjects.txt ${FBM}temp_mean.txt ${FBM}temp_median.txt > ../${FBM}_stats.txt
  cd ..

done
