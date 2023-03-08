#!/usr/bin/tcsh -f
## on costinquad, bob
#!/bin/tcsh -f
## on amict
# FreeSurfer & FSL data importer
# ------------------------------------------------------------------

## ------------------------------------------------------------------
## setup block
## ------------------------------------------------------------------
## the project name
setenv prjname  '0110fw2';
## the machine 
setenv TCMACHINE amict;
setenv TCMACHINE costinquad;
setenv TCMACHINE `hostname`;
setenv TCMACHINE bob51;


setenv cr '/nfs/ep2/inflammation/code';
set src = ${cr};
setenv SRoot '/nfs/ep2/inflammation';
set dcm = ${SRoot}/raw_data/scans/00_MONTH;
set rez = ${SRoot}/dti/preprocessed/00_MONTH/Results;
set nii = ${SRoot}/dti/preprocessed/00_MONTH/Nifti;
set sub = ${SRoot}/mprage/preprocessed/00_MONTH

setenv ulb  ${cr}/proc2

prolog

foreach i ( $dcm  $nii $rez $sub )
	 #echo "$i";
end

## ------------------------------------------------------------------
## the sorted trees
## ------------------------------------------------------------------
if ( $#argv != 1 )  then 
	 echo "Usage: $0 ep_subject_name";
	 exit;
endif

### subject specific folders 
set snam = ${1}; 
echo "\nsnam = $snam"
set sdcm = ${dcm}/${snam}; 
set snii = ${nii}/${snam}; #mkdir -p ${snii}; 
set srez = ${rez}/${snam}; #mkdir -p ${srez};
set ssub = ${sub}/${snam}; #mkdir -p ${ssub};
 
foreach i ( sdcm snii srez )
	 eval set j = \$$i;
	 if ( ! -d $j ) then
		  echo "not a folder $j"; exit
	 endif
end


## ------------------------------------------------------------------
## DT native space
## ------------------------------------------------------------------

## set the fw folder
set fwf = $srez/dti.REC.bH_900_1400_bl_500_woBadVols;
## set the dt_fsl_proc folder
set dtfsl = $srez/dtt2ants; mkdir -p $dtfsl;
pushd $dtfsl > /dev/null;

echo "\nCopy the Dipy computed images to work on";
[ -f F.nii.gz ]|| fslmaths ${rez}/Dipy/${snam}/F.nii.gz F;
[ -f tF.nii.gz ]|| fslmaths ${rez}/Dipy/${snam}/tF.nii.gz tF;
[ -f balldata0_head.nii.gz ] || fslroi ${fwf}/dipy_prep/balldata_eddy.nii.gz balldata0_head.nii.gz 0 1;
[ -f balldata0_brain.nii.gz ] || fslmaths balldata0_head.nii.gz -mas ${fwf}/topup/both_b0/unwarp/Tiout_b0_${snam}_brain_mask.nii.gz balldata0_brain.nii.gz;
[ -f FA.nii.gz ]|| fslmaths ${rez}/Dipy/${snam}/FA FA;
[ -f tFA.nii.gz ]|| fslmaths ${rez}/Dipy/${snam}/tFA tFA;

set strucroot = '/nfs/inflammation/mprage/unc_ucd_nhp_final_structs_Dec2019'
set nam = `echo ${snam} | awk -F_ '{print $1}'`
if ( ! -f raw_reorient.mat ) then 
 echo "-1 0 0 140\n0 -1 0 140\n0 0 1 0\n0 0 0 1\n" > raw_reorient.mat
endif
#copy over structurals that have been segmented and parcellated
echo Copying ${snam} structurals
[ -f t2w.nii.gz ] ||\
	 cp ${strucroot}/${nam}/${snam}*/sMRI/atlasReg/${snam}_*_T2w_atlasT1Reg_norm.nii.gz t2w.nii.gz;
[ -f t1w.nii.gz ] ||\
	 cp ${strucroot}/${nam}/${snam}*/sMRI/atlasReg/${snam}_*_T1w_atlasReg_norm.nii.gz t1w.nii.gz;
[ -f t2w_brain.nii.gz ] ||\
	 cp ${strucroot}/${nam}/${snam}*/sMRI/Bias/${snam}_*_T2w_atlasT1Reg_N4GM_final_strip.nii.gz t2w_brain.nii.gz;
[ -f t2w_brainmask.nii.gz ] ||\
	 cp ${strucroot}/${nam}/${snam}*/sMRI/Strip/${snam}_*_T2w_atlasT1Reg_MaskManual.nii.gz t2w_brainmask.nii.gz;

set moveimage = 'balldata0_head'
set refimage = 't2w'

# align b0 to t2
echo Rigid weighted alignment of b0 to t2
[ -f b02t2_6dof.nii.gz ] ||\
	 flirt -in  ${moveimage} -ref ${refimage} -out b02t2_6dof -dof 6 -omat b02t2_6dof.mat -searchrz -180 180 -searchry -180 180 -searchrx -180 180
# carry over flirt registered dipy results here, this is mostly done for QC purposes
echo Applying affine alignment to free water images
[ -f F2t2_6dof.nii.gz ] ||\
     flirt -in F -ref ${refimage} -out F2t2_6dof -applyxfm -init b02t2_6dof.mat;
[ -f FA2t2_6dof.nii.gz ] ||\
     flirt -in FA -ref ${refimage} -out FA2t2_6dof -applyxfm -init b02t2_6dof.mat;

#compute inverse transformation to get t2 to b0 if needed
[ -f t22b0_6dof.mat ] || convert_xfm -omat t22b0_6dof.mat -inverse b02t2_6dof.mat;

#convert fslmat to ANTS format
[ -f b02t2_6dof.tfm ] ||\
	 /nfs/inflammation/code/primate-MRI-processing-master/convert3D/c3d-1.0.0-Linux-x86_64/bin/c3d_affine_tool -ref ${refimage}.nii.gz -src \
	 ${moveimage}.nii.gz b02t2_6dof.mat -fsl2ras -oitk b02t2_6dof.tfm

## nonlinear DTI 2 orig brain (b0)
## higher SyN values, more iterations leads to greater deformations (default is SyN of .25)
## Gauss is the smoothing sigma of the deformation field (default is 3.0)
## WARNING these settings often need some adjustment depending upon data type do not assume they are the best for you
echo Compute nonlinear transformation of b0 to t2
if ( ! -f ${dtfsl}/b02t2_nlAffine.txt ) then
	 ANTS 3 -m CC\[${refimage}.nii.gz,${moveimage}.nii.gz,1,5\] -o ${dtfsl}/b02t2_nl.nii.gz -i 120x80x30 -r Gauss\[1,0\] -t SyN\[0.1\] -a b02t2_6dof.tfm -x t2w_brainmask.nii.gz --continue-affine false --use-Histogram-Matching;
endif

#threshold free water data to eliminate voxels failing nonlinear model fit (when the WLS model fails to fit the value is much higher, closer to .001 and the threshold would need to be adjusted)

[ -f tFraw.1.7.nii.gz ]||\
	 fslmaths tF.nii.gz -thr .00001 tFraw.1.7.nii.gz;
	 
echo Apply nonlinear transforms to diffusion images
#this step will apply the affine and nonlinear warps to your native space diffusion data and bring them to the T2 space of the subject
[ -f b02t2_nl.nii.gz ]||\
	 WarpImageMultiTransform 3 \
	 ${moveimage}.nii.gz b02t2_nl.nii.gz b02t2_nlWarp.nii.gz -R ${refimage}.nii.gz b02t2_nlAffine.txt;
[ -f tFA2t2_nl.nii.gz ]||\
	 WarpImageMultiTransform 3 \
	 tFA.nii.gz tFA2t2_nl.nii.gz b02t2_nlWarp.nii.gz -R ${refimage}.nii.gz b02t2_nlAffine.txt;
[ -f tFraw.1.7tot2_nl.nii.gz ]||\
	 WarpImageMultiTransform 3 \
	 tFraw.1.7.nii.gz tFraw.1.7tot2_nl.nii.gz b02t2_nlWarp.nii.gz -R ${refimage}.nii.gz b02t2_nlAffine.txt;

popd > /dev/null; ##pushd $dtfsl > /dev/null;


