#!/usr/bin/tcsh -f
#!/bin/tcsh -f

## provide name of directory where analyses are run as the dirname
set slicedir = /home/lesh/dtitools/Slicer-4.8.1-linux-amd64;
set dtiprepdir = /home/lesh/dtitools/DTIPrepTools-0.1.1-Linux/bin;
set rawdir = /nfs/inflammation/raw_data/scans/00_MONTH;
set dtioutdir = /nfs/inflammation/dti/preprocessed/00_MONTH/dtiprep_15_nii;
set niftidir = /nfs/inflammation/dti/preprocessed/00_MONTH/Nifti;
set startdir = `pwd`

##indicate what protocol xml file you would like to use here
#if starting from scratch please create new xml using dtiprep instructions
set protocol = /home/lesh/Documents/DTIPrep_minimal_nhp_15p_nifti.xml

echo Processing DTI QA using DTIPrep protocol ${protocol}
cd ${rawdir}

foreach sub (MMU46046_45M MMU46052_45M MMU46054_45M MMU46060_45M MMU46064_45M MMU46068_45M MMU46071_45M MMU46072_45M MMU46073_45M)
if ( -d ${niftidir}/${sub} ) then
	 #if raw data exists then proceed
	 echo Nifti data for ${sub} exist
	 if ( -d ${dtioutdir}/${sub} ) then
		  echo ${sub} directory already created
	 else
		  mkdir ${dtioutdir}/${sub}
	 endif
	 if ( -f ${dtioutdir}/${sub}/${sub}_DWI.nrrd ) then
		  echo ${sub} nrrd already created
	 else
		  echo Converting ${sub} Nifti to NRRD
		  ${slicedir}/Slicer --launch DWIConvert --inputVolume ${niftidir}/${sub}/multd/orig.nii.gz -o ${dtioutdir}/${sub}/${sub}_DWI.nrrd --inputBVectors ${niftidir}/${sub}/multd/orig_bvecs --inputBValues ${niftidir}/${sub}/multd/orig_bvals --conversionMode FSLToNrrd
	 endif
	 if ( -f ${dtioutdir}/${sub}/${sub}_DWI_XMLQCResult.xml ) then
		  echo ${sub} QC already completed
	 else
		  echo Performing DTIPrep QA for ${sub}
		  ${dtiprepdir}/DTIPrep -w ${dtioutdir}/${sub}/${sub}_DWI.nrrd -c -p ${protocol} --numberOfThreads 4
	 endif
	 if ( -f ${dtioutdir}/${sub}/${sub}_DWI_artifacts.txt ) then
		  echo Artifact file already created
	 else
		  echo Creating artifact file
		  cat  ${dtioutdir}/${sub}/${sub}_DWI_QCReport.txt | grep Slicewise_artifacts > ${dtioutdir}/${sub}/${sub}_DWI_artifacts.txt
		  cat  ${dtioutdir}/${sub}/${sub}_DWI_QCReport.txt | grep Interlacewise_artifacts >> ${dtioutdir}/${sub}/${sub}_DWI_artifacts.txt
	 endif
else
	 echo Subject ${sub} inflammation nifti data directory does not exist
endif
end

cd ${startdir}
