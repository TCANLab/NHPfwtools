#!/usr/bin/tcsh

setenv cr /nfs/inflammation/code;
###################
set SRoot = $cr;

set Res = "/nfs/inflammation/dti/preprocessed/00_MONTH/Results"
pushd $cr > /dev/null
source set_FS5.3bob.sh

foreach age ( '6M' '12M' '24M' '36M' '45M' )

##uncomment sufx depending on whether wanting white matter or gray matter or could add another loop
#foreach sufx ( 'UNCParc_1' 'UNCParc_2' 'cleanWM')
foreach sufx (  'cleanWMGM_nosubcort' 'UNCParc_MultiAtlas_2' 'UNCParc_MultiAtlas_1')
# ------------------------------------------------------------------
#find all subjects
#choose a set all to run all subjects of specific ones
#set all = `( cd /nfs/inflammation/dti/preprocessed/00_MONTH/Results; ls -1d MMU?????_${age})`
set all = "MMU46046_${age} MMU46052_${age} MMU46054_${age} MMU46060_${age} MMU46064_${age} MMU46068_${age} MMU46071_${age} MMU46072_${age} MMU46073_${age}"

#process i subjects
foreach i ( $all )
	 echo "Running $i";
	 set j = `echo $i | sed 's/_'${age}'//'`;
	 set structdir = "/nfs/inflammation/mprage/unc_ucd_nhp_final_structs_Dec2019";
	 set segdir = `dirname ${structdir}/$j/${i}_*/sMRI/Seg/*_MultiAtlas_2.nii.gz`;
	 set prefixname = `basename ${structdir}/$j/${i}_*/sMRI/Seg/*_MultiAtlas_2.nii.gz | awk -F"_T1w" '{print $1}'`;
	 echo Removing subcortical regions from $i
	  [ -f ${segdir}/${i}_L_subcort_mask_inv.nii.gz ]|| fslmaths ${segdir}/${prefixname}_T1w_atlasReg_UNCParc_MultiAtlas_2.nii.gz -uthr 16.1 -thr 15.9 -binv ${segdir}/${i}_L_subcort_mask_inv.nii.gz;
	  [ -f ${segdir}/${i}_R_subcort_mask_inv.nii.gz ]|| fslmaths ${segdir}/${prefixname}_T1w_atlasReg_UNCParc_MultiAtlas_2.nii.gz -uthr 3.1 -thr 2.9 -binv ${segdir}/${i}_R_subcort_mask_inv.nii.gz;
	  [ -f ${segdir}/${i}_T1w_atlasReg_label_cleanWMGM_nosubcort.nii.gz ]|| fslmaths ${segdir}/${prefixname}_T1w_atlasReg_label_cleanWMGM.nii.gz -mas ${segdir}/${i}_R_subcort_mask_inv.nii.gz -mas ${segdir}/${i}_L_subcort_mask_inv.nii.gz ${segdir}/${i}_T1w_atlasReg_label_cleanWMGM_nosubcort.nii.gz;
	 set matlas = ${structdir}/$j/${i}_*/sMRI/Seg/*_${sufx}.nii.gz;
	 set onam = fmap_${i}_${sufx}.rez;
	 #echo $matlas;
	 set fmap = "${Res}/${i}/dtt2ants/tF2t2_nl.nii.gz";
	 #make directories to output masks and text output
	 mkdir -p /nfs/inflammation/mprage/masks/${i};
	 pushd /nfs/inflammation/mprage/masks/${i} > /dev/null;
	 # loop over maps
	 rm -f $onam;
	 touch $onam;
	 #eliminate values below .00001 where model fitting failed if using NLS
	 [ -f ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz ]|| fslmaths ${fmap} -thr .00001 ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz;
	 if ( ${sufx} == "cleanWMGM_nosubcort" ) then
		  set k = "";
		  echo Computing whole-brain masks for ${i}
		  foreach k ( 1 2 )
				#echo $k;
				set km = ` echo ${k} - .1 | bc -l `;
				# echo $km;
				set kM = ` echo ${k} + .1 | bc -l `;
				# echo $kM;
				# extract the msk_1, etc.
				echo Computing whole-brain mask ${k} for ${i}
				[ -f msk_${k}_${sufx}.nii.gz ]||\
					 fslmaths $matlas -thr $km -uthr $kM -bin msk_${k}_${sufx};
				echo "${j}, ${age}, ${k}, ${sufx}, `fslstats ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz -k msk_${k}_${sufx} -M`" >> $onam;
		  end
	 else
		  set k = "";
		  echo Computing unilateral masks for ${i} ${sufx}
		  foreach k ( `seq 1 1 28` )
				#echo $k;
				set km = ` echo ${k} - .1 | bc -l `;
				# echo $km;
				set kM = ` echo ${k} + .1 | bc -l `;
				# echo $kM;
				# extract the msk_1, etc.
				[ -f msk_${k}_${sufx}.nii.gz ]||\
					 fslmaths $matlas -thr $km -uthr $kM -bin msk_${k}_${sufx};
				echo "${j}, ${age}, ${k}, ${sufx}, `fslstats ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz -k msk_${k}_${sufx} -M`" >> $onam;
		  end
		  #compute bilateral masks
		  echo Computing bilateral masks for ${i} ${sufx}
		  [ -f msk_prefront_${sufx}.nii.gz ]|| fslmaths msk_9_${sufx} -add msk_22_${sufx} -bin msk_prefront_${sufx};
		  echo "${j}, ${age}, prefront, ${sufx}, `fslstats ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz -k msk_prefront_${sufx} -M`" >> $onam;
		  [ -f msk_midfront_${sufx}.nii.gz ]|| fslmaths msk_4_${sufx} -add msk_17_${sufx} -bin msk_midfront_${sufx};
		  echo "${j}, ${age}, midfront, ${sufx}, `fslstats ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz -k msk_midfront_${sufx} -M`" >> $onam;
		  [ -f msk_cingulate_${sufx}.nii.gz ]|| fslmaths msk_7_${sufx} -add msk_20_${sufx} -bin msk_cingulate_${sufx};
		  echo "${j}, ${age}, cingulate, ${sufx}, `fslstats ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz -k msk_cingulate_${sufx} -M`" >> $onam;
		  [ -f msk_medtemp_${sufx}.nii.gz ]|| fslmaths msk_12_${sufx} -add msk_25_${sufx} -bin msk_medtemp_${sufx};
		  echo "${j}, ${age}, medtemp, ${sufx}, `fslstats ${Res}/${i}/dtt2ants/tFraw.1.7tot2_nl.nii.gz -k msk_medtemp_${sufx} -M`" >> $onam;
	 endif
popd > /dev/null; #pushd /nfs/inflammation/mprage/masks/${i} > /dev/null;
end

# ------------------------------------------------------------------
end
end

pushd $cr > /dev/null
