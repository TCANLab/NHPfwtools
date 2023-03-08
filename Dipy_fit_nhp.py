#!/nfs/pkg64/anaconda3/bin/python -i
#
#!/usr/bin/python3
# ------------------------------------------------------------------
## imports
## some imports may be unneccessary or legacy please update to your needs
import sys
import numpy as np
import matplotlib.pyplot as plt
import dipy
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import dipy.reconst.fwdti as fwdti
import dipy.reconst.dti as dti
from dipy.segment.mask import median_otsu
import os.path
import nibabel as nib
import pandas as pd
from numpy.random import randn
from scipy import stats
import matplotlib as mpl
import seaborn as sns

#Use script by running script name followed by the ID of the subject

## specify the location of the processed diffusion data depending on the server being used
## if not running in TCAN Lab just hard code for your own directory structure
if 'amict' == os.uname().nodename:
    drep = os.path.join('/media', os.environ['USER'], 'ldat/DipyFW/dat');
if 'bob' == os.uname().nodename:
    drep = os.path.join('/nfs/inflammation/dti/preprocessed/00_MONTH/Results/Dipy');
if 'bill' == os.uname().nodename:
    drep = os.path.join('/nfs/inflammation/dti/preprocessed/00_MONTH/Results/Dipy');
if 'mabo' == os.uname().nodename:
    drep = os.path.join('/nfs/inflammation/dti/preprocessed/00_MONTH/Results/Dipy');

## pnam is the ID number of the subject and specified in the command line
pnam = sys.argv[1];
print('    Processing ', pnam);
##
# this is the path where all diffusion data and valsvecs reside
dname = os.path.join(drep, pnam);

#specify name of the processed 4D diffusion data
fdwi = os.path.join(dname, 'balldata_eddy_dilated_brain.nii.gz');

#specify the name of the bval file
fbval = os.path.join(dname, 'bvals');

#specify the name of the bvec file
fbvec = os.path.join(dname, 'rotatedbvecs');

bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
gtab = gradient_table(bvals, bvecs)

## load diffusion data
img = nib.load(fdwi);
data = img.get_data();
print("Data shape:", data.shape);
affine = img.affine
print("affine", affine);

## fw model
maskdata, mask = median_otsu(data, vol_idx= [0, 1], median_radius = 3, numpass = 3, autocrop=False, dilate=2);

## quick select all slices 
mask_roi = np.zeros(data.shape[:-1], dtype=bool);
axial_slices = [ i for i in range(data.shape[2]) ];
for i in axial_slices:
    mask_roi[:, :, i] = mask[:, :, i];

# set up file placeholders
## the FW DTI fit 
FA = np.zeros(mask_roi.shape, dtype = float);
F = np.copy(FA);

## model compute, by slice to reduce the memory footprint
# be aware you can choose the fit method of NLS or WLS
for i in axial_slices:
    fwdtimodel = fwdti.FreeWaterTensorModel(gtab, fit_method='WLS');
    fwdtifit = fwdtimodel.fit(maskdata[:, :, i,:], mask=mask_roi[:, :, i]);
    # extract quantities
    F[:, :, i] = np.copy(fwdtifit.f);
    FA[:, :, i] = np.copy(fwdtifit.fa);
    print( "{:2d} ".format(i), end = '..');
print('...fit(s) done');

## eliminate NaNs from the images
F[np.isnan(F)] = 0;
FA[np.isnan(FA)] = 0;

## eliminate voxels potentially contaminated by CSF
tFA = np.copy(FA);
tFA[ F > 0.7 ] = 0;
tF = np.copy(F);
tF[ F > 0.7 ] = 0;

## save the raw and truncated data
nib.save(nib.Nifti1Image(F, affine), os.path.join(dname, 'F.nii.gz'))
nib.save(nib.Nifti1Image(FA, affine), os.path.join(dname, 'FA.nii.gz'))
## 
nib.save(nib.Nifti1Image(tF, affine), os.path.join(dname, 'tF.nii.gz'))
nib.save(nib.Nifti1Image(tFA, affine), os.path.join(dname, 'tFA.nii.gz'))
