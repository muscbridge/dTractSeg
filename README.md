# dTractSeg
A TractSeg pipeline compatible with PyDesigner outputs

## Introduction
PyDesigner is an open-source hands-off DTI/DKI preprocessing CLI provided by MUSC.
It was designed specifically for estimating accurate tensors and scalar maps.
TractSeg, on the other hand, is a white matter bundle segmentation and tractometry
analysis pipeline.

This Python scripts package prepares bridges the two pipelines together to allow 
tractometry on PyDesigner's scalar maps.

## Installation
This package has the same requirements as TractSeg and PyDesigner:

1. [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
2. [MRtrix3](https://www.mrtrix.org/)
3. [PyTorch](https://pytorch.org/)
4. Python 3

### Install TractSeg
Install TractSeg using the instructions posted on their page. The package is also provided in the ``TractSeg`` folder within this directory, that can be installed after ``cd`` with:

```
pip insall .
```

## Running the Pipeline
The TractSeg pipeline can be run in one of the following two ways:

### Option 1: Run with one function
This function will run all the steps highlighted in 2. in one single function. Attempt running
this first, then move to 2. if this fails.

Import the ``tscompatibility`` module with:

```
import tscompatibility as ts
```

Execute the  ``ts.runtractseg`` function with:

```
ts.runtractseg(input='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158',
            output='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158')
```

Where ``input`` is the path to PyDesigner output directory, and ``output`` is the directory where
TractSeg outputs are saved.


### Option 2: Run with individual functions
This is the recommended way to run the pipeline as it is prone to less errors, while providing refined
control over the analysis. The subject being processed here is ``IAM_1158``

Import the ``tscompatibility`` module with:

```
import tscompatibility as ts
```

#### Remove NaNs from scalar image:

FSL functions like ``flirt`` do not work if there are NaNs present in input
volumes, which is definitely the case in PyDesigner outputs. Start by removing
NaNs on any scalar (FA in this example) with:

```
ts.nan_to_zero(
    input='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/metrics/fa.nii'
    output='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_NO_NAN.nii.gz'
)
```

#### Compute transformation affine matrix:

Compute the transformation affine matrix to bring all scalar map and DWI into MNI space with.
This method computes a temporary FA map that is used to register to an MNI image. TractSeg
recommends registering to their ``MNI_FA_template.nii.gz`` template map.

```
print('Computing transformation affine matrix')
createAffineFA(
    dwi = '/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/dwi_preprocessed.nii',
    bval='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/dwi_preprocessed.bval',
    bvec='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/dwi_preprocessed.bvec',
    omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_2_MNI.mat',
    template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
    mask='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/brain_mask.nii'
)
```

The argument ``omat`` defines the path to write affine matrix; please specify this with a
``.mat`` extension

#### Transform volumes into MNI space:

Using the affine transformation matrix file, transform the DWI, scalar map, and brain mask with:

```
# Transform FA without NaNs
print('Transform FA into MNI space...')
applyTransform(
    moving='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_NO_NAN.nii.gz',
    template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
    omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_2_MNI.mat',
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_MNI.nii.gz'
)
# Remove FA_NO_NAN.nii.gz because it is not needed anymore
os.remove(path_fa_nan)

# Transform DWI
print('Transform DWI into MNI space...')
applyTransform(
    moving=/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/dwi_preprocessed.nii,
    template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
    omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_2_MNI.mat',
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/DWI_MNI.nii.gz'
)

# Transform brain mask
print('Transform brain mask into MNI space...')
applyTransform(
    moving=/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/brain_mask.nii,
    template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
    omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_2_MNI.mat',
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/nodif_brain_mask_MNI.nii.gz'
    interp='nearestneighbour',
    docker=docker
)
```

Transforming the he brain mask requires ``'nearestneighbour'`` interpolation to
retain its binary composition.

#### Rotate BVEC and copy BVAL:

The gradient vectors (BVECs) have to be rotated accordingly to represent the transformed
DWI. BVALs also need to be copied over and both these steps can be performed with:

```
# Rotate BVEC
print('Rotating BVECs into MNI space...')
rotatebvec(
    bvec='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/dwi_preprocessed.bvec',
    omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_2_MNI.mat',
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/DWI_MNI.bvec',
)

# Copy BVAL
print('Copying BVALs...')
shutil.copyfile('/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1158/dwi_preprocessed.bvec',
                '/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/DWI_MNI.bval'
)
```

**Note**: run ``import shutil`` to import to be able to use copy command

#### Create segmentation bundles
Begin creating segmentation bundles using TractSeg wrappers with the command:

```
print('Creating segmentation bundles...')
peaks_dir, bundle_dir = segBundle(
    dwi='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/DWI_MNI.nii.gz',
    bval=/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/DWI_MNI.bval,
    bvec='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/DWI_MNI.bvec',
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158',
    mask='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/nodif_brain_mask_MNI.nii.gz
)
```
 This creates two variables ``peaks_dir`` and ``bundle_dir``, to provide the absolute
 paths to ``peaks.nii.gz`` and ``bundle_segmentations`` directories within the output
 folder.

**Note**: execute this on transformed DWI (MNI space)

#### Create ending segmentation of bundles
Next, create endings segmentation with:

```
print('Creating enging segmentation bundles...')
end_dir = segStartEnd(
    peaks=peaks_dir,
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158',
)
```

The variable ``end_dir`` holds the absolute path to ``endings_segmentations`` folder
within the output directory.

#### Create tract orientation maps (TOMs)
Once ending segmentation is done, perform TOM with:

```
print('Creating TOMs...')
tom_dir = segCreateTOM(
    peaks=peaks_dir,
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158'
)
```

The variable ``tom_dir`` holds the absolute path to ``TOM`` folder within
the output directory.

#### Create TOM tractograms
Next, create TOM tractograms with:

```
print('Creating TOM tractograms...')
tracking_dir = segTracking(
    peaks=peaks_dir,
    out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158',
    nr_fibers=5000
)
```

The varible ``tracking_dir`` holds the absolute path to ``TOM_trackings`` folder
withing the output directory.

#### Perform tractometry
Finally, run tractometry on MNI-space scalar (``FA.nii.gz)`` with:

```
print('Running tracotometry...')
segTractometry(
    metric='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/FA_MNI.nii.gz',
    tracking_dir=tracking_dir,
    end_dir=end_dir,
    out=/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1158/Tractometry.csv
)
```

This will write the CSV file ``Tractometry.csv`` containing fiber node values per ROI.

### Batch Processing
Users may run an entire collection of subjects easily with option 1. Refer to
the Jupyter Notebook file ``TractSeg_Pipeline.ipynb`` in this repo to see an example.

## Group Based Analysis
Running the entire pipeline produces all outputs required for group-based analysis.
Please refer to TractSeg's documentation view the instructions for this.
