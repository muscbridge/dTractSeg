#!/usr/bin/env python
# -*- coding : utf-8 -*-

"""
Functions to make PyDesigner and TractSeg compatible
"""

import os
import os.path as op
import shutil
import glob
import numpy as np
import subprocess
import nibabel as nib

def get_image_spacing(img_path):
    """
    Fetches images spacing
    
    img_path : str
        Path to input volume
        
    Returns
    -------
    str
        Image spacing
        
    See Also
    --------
    applyTransform()
    """
    img = nib.load(img_path)
    affine = img.affine
    return str(abs(round(affine[0, 0], 2)))

def nan_to_zero(input, output, docker=None):
    """
    Convert all NaNs in input volume to zeros

    Parameters
    ----------
    input : str
        Path to input volume
    output : str
        Path to output volume
    docker : str, optional
        Name of docker container to run

    Returns
    -------
    None; writes out file
    """
    if not op.exists(input):
        raise OSError('Input file {} not found'.format(input))
    if not op.exists(op.dirname(output)):
        raise OSError('Directory {} for writing output file does not exist'.format(op.dirname(output)))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'fslmaths',
        input,
        '-nan',
        output
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0: 
        print('Error in converting NaNs to zeros: \{}'.format(error))

def zeroNegative(input, output, docker=None):
    """
    Thresholds all nubers below zero to zero

    Parameters
    ----------
    input : str
        Path to input volume
    output : str
        Path to output volume
    docker : str, optional
        Name of docker container to run

    Returns
    -------
    None; writes out file
    """
    if not op.exists(input):
        raise OSError('Input file {} not found'.format(input))
    if not op.exists(op.dirname(output)):
        raise OSError('Directory {} for writing output file does not exist'.format(op.dirname(output)))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'fslmaths',
        input,
        '-thr', '0',
        output
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0: 
        print('Error in removing negative values: \{}'.format(error))


def createTransform(moving, template, out, omat, docker=None):
    """
    Computes and saves the 4 x 4 affine transformation matrix for transforming
    a moving image to template with 6 DOF rigid-body registration
    
    Parameters
    ----------
    moving : str
        Path to input moving volume
    template : str
        Path to reference volume
    out : str
        Path to output registered volume
    omat : str
        Path to save 4 x 4 affine matrix in .mat extension
    docker : str, optional
        Name of docker container to run
        
    Returns
    -------
    None
        Writes out file
    
    See Also
    --------
    createTransform()

    Examples
    --------
    >>> createTransform(
            moving='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1157/metrics/fa.nii',
            template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
            out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/FA_MNI.nii.gz',
            omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/FA_2_MNI.mat'
        )
        
    Writes both the transformed volume and its 4 x 4 affine matrix
    """
    if not op.exists(moving):
        raise OSError('Moving file {} not found'.format(moving))
    if not op.exists(template):
        raise OSError('Template file {} not found'.format(template))
    if not op.exists(op.dirname(out)):
        raise OSError('Directory {} for writing output file does not exist'.format(op.dirname(out)))
    if not op.exists(op.dirname(omat)):
        raise OSError('Directory {} for writing omat file does not exist'.format(op.dirname(omat)))
    if op.splitext(omat)[-1] != '.mat':
        raise OSError('Affine matrix {} needs to be have .mat extension'.format(omat))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'flirt',
        '-ref', template,
        '-in', moving,
        '-omat', omat,
        '-out', out,
        '-dof', '6',
        '-cost', 'mutualinfo',
        '-searchcost', 'mutualinfo'
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0: 
        print('Error in registering {} to {}: \n{}'.format(moving, template, error))

def applyTransform(moving, template, omat, out, interp='spline', docker=None):
    """
    Transforms moving image to template using a precomputed 4 x 4 affine matrix
    
    Parameters
    ----------
    moving : str
        Path to input moving volume
    template : str
        Path to reference volume
    omat : str
        Path to 4 x 4 affine matrix in .mat extension
    out : str
        Path to save transformed volume
    interp : str, {'spline', 'trilinear', 'nearestneighbour', 'sinc'}, optional
        Interpolation method
    docker : str, optional
        Name of docker container to run
        
    Returns
    -------
    None
        Writes out file
    
    See Also
    --------
    createTransform()
    
    Examples
    --------
    >>> applyTransform(
            moving='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1157/dwi_preprocessed.nii',
            template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
            omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/FA_2_MNI.mat',
            out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/DWI_MNI.nii.gz'
        )
    
    Writes transformed volume with default ``''spline'`` interpolation
    
    >>> applyTransform(
            moving='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1157/brain_mask.nii',
            template='/Users/dataprocessing/Documents/IAM/TractSeg/MNI_FA_template.nii.gz',
            omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/FA_2_MNI.mat',
            out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/brain_mask_MNI.nii.gz',
            interp='nearestneighbour'
        )
    Writes transformed brain mask with ``'nearestneighbours'`` interpolation
    """
    if not op.exists(moving):
        raise OSError('Moving file {} not found'.format(moving))
    if not op.exists(template):
        raise OSError('Template file {} not found'.format(template))
    if not op.exists(omat):
        raise OSError('Omat file {} not found'.format(omat))
    if op.splitext(omat)[-1] != '.mat':
        raise OSError('Affine matrix {} needs to be have .mat extension'.format(omat))
    if not op.exists(op.dirname(out)):
        raise OSError('Directory {} for writing save file does not exist'.format(op.dirname(out)))
    if not any([interp==x for x in ['spline', 'trilinear', 'nearestneighbour', 'sinc']]):
        raise Exception('Invalid interpolation method selection. '
                        'Valid options are "spline", "trilinear", "nearestneighbour", "sinc"')
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    dwi_spacing = get_image_spacing(moving)
    arg = [
        'flirt',
        '-ref', template,
        '-in', moving,
        '-out', out,
        '-init', omat,
        '-dof', '6',
        '-interp', interp,
        '-applyisoxfm', dwi_spacing
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0: 
        print('Error occured applying affine matrix registration of {} to {}: \n{}'.format(moving, template, error))
        
def rotatebvec(bvec, omat, out, docker=None):
    """
    Rotates .bvec file based on 4 x 4 affine transformation matrix
    
    Parameters
    ----------
    bvec : str
        Path to .bvec file
    omat : str
        Path to 4 x 4 affine matrix file
    out : str
        Path to output .bval file
    docker : str, optional
        Name of docker container to run
    
    Return
    ------
    None
        Writes out file
    
    See Also
    --------
    applyTransform()
    
    Examples
    --------
    >>> rotatebvec(
            bvec='/Users/dataprocessing/Documents/IAM/TractSeg_Subs/IAM_1157/dwi_preprocessed.bvec',
            omat='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/FA_2_MNI.mat',
            out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/DWI_MNI.bvec'
        )
        
    Writes rotates bvec file with 4 x 4 affine transformation
    """
    if not op.exists(bvec):
        raise OSError('BVEC fle {} not found'.format(bvec))
    if not op.exists(omat):
        raise OSError('Omat fle {} not found'.format(omat))
    if not op.exists(op.dirname(out)):
        raise OSError('Directory {} for writing save file does not exist'.format(op.dirname(out)))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'rotate_bvecs',
        '-i', bvec,
        '-t', omat,
        '-o', out
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.Popen(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0: 
        print('Unable to rotate bvec file: \n{}'.format(error))

def segBundle(dwi, bval, bvec, out, mask=None, docker=None):
    """
    Creates segmentation of bundles
    
    Parameters
    ----------
    dwi : str
        Path to DWI file
    bval : str
        Path to .bval file accompanying DWI
    bvec : str
        Path to .bvec file accompanying DWI
    out : str
        Path to output directory
    mask : str, optional
        Path to brain mask
    docker : str, optional
        Name of docker container to run
    
    Returns
    -------
    str
        Path to peaks
    str
        Path to bundles
    None
        Writes out files
    
    See Also
    --------
    segStartEnd()
    
    Examples
    --------
    >>> segBundle(
            dwi='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/DWI_MNI.nii.gz',
            bval='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/DWI_MNI.bval',
            bvec='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/DWI_MNI.bvec',
            out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157',
            mask='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/nodif_brain_mask_MNI.nii.gz'
        )
        
    Writes out CSD peaks in output directory, and segmentation bundles in
    ``out/bundle_segmentations``
    """
    if not op.exists(dwi):
        raise OSError('DWI file {} not found'.format(dwi))
    if not op.exists(bval):
        raise OSError('BVAL file {} not found'.format(bval))
    if not op.exists(bvec):
        raise OSError('BVEC file {} not found'.format(bvec))
    if not op.exists(out):
        raise OSError('Output path {} not found'.format(out))
    if mask is not None:
        if not op.exists(mask):
            raise OSError('Brain mask file {} not found'.format(mask))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'TractSeg',
        '-i', dwi,
        '--bvals', bval,
        '--bvecs', bvec,
        '-o', out,
        '--raw_diffusion_input',
        '--output_type', 'tract_segmentation'
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    if mask is not None:
        arg.extend(['--brain_mask', mask])
    p = subprocess.run(arg)
    if p.returncode != 0: 
        print('Unable to run segBundle. '
              'See above for errors')
    return op.join(out, 'peaks.nii.gz'), op.join(out, 'bundle_segmentations')
        
def segStartEnd(peaks, out, docker=None):
    """
    Create segmentation of start and end regions of bundles
    
    Parameters
    ----------
    peaks : str
        Path to CSD peaks volume
    out : str
        Path to output directory
    docker : str, optional
        Name of docker container to run
    
    Returns
    -------
    str
        Path to bundle ending segmentation
    None
        Writes out file
    
    See Also
    --------
    segBundle(), segTracking()
    
    Examples
    --------
    >>> segStartEnd(
        peaks='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/peaks.nii.gz',
        out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157'
        )
        
    Writes out segmentation of start and end regions of bundles in directory
    ``out/endings_segmentation``
    
    
    """
    if not op.exists(peaks):
        raise OSError('Peaks file {} not found'.format(peaks))
    if not op.exists(out):
        raise OSError('Output path {} not found'.format(out))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'TractSeg',
        '-i', peaks,
        '-o', out,
        '--output_type', 'endings_segmentation'
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.run(arg)
    if p.returncode != 0: 
        print('Unable to run segStartEnd. '
              'See above for errors')
    return op.join(out, 'endings_segmentations')

def segCreateTOM(peaks, out, docker=None):
    """
    Create Tract Orientation Maps
    
    Parameters
    ----------
    peaks : str
        Path to CSD peaks volume
    out : str
        Path to output directory
    docker : str, optional
        Name of docker container to run
        
    Returns
    -------
    str
        Path to TOMs
    None
        Writes out file
    
    See Also
    --------
    segBundle(), segStartEnd()
    
    Examples
    --------
    >>> tom_dir = segCreateTOM(
            peaks=peaks_dir,
            out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157',
    )
    
    """
    if not op.exists(peaks):
        raise OSError('Peaks file {} not found'.format(peaks))
    if not op.exists(out):
        raise OSError('Output path {} not found'.format(out))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'TractSeg',
        '-i', peaks,
        '-o', out,
        '--output_type', 'TOM'
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.run(arg)
    if p.returncode != 0: 
        print('Unable to run segCreateTom. '
              'See above for errors')
    return op.join(out, 'TOM')
        
def segTracking(peaks, out, nr_fibers, docker=None):
    """
    Do bundle-specific tracking
    
    Parameters
    ----------
    peaks : str
        Path to CSD peaks volume
    out : str
        Path to output directory
    nr_fibers : int
        Number of fibers to create (default is 2000)
    docker : str, optional
        Name of docker container to run
        
    Returns
    -------
    str
        Path to bundle-specific tracking
    None
        Writes out file
    
    See Also
    --------
    segBundle(), segStartEnd()
    
    Examples
    --------
    >>>segTracking(
        peaks='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157/peaks.nii.gz',
        out='/Users/dataprocessing/Documents/IAM/TractSeg/IAM_1157',
        nr_fibers=1000
        )
    """
    if not op.exists(peaks):
        raise OSError('Peaks file {} not found'.format(peaks))
    if not op.exists(out):
        raise OSError('Output path {} not found'.format(out))
    if not isinstance(nr_fibers, int):
        raise Exception('Please specify number of fibers as an integer')
    if nr_fibers < 0:
        raise Exception('Number of fibers needs to be specified as a positive real integer')
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'Tracking',
        '-i', peaks,
        '-o', out,
        '--nr_fibers', str(nr_fibers)
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.run(arg)
    if p.returncode != 0: 
        print('Unable to run segTracking. '
              'See above for errors')
    return op.join(out, 'TOM_trackings')
    
def segTractometry(metric, tracking_dir, end_dir, out, docker=None):
    """
    Run ROI-based tractometry
    
    Parameters
    ----------
    metric : str
        Path to metric to run tractometry on. Usually an FA map.
    tracking_dir : str
        Path to TOM tracking directory (output of  segTracking)
    end_dir : str
        Path to bundle ending directory
    out : str
        Path to output CSV file
    docker : str, optional
        Name of docker container to run
        
    Returns
    -------
    None
        Writes out file
    """
    if not op.exists(metric):
        raise OSError('Metric file {} not found'.format(metric))
    if not op.exists(tracking_dir):
        raise OSError('TOM tracking directory {} not found'.format(tracking_dir))
    if not op.exists(end_dir):
        raise OSError('Endings segmentation directory {} not found'.format(end_dir))
    if not op.exists(op.dirname(out)):
        raise OSError('Directory {} for writing output CSV file does not exist'.format(op.dirname(out)))
    if op.splitext(out)[-1] != '.csv':
        raise OSError('Output CSV file {} needs to be have .csv extension'.format(out))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    arg = [
        'Tractometry',
        '-s', metric,
        '-i', tracking_dir,
        '-e', end_dir,
        '-o', out
    ]
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    print(' '.join(arg))
    p = subprocess.run(arg)
    if p.returncode != 0: 
        print('Unable to run segTractometry. '
              'See above for errors')

def createAffineFA(dwi, bval, bvec, omat, template, mask=None, docker=None):
    """
    Creates 4 x 4 from TractSeg's FA map for registration of custom scalars 
    
    Parameters
    ----------
    dwi : str
        Path to DWI
    bval : str
        Path to .bval accompanying DWI
    bvec : str
        Path to .bvec accompanying DWI
    omat : path to affine in .mat extension
    mask : str, optional
        Path to brain mask
    docker : str, optional
        Name of docker container to run

    Returns
    -------
    None; writes out files
    """
    if not op.exists(dwi):
        raise OSError('DWI file {} not found'.format(dwi))
    if not op.exists(template):
        raise OSError('Template file {} not found'.format(template))
    if not op.exists(bval):
        raise OSError('BVAL file {} not found'.format(bval))
    if not op.exists(bvec):
        raise OSError('BVEC file {} not found'.format(bvec))
    if not op.exists(op.dirname(omat)):
        raise OSError('Directory {} for writing save file does not exist'.format(op.dirname(omat)))
    if op.splitext(omat)[-1] != '.mat':
        raise OSError('Affine matrix {} needs to be have .mat extension'.format(omat))
    if mask is not None:
        if not op.exists(mask):
            raise OSError('Brain mask file {} not found'.format(mask))
    if docker is not None:
        if not isinstance(docker, str):
            raise Exception('Please provide name of Docker container as a string')
    FA = op.join(op.dirname(omat), 'FA.nii.gz')
    FA_ = op.join(op.dirname(omat), 'FA_MNI.nii.gz')
    arg = [
        'calc_FA',
        '-i', dwi,
        '--bvals', bval,
        '--bvecs', bvec,
        '-o', FA
    ]
    if mask is not None:
        arg.extend(['--brain_mask', mask])
    if docker is not None:
        arg.insert(['docker', 'run', '-it', '--rm', docker])
    p = subprocess.run(arg)
    if p.returncode != 0: 
        print('Unable to run segTractometry. '
              'See above for errors')
    createTransform(
        moving=FA,
        template=template,
        out=FA_,
        omat=omat,
        docker=docker
    )
    os.remove(FA)
    os.remove(FA_)
    
def runtractseg(input, output, template, docker=None):
    """
    Exevutes the entire TractSeg pipeline from start to finish
    
    Parameters
    ----------
    input : str
        Path to subject input folder
    output : str
        Path to subject output folder
    docker : str
        Name of Docker container to run
        
    Returns
    -------
    None
        Writes out files
    """
    subID = op.basename(input)
    print('Processing {}'.format(subID))
    if not op.isdir(output):
        os.makedirs(output, exist_ok=True)

        
    path_fa = op.join(input, 'metrics', 'fa.nii')
    path_fa_nan = op.join(output, 'FA_NO_NAN.nii.gz')
    path_fa_zero = op.join(output, 'FA_NO_NEG.nii.gz')
    path_dwi = op.join(input, 'dwi_preprocessed.nii')
    path_bvec = op.join(input, 'dwi_preprocessed.bvec')
    path_bval = op.join(input, 'dwi_preprocessed.bval')
    path_mask = op.join(input, 'brain_mask.nii')
    print('')
    print('----- Files being processed -----')
    if op.exists(path_fa):
        print('FA: {}'.format(path_fa))
    if op.exists(path_dwi):
        print('DWI: {}'.format(path_dwi))
    if op.exists(path_bvec):
        print('BVEC: {}'.format(path_bvec))
    if op.exists(path_bval):
        print('BVAL: {}'.format(path_bval))
    if op.exists(path_mask):
        print('Mask: {}'.format(path_mask))
    print('---------------------------------')
    print('')
    path_mni_template = template
    path_mni_fa = op.join(output, 'FA_MNI.nii.gz')
    path_mni_dwi = op.join(output, 'DWI_MNI.nii.gz')
    path_mni_dwi_mif = op.join(output, 'DWI_MNI.mif')
    path_mni_bvec = op.join(output, 'DWI_MNI.bvec')
    path_mni_bval = op.join(output, 'DWI_MNI.bval')
    path_mni_mask = op.join(output, 'nodif_brain_mask_MNI.nii.gz')
    path_mni_omat = op.join(output, 'FA_2_MNI.mat')
    
    print('STAGE 1: Image registration into MNI space')
    print('')

    print('Removing NaNs from scalar image')
    nan_to_zero(path_fa, path_fa_nan)

    print('Computing transformation affine matrix')
    createAffineFA(
        dwi = path_dwi,
        bval=path_bval,
        bvec=path_bvec,
        omat=path_mni_omat,
        template=path_mni_template,
        mask=path_mask,
        docker=docker
    )

    print('Transform FA into MNI space...')
    applyTransform(
        moving=path_fa_nan,
        template=path_mni_template,
        omat=path_mni_omat,
        out=path_mni_fa,
        docker=docker
    )
    # Remove negative values from FA
    print('Removing negative values from scalar image')
    zeroNegative(path_mni_fa, path_fa_zero)

    # Remove obsolete files
    os.remove(path_fa_nan)
    os.remove(path_mni_fa)

    # Rename zero-corrected file
    os.rename(path_fa_zero, path_mni_fa)
    
    print('Transform DWI into MNI space...')
    applyTransform(
        moving=path_dwi,
        template=path_mni_template,
        omat=path_mni_omat,
        out=path_mni_dwi,
        docker=docker
    )

    print('Transform brain mask into MNI space...')
    applyTransform(
        moving=path_mask,
        template=path_mni_template,
        omat=path_mni_omat,
        out=path_mni_mask,
        interp='nearestneighbour',
        docker=docker
    )
    
    print('Rotating BVECs into MNI space...')
    rotatebvec(
        bvec=path_bvec,
        omat=path_mni_omat,
        out=path_mni_bvec,
        docker=docker
    )
    
    print('Copying BVALs...')
    shutil.copyfile(path_bval,
                    path_mni_bval
    )
    
    print('')
    print('STAGE 2: Processing with TractSeg')
    print('')
    
    print('Creating segmentation bundles...')
    peaks_dir, bundle_dir = segBundle(
        dwi=path_mni_dwi,
        bval=path_mni_bval,
        bvec=path_mni_bvec,
        out=output,
        mask=path_mni_mask,
        docker=docker
    )
    
    print('Creating enging segmentation bundles...')
    end_dir = segStartEnd(
        peaks=peaks_dir,
        out=output,
        docker=docker
    )
    
    print('Creating TOMs...')
    tom_dir = segCreateTOM(
        peaks=peaks_dir,
        out=output,
        docker=docker
    )
    
    print('Creating TOM tractograms...')
    tracking_dir = segTracking(
        peaks=peaks_dir,
        out=output,
        nr_fibers=1000,
        docker=docker
    )
    
    print('Running tracotometry...')
    segTractometry(
        metric=path_mni_fa,
        tracking_dir=tracking_dir,
        end_dir=end_dir,
        out=op.join(output, 'Tractometry.csv'),
        docker=docker
    )
    
    print('')
    print('')
