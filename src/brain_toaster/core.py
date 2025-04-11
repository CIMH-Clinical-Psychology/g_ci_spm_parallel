import os
import numpy as np
import nibabel as nib
from .utils import prepare_contrast, compute_effect_size, get_masked_t_values, build_ci_maps, save_nifti_images
from .parallel import compute_ci_for_all

def es_ci_spm_parallel(t_map, con, X, mask_img, confLevel, out_name):
    """
    Estimates effect size g and its CI from SPM t maps.

    Parameters
    ----------
    t_map : str
        Path to the SPM results NIfTI file with t values
    con : array_like
        Contrast vector used to estimate the t map
    X : array_like
        Filtered and pre-whitened SPM design matrix
    mask_img : str
        Path to the SPM mask file
    confLevel : float
        Confidence level (e.g., 0.95)
    out_name : str
        Prefix for output NIfTI files

    Returns
    -------
    dict
        Dictionary containing paths to output files
    """
    # Block 1: Preparations
    con, DoF, td_fac = prepare_contrast(con, X)

    # Block 2: Read t map and compute effect size g
    g, t_values, img_t = compute_effect_size(t_map, td_fac, DoF)

    # Block 3: Load mask image and extract t values within the mask
    ts, mask_indices, Y_m = get_masked_t_values(mask_img, t_values)

    # Block 4: Compute confidence intervals for each voxel in parallel
    # 'ci' will be a (2, number_of_voxels) array: first row lower CI, second row upper CI.
    ci = compute_ci_for_all(ts, td_fac, DoF, confLevel)

    # Block 5: Reconstruct the 3D confidence interval maps from the vectorized results
    ci_l, ci_u = build_ci_maps(Y_m, mask_indices, ci)

    # Block 6: Save the effect size map and CI maps as NIfTI files.
    save_nifti_images(g, ci_l, ci_u, img_t, out_name)
    pass