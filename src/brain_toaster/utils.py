import numpy as np
from numpy.linalg import matrix_rank, pinv
import nibabel as nib
import scipy.io as sio  

def prepare_contrast(con, X):    
    # Ensure con is a 2D column vector. If con is 1D, convert to shape (n, 1).    
    con = np.atleast_2d(con)    
    if con.shape[0] < con.shape[1]:    
        con = con.T    
  
    # Pad contrast vector with zeros if not long enough.    
    if con.shape[0] < np.atleast_2d(X).shape[0]:    
        padding = np.zeros((X.shape[0] - con.shape[0], 1))    
        con = np.vstack([con, padding])    
  
    # Degrees of freedom    
    DoF = X.shape[0] - matrix_rank(X)    
  
    # Conversion factor from t to d    
    # Fix: Ensure X.T @ X is at least 2D before using pinv  
    XTX = X.T @ X  
    if XTX.ndim < 2:  
        XTX = np.atleast_2d(XTX)  
      
    td_fac = np.sqrt((con.T @ pinv(XTX) @ con).item())    
  
    return con, DoF, td_fac
    
def compute_effect_size(t_map, td_fac, DoF):  
    # Load the t map image (t_map is a filepath)  
    img_t = nib.load(t_map)  
    Y = img_t.get_fdata()        # This holds the t-statistic values  

    # Save a copy of t values for later CI estimation.  
    t_values = Y.copy()  

    # Compute effect size: d = t * td_fac  
    g = Y * td_fac  

    # Apply bias correction factor:  
    J = 1 - (3.0 / (4 * DoF - 1))  
    g = g * J  

    return g, t_values, img_t  

# Example usage:  
# g, t_values, img_t = compute_effect_size(t_map, td_fac, DoF)

def get_masked_t_values(mask_img, t_values):  
    # Load the mask image and get its data  
    img_mask = nib.load(mask_img)  
    Y_m = img_mask.get_fdata()  

    # Get indices where mask is nonzero (or >0)  
    mask_indices = np.where(Y_m > 0)  

    # Extract corresponding t values  
    ts = t_values[mask_indices]  
    return ts, mask_indices, Y_m  

def load_design_matrix(spm_path):  
    """  
    Load the design matrix X from SPM.mat using scipy.io.  

    Parameters  
    ----------  
    spm_path : str  
        Path to the SPM.mat file.  

    Returns  
    -------  
    numpy.ndarray  
        The design matrix X.  
    """  
    try:  
        # Load the SPM.mat file  
        spm_data = sio.loadmat(spm_path, squeeze_me=True, struct_as_record=False)  

        # Extract the design matrix  
        X = spm_data['SPM'].xX.xKXs.X  

        # Ensure it's a numpy array  
        if not isinstance(X, np.ndarray):  
            X = np.array(X)  

        return X  

    except Exception as e:  
        print(f"Error loading SPM.mat file: {e}")  
        raise

# Example usage:  
# ts, mask_indices, Y_m = get_masked_t_values(mask_img, t_values)

def build_ci_maps(Y_m, mask_indices, ci):  
    # Initialize lower and upper CI maps with NaN (same shape as mask)  
    ci_l = np.full(Y_m.shape, np.nan)  
    ci_u = np.full(Y_m.shape, np.nan)  

    # Assign computed CI values to the voxels inside the mask  
    ci_l[mask_indices] = ci[0, :]  
    ci_u[mask_indices] = ci[1, :]  

    return ci_l, ci_u

def save_nifti_images(g, ci_l, ci_u, img_t, out_name):  
    # Use the affine and header from the input t-map  
    out_g = nib.Nifti1Image(g, img_t.affine, img_t.header)  
    out_g_ci_l = nib.Nifti1Image(ci_l, img_t.affine, img_t.header)  
    out_g_ci_u = nib.Nifti1Image(ci_u, img_t.affine, img_t.header)  

    # Save images to disk.  
    nib.save(out_g, out_name + '_g.nii')  
    nib.save(out_g_ci_l, out_name + '_g_ci_l.nii')  
    nib.save(out_g_ci_u, out_name + '_g_ci_u.nii')  

    print("Calculation completed, images saved as:")  
    print(out_name + '_g.nii')  
    print(out_name + '_g_ci_l.nii')  
    print(out_name + '_g_ci_u.nii')