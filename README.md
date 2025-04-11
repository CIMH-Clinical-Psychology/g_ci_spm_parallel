# g_ci_spm_parallel

Please check the [original repo](https://github.com/Fungisai/g_ci_spm) from [@Fungisai](https://github.com/Fungisai) for the original serialised MATLAB version of the library.

This fork is made as a port of CI estimation for the g values in Python.

## Features  
  
- Extract design matrices from SPM `.mat` files.  
- Compute effect sizes (e.g., Hedges' g) from SPM t-maps.  
- Generate confidence intervals for effect sizes.  
- Parallelized computation for large datasets.  
- Integration with neuroimaging tools like `nilearn` for visualization.  
  
---

## Installation  
  
### Prerequisites  

- Python 3.8 or higher  
- Required libraries: `numpy`, `scipy`, `nibabel`, `joblib`, `nilearn`


### Clone the repository  

git clone https://github.com/yourusername/brain_toaster.git  
  
### Navigate to the project directory  

cd brain_toaster  
  
### Install in development mode 

pip install -e .

---

## Usage example

Use by `es_ci_spm_parallel(t_map,con,X,mask_img,confLevel,out_name)` with
```
t_map:      string with the name of a SPM results nifti file with t values (usually something like 'spmT_000X.nii');
con:        contrast vector used to estimate the t map;
X:          filtered and pre-whitened SPM design matrix (This should be SPM.xX.xKXs.X in the respective SPM.mat file. Please note that SPM is mean centering covariates); 
mask_img:   string with the name of the SPM mask file for the analysis (usually 'mask.nii');
confLevel:  confidence level of the estimated condfidence interval (usually something like .90 or .95);
out_name:   string with a prefix to add to the name of the results nifti files 
```

```python
from brain_toaster import es_ci_spm_parallel
from brain_toaster.utils import load_design_matrix, prepare_contrast
from pathlib import Path
import numpy as np

def main():
    # Example usage of the package
    t_map = str(Path('./test_files/masked_spmT_0001.nii'))
    mask_img = str(Path('./test_files/mask.nii'))
    con = np.array([1])
    X = load_design_matrix(str(Path('./test_files/SPM.mat')))
    confLevel = 0.95
    out_name = 'results/exp1'

    es_ci_spm_parallel(t_map, con, X, mask_img, confLevel, out_name)

if __name__ == "__main__":
    main()
```

---

## Visualisation example

```python
plotting.view_img(
    stat_map_img='results/exp1_g.nii', # Input the resulting g calculation
    bg_img=anat_filename,  # Set background image / likely an anatomical file
    cmap='jet',  # Set colormap
    black_bg=False  # Use a black background
)
```

---

## Some notes inherited from the original repository

Please note that the calculations are computationally expensive and can take up to several hours. (Even though it is parallelised now)

If this function is of help for your own work, please cite _Gerchen, M.F., Kirsch, P., & Feld, G.B. (2021) Brain-Wide Inferiority and Equivalence Tests in fMRI Group Analyses: Selected Applications. Human Brain Mapping. DOI: 10.1002/hbm.25664_