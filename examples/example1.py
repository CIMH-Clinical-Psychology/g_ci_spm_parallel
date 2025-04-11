from brain_toaster import es_ci_spm_parallel
from brain_toaster.utils import load_design_matrix, prepare_contrast
from pathlib import Path
import numpy as np

def main():
    # Example usage of the package
    t_map = str(Path('./test_files/spmT_0001.nii'))
    mask_img = str(Path('./test_files/mask.nii'))
    con = np.array([1])
    X = load_design_matrix(str(Path('./test_files/SPM.mat')))
    confLevel = 0.95
    out_name = 'results/exp1'

    es_ci_spm_parallel(t_map, con, X, mask_img, confLevel, out_name)

if __name__ == "__main__":
    main()

