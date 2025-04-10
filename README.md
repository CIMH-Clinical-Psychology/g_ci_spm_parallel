# g_ci_spm_parallel
Please check the [original repo](https://github.com/Fungisai/g_ci_spm) from [@Fungisai](https://github.com/Fungisai) for the serialised original version of the function.

This fork is made for parallelisation of CI estimation for the g values.
You can either run the function on it's own by calling `es_ci_spm_parallel(t_map,con,X,mask_img,confLevel,out_name)` or running `equiv_test.m` script by making the required changes.
The usage of `equiv_test.m` is not bug-free as of 14.12.2023; but soon it will be :)

### Update on 10.04.2025: 

*Currently `equiv_test.m` works in following ways:*

- **No external masking**: Select no external masking when prompted, then the calculation will resume as usually
- **Multiple external masking**: Select multiple external masks when prompted, the 2nd level spmT files will be masked and consequently individual Hedge's g will be calculated''

### *Future work*:

- **Single external masking**: Select a single external mask when prompted, the 2nd level spmT files will be masked and consequently the Hedge's g will be calculated''

## Matlab function for the estimation of effect size g and its confidence interval from SPM t maps

Use by `es_ci_spm_parallel(t_map,con,X,mask_img,confLevel,out_name)` with
```
t_map:      string with the name of a SPM results nifti file with t values (usually something like 'spmT_000X.nii');
con:        contrast vector used to estimate the t map;
X:          filtered and pre-whitened SPM design matrix (This should be SPM.xX.xKXs.X in the respective SPM.mat file. Please note that SPM is mean centering covariates); 
mask_img:   string with the name of the SPM mask file for the analysis (usually 'mask.nii');
confLevel:  confidence level of the estimated condfidence interval (usually something like .90 or .95);
out_name:   string with a prefix to add to the name of the results nifti files 
```
The results are saved as three nifti files `'out_name_g.nii'` containing the g values map, and `'out_name_g_ci_l.nii'` and `'out_name_g_ci_u.nii'` containing the lower and upper CI limit maps respectively.

Please note that the calculations are computationally expensive and can take up to several hours. (Even though it is parallelised now)

This function is dependent on SPM 12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and the Measures of Effect Size (MES) toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox).

`Examples.zip` contains the code for the simulations in Figure 1 and the second level maps and necessary information to run `g_ci_spm` for the data presented in Figures 2-5 of the manuscript.

If this function is of help for your own work, please cite _Gerchen, M.F., Kirsch, P., & Feld, G.B. (2021) Brain-Wide Inferiority and Equivalence Tests in fMRI Group Analyses: Selected Applications. Human Brain Mapping. DOI: 10.1002/hbm.25664_
