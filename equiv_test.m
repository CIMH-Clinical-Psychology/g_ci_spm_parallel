%% Brain-wide inferiority and equivalence tests
% With using functions generated by Gerchen et al. 2021

%% First load the necessary tools
% might change the locations in the future
gCI_path = uigetdir(title='Select the directory where the estimation of effect size g and its confidence interval toolbox is located');
while gCI_path == 0
    gCI_path = uigetdir(title='Select the directory where the estimation of effect size g and its confidence interval toolbox is located');
end
mes_path = uigetdir(title='Select the directory where the measures of effect size toolbox is located');
while mes_path == 0
    mes_path = uigetdir(title='Select the directory where the measures of effect size toolbox is located');
end

addpath(gCI_path) % g-CI toolbox
addpath(mes_path) % MES toolbox

%% Now run the function
% with using es_ci_spm(t_map,con,X,mask_img,confLevel,out_name)
% t_map: string with the name of a SPM results nifti file with t values (usually something like 'spmT_000X.nii'); 
% con: contrast vector used to estimate the t map; 
% X: filtered and pre-whitened SPM design matrix (This should be SPM.xX.xKXs.X in the respective SPM.mat file. Please note that SPM is mean centering covariates); 
% mask_img: string with the name of the SPM mask file for the analysis (usually 'mask.nii'); 
% confLevel: confidence level of the estimated condfidence interval (usually something like .90 or .95); 
% out_name: string with a prefix to add to the name of the results nifti files

%% In a loop, run the function
% but first set the directories
main_dir = '/zi/flstorage/group_klips/data/data/Cagatay/TEDDS_OCRT2/OCRT2';
% main_folds = dir([main_dir filesep 'mult_second*']); 
main_folds = dir([main_dir filesep 'mult_sep_secondlevels']);
main_folds = strcat({main_folds.folder},filesep,{main_folds.name})';
% loop over them
for q = 1:numel(main_folds)
    if ~contains(main_folds{q}, '/.')s
        cont_folds = dir([main_folds{q} filesep 'Contrast*']);
        cont_folds = strcat({cont_folds.folder},filesep,{cont_folds.name})';
        for qq = 1:numel(cont_folds)
            t_map = [cont_folds{qq} filesep 'spmT_0001.nii'];
            con = [1];
            SPM = load([cont_folds{qq} filesep 'SPM.mat']);
            X = SPM.SPM.xX.xKXs.X;
            mask_img = [cont_folds{qq} filesep 'mask.nii'];
            confLevel = .95;
            out_name = [cont_folds{qq} filesep];
            % es_ci_spm(t_map, con, X, mask_img, confLevel, out_name)
            es_ci_spm_edited(t_map, con, X, mask_img, confLevel, out_name, pool)
        end
    end
end


