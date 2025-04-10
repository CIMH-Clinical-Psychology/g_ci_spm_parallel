function pool = es_ci_spm_parallel(t_map,con,X,mask_img,confLevel,out_name, pool)
%es_ci_spm.m Estimates effect size g and its CI from SPM t maps
%
%   Use by es_ci_spm(t_map,con,X,mask_img,confLevel,out_name) with
%   t_map:      string with the name of a SPM results nifti file with t values (usually something like 'spmT_000X.nii')
%   con:        contrast vector used to estimate the t map
%   X:          filtered and pre-whitened SPM design matrix (This should be SPM.xX.xKXs.X in the respective SPM.mat file. Please note that SPM is mean centering covariates) 
%   mask_img:   string with the name of the SPM mask file for the analysis (usually 'mask.nii')
%   confLevel:  confidence level of the estimated condfidence interval (usually something like .90 or .95)
%   out_name:   string with a prefix to add to the name of the results nifti files 
%
%   The results are saved as three nifti files:[out_name '_g.nii']
%   containing the g values map, and [out_name '_g_ci_l.nii'] and
%   [out_name '_g_ci_u.nii'] containing the lower and upper CI limit maps
%   respectively.
%
% Please note that the calculations are computationally expensive and can take up to several
% hours.
%
% This function is dependent on SPM 12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
% and the Measures of Effect Size (MES) toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox).
%
% This work is is distributed under the terms of 
% the GNU General Public Licence as published by the Free Software Foundation
% (either version 3, or at your option, any later version).
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% If this file is of help for your own work, please cite Gerchen, M.F., Kirsch, P., & Feld, G.B. (2021) Brain-Wide Inferiority 
% and Equivalence Tests in fMRI Group Analyses: Selected Applications. Human Brain Mapping. DOI: 10.1002/hbm.25664
%
% Copyright (C) 2021
%
% Martin Fungisai Gerchen & Gordon Feld
% Department of Clinical Psychology
% Medical Faculty Mannheim/Heidelberg University 
% Central Institute of Mental Health
% Mannheim, Germany
%% preparations

% make contrast a column vector, if not already
if size(con,2)>size(con,1)
    con=con';
end

% pad contrast vector with zeros if not of sufficient length
if size(con,1) < size(X,2)
    con = [con; zeros(size(X,2)-size(con,1),1)];
end

% estimate degrees of freedom
DoF = size(X,1)-rank(X);

% estimate conversion factor from t to d
td_fac = sqrt(con'*pinv(X'*X)*con);

%% estimate effect size g
% load t map
V=spm_vol(t_map);
Y=spm_read_vols(V);

% keep t values for CI estimation 
t=Y;

% d
Y = Y.*td_fac;

% bias correction
J=(1-(3./(4*DoF-1))); % bias correction factor
Y = Y.*J;

%% estimate confidence interval for g 
    
% load mask file to extract voxel values
V_m=spm_vol(mask_img);
Y_m=spm_read_vols(V_m);

ts = zeros(length(find(Y_m)),1); % vector of ts
ts(:,1) = t(logical(Y_m));

% estimate exact CI
ci = zeros(2,numel(ts));
%% parallel loop checks
%% Arg check
switch nargin
    case 6
        %% Check if parpool is on
        myCluster = parcluster('local');
        % if the num of ptcps are more than available cores, adjust
        if numel(ts) > myCluster.NumWorkers
            pool = parpool(myCluster.NumWorkers);
        else
            pool = parpool(length(ptcp_array));
        end
end
% %% This block is to activate the parallel pooling
% % Check if parpool is on
% myCluster = parcluster('local');
% % if the num of ptcps are more than available cores, adjust
% if numel(ts) > myCluster.NumWorkers
%     pool = parpool('local', myCluster.NumWorkers-2);
% else
%     pool = parpool('local', numel(ts));
% end
%%
% loop over t values and estimate CI for g
% for i=1:numel(ts)
% do it with parallel pooling, it does require a special logic of dividing
% the numel into unique parts
total_iterations = 25;
inds = floor(1:numel(ts)/total_iterations:numel(ts)); % floor here so always integer indices
%% parallel loop
try
%     parfor ind = 1:25
%         license('test', 'Statistics_Toolbox', 'enable')
%         if ind == 1
%             rng = inds(ind):inds(ind+1);
%         elseif ind == 25
%             rng = inds(ind)+1:numel(ts);
%         else
%             rng = inds(ind)+1:inds(ind+1);
%         end
%         aux_ci = zeros(2,numel(rng)); % create the auxillary matrix
%         for i = rng
%             t_tmp = ts(i);
%             ci_tmp = ncpci(t_tmp,'t',DoF,'confLevel',confLevel)'*td_fac; % uses the 'ncpci.m' function from the MES toolbox
%             aux_ci(:,i) = ci_tmp;
%         end
% %         ci(:,rng(1):rng(end)) = aux_ci; %
%         % Send update - adjust frequency if needed (e.g., every 1% or 5%)  
%         if mod(i, max(5, floor(ind/100))) == 0  % Update every 5%  
%             send(D, 1);  
%             p = i;  
%         end
%         send(progress, i);
%     end
%     parfor ind = 1:total_iterations
%         license('test', 'Statistics_Toolbox', 'enable');
%         
%         % Define the range for this iteration
%         if ind == 1
%             rng = inds(ind):inds(ind+1);
%         elseif ind == total_iterations
%             rng = inds(ind)+1:numel(ts);
%         else
%             rng = inds(ind)+1:inds(ind+1);
%         end
%         
%         % Create the auxiliary matrix
%         aux_ci = zeros(2, numel(rng));
%         
%         % Process each element in the range
%         for idx = 1:numel(rng)  % Use local indexing
%             i = rng(idx);  % Get the actual index from rng
%             t_tmp = ts(i);
%             %             ci_tmp = ncpci(t_tmp,'t',DoF,'confLevel',confLevel)'*td_fac;
%             %             aux_ci(:,idx) = ci_tmp;  % Use idx instead of i
%             % Check if t_tmp is zero or very close to zero
%             if t_tmp == 0 || abs(t_tmp) < eps
%                 aux_ci(:,idx) = [0; 0];  % if t=0
%             else
%                 ci_tmp = ncpci(t_tmp,'t',DoF,'confLevel',confLevel)'*td_fac;
%                 aux_ci(:,idx) = ci_tmp;
%             end
%         end
%     end
% Preallocate a cell array to store results from each iteration
    ci_parts = cell(total_iterations, 1);

    parfor ind = 1:total_iterations
        license('test', 'Statistics_Toolbox', 'enable');

        % Define the range for this iteration
        if ind == 1
            rng = inds(ind):inds(ind+1);
        elseif ind == total_iterations
            rng = inds(ind)+1:numel(ts);
        else
            rng = inds(ind)+1:inds(ind+1);
        end

        % Create the auxiliary matrix
        aux_ci = zeros(2, numel(rng));

        % Process each element in the range
        for idx = 1:numel(rng)  % Use local indexing
            i = rng(idx);  % Get the actual index from rng
            t_tmp = ts(i);

            % Check if t_tmp is zero or very close to zero
            if t_tmp == 0 || abs(t_tmp) < eps
                aux_ci(:, idx) = [0; 0];  % if t=0
            else
                ci_tmp = ncpci(t_tmp, 't', DoF, 'confLevel', confLevel)' * td_fac;
                aux_ci(:, idx) = ci_tmp;
            end
        end

        % Store the results for this iteration in the cell array
        ci_parts{ind} = {rng, aux_ci};
    end

    % Combine the results from all iterations
    for ind = 1:total_iterations
        rng = ci_parts{ind}{1};
        aux_ci = ci_parts{ind}{2};
        ci(:, rng) = aux_ci;
    end
catch err_fold
    % Display detailed error information  
    disp('Error occurred:');  
    disp(err_fold.message);  
      
    % Display stack trace with line numbers  
    for i = 1:length(err_fold.stack)  
        disp(['File: ' err_fold.stack(i).file]);  
        disp(['Function: ' err_fold.stack(i).name]);  
        disp(['Line: ' num2str(err_fold.stack(i).line)]);  
        disp('---');  
    end  
      
    % use rethrow to see the error in standard format  
    % rethrow(err_fold);  
      
%     % End pool  
%     delete(gcp)
end

% put lower and upper CI limit into 3D maps
ci_l = nan(size(Y_m));
ci_l(logical(Y_m)) = ci(1,:);

ci_u = nan(size(Y_m));
ci_u(logical(Y_m)) = ci(2,:);    

%% save images
fprintf('\n====================================\nCalculation completed, saving images:\n--> %s\n====================================\n', [out_name '_g.nii']);
V_out=V;

% g map
V_out.fname=[out_name '_g.nii'];
spm_write_vol(V_out,Y);

% g CI lower limit
V_out.fname=[out_name '_g_ci_l.nii'];
spm_write_vol(V_out,ci_l);

% g CI upper limit
V_out.fname=[out_name '_g_ci_u.nii'];
spm_write_vol(V_out,ci_u);

% %% End parallel pool
% delete(gcp)
end

