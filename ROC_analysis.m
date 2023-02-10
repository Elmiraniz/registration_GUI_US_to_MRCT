function [AUC_TPR_FPR, AUC_PPV_TPR, ppv, tpr, fpr, thresh, trial_count] = ROC_analysis(trials, input_directory_regist_decorrs_IBS, input_directory_ablation_masks, output_directory, varargin)


if ~strcmp(trials, 'All trials')
    registered_local_decorr = varargin{1}.registered_local_decorr_unscaled;
    registered_local_motcorr_decorr = varargin{1}.registered_local_motcorr_decorr_unscaled;
    
    registered_global_decorr = varargin{1}.registered_global_decorr_unscaled;
    registered_global_motcorr_decorr = varargin{1}.registered_global_motcorr_decorr_unscaled;
    
    registered_combined_decorr = varargin{1}.registered_combined_decorr_unscaled;
    registered_combined_motcorr_decorr = varargin{1}.registered_combined_motcorr_decorr_unscaled;
    
    registered_IBS = varargin{1}.registered_IBS;
    
    mask3D = varargin{2};
    
    inside_echo_ind = find(registered_combined_decorr ~= 0);
    trial_count = 1;
    
    [ppv(:,1), tpr(:,1), fpr(:,1), thresh(:,1)] = ...
        prec_rec(registered_local_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    [ppv(:,2), tpr(:,2), fpr(:,2), thresh(:,2)] = ...
        prec_rec(registered_local_motcorr_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    
    [ppv(:,3), tpr(:,3), fpr(:,3), thresh(:,3)] = ...
        prec_rec(registered_global_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    [ppv(:,4), tpr(:,4), fpr(:,4), thresh(:,4)] = ...
        prec_rec(registered_global_motcorr_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    
    [ppv(:,5), tpr(:,5), fpr(:,5), thresh(:,5)] = ...
        prec_rec(registered_combined_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    [ppv(:,6), tpr(:,6), fpr(:,6), thresh(:,6)] = ...
        prec_rec(registered_combined_motcorr_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    
    [ppv(:,7), tpr(:,7), fpr(:,7), thresh(:,7)] = ...
        prec_rec(registered_IBS(inside_echo_ind), mask3D(inside_echo_ind));

    save_directory = [output_directory,trials,'_allinfo_for_ROCcurve_and_imaging_vectors.mat'];

else
    registered_local_decorr = [];
    registered_local_motcorr_decorr = [];
    
    registered_global_decorr = [];
    registered_global_motcorr_decorr = [];
    
    registered_combined_decorr = [];
    registered_combined_motcorr_decorr = [];
    
    registered_IBS = [];
    mask3D = [];
    
    all_local_decorrs = dir([input_directory_regist_decorrs_IBS, '*local_decorr.mat']);
    all_local_motcorr_decorrs = dir([input_directory_regist_decorrs_IBS, '*local_motcorr_decorr.mat']);
    
    all_global_decorrs = dir([input_directory_regist_decorrs_IBS, '*global_decorr.mat']);
    all_global_motcorr_decorrs = dir([input_directory_regist_decorrs_IBS, '*global_motcorr_decorr.mat']);
    
    all_combined_decorrs = dir([input_directory_regist_decorrs_IBS, '*combined_decorr.mat']);
    all_combined_motcorr_decorrs = dir([input_directory_regist_decorrs_IBS, '*combined_motcorr_decorr.mat']);
    
    all_IBS = dir([input_directory_regist_decorrs_IBS, '*_IBS.mat']);
    
    all_masks = dir([input_directory_ablation_masks, '*ablation_mask3D.mat']);
    
    trial_count = length(all_local_decorrs);
    
    for all_trials_counter = 1:trial_count
        
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_local_decorrs(all_trials_counter).name));
        temp11 = temp1.registered_vol;
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_local_motcorr_decorrs(all_trials_counter).name));
        temp12 = temp1.registered_vol;
        
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_global_decorrs(all_trials_counter).name));
        temp13 = temp1.registered_vol;
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_global_motcorr_decorrs(all_trials_counter).name));
        temp14 = temp1.registered_vol;
        
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_combined_decorrs(all_trials_counter).name));
        temp15 = temp1.registered_vol;
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_combined_motcorr_decorrs(all_trials_counter).name));
        temp16 = temp1.registered_vol;
        
        temp1 = load(strcat(input_directory_regist_decorrs_IBS, all_IBS(all_trials_counter).name));
        temp17 = temp1.registered_vol;
        
        temp2 = load(strcat(input_directory_ablation_masks, all_masks(all_trials_counter).name));
        temp21 = temp2.ablation_mask3D;
        
        inside_echo_ind = find(temp11 ~= 0);
        
        if all_trials_counter == 1
            registered_local_decorr = temp11(inside_echo_ind);
            registered_local_motcorr_decorr = temp12(inside_echo_ind);
            
            registered_global_decorr = temp13(inside_echo_ind);
            registered_global_motcorr_decorr = temp14(inside_echo_ind);
            
            registered_combined_decorr = temp15(inside_echo_ind);
            registered_combined_motcorr_decorr = temp16(inside_echo_ind);
            
            registered_IBS = temp17(inside_echo_ind);
            
            mask3D = temp21(inside_echo_ind);
        else
            registered_local_decorr = cat(1,registered_local_decorr,temp11(inside_echo_ind));
            registered_local_motcorr_decorr = cat(1,registered_local_motcorr_decorr,temp12(inside_echo_ind));
            
            registered_global_decorr = cat(1,registered_global_decorr,temp13(inside_echo_ind));
            registered_global_motcorr_decorr = cat(1,registered_global_motcorr_decorr,temp14(inside_echo_ind));
            
            registered_combined_decorr = cat(1,registered_combined_decorr,temp15(inside_echo_ind));
            registered_combined_motcorr_decorr = cat(1,registered_combined_motcorr_decorr,temp16(inside_echo_ind));
            
            registered_IBS = cat(1,registered_IBS,temp17(inside_echo_ind));
            
            mask3D = cat(1,mask3D,temp21(inside_echo_ind));
        end
        save_directory = [output_directory,trials,num2str(trial_count),'_allinfo_for_ROCcurve_and_imaging_vectors.mat'];
    
    end
    
    [ppv(:,1), tpr(:,1), fpr(:,1), thresh(:,1)] = ...
        prec_rec(registered_local_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    [ppv(:,2), tpr(:,2), fpr(:,2), thresh(:,2)] = ...
        prec_rec(registered_local_motcorr_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    
    [ppv(:,3), tpr(:,3), fpr(:,3), thresh(:,3)] = ...
        prec_rec(registered_global_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    [ppv(:,4), tpr(:,4), fpr(:,4), thresh(:,4)] = ...
        prec_rec(registered_global_motcorr_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    
    [ppv(:,5), tpr(:,5), fpr(:,5), thresh(:,5)] = ...
        prec_rec(registered_combined_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    [ppv(:,6), tpr(:,6), fpr(:,6), thresh(:,6)] = ...
        prec_rec(registered_combined_motcorr_decorr(inside_echo_ind), mask3D(inside_echo_ind));
    
    [ppv(:,7), tpr(:,7), fpr(:,7), thresh(:,7)] = ...
        prec_rec(registered_IBS(inside_echo_ind), mask3D(inside_echo_ind));
    
end

nthresh = size(thresh,1);
AUC_TPR_FPR = sum(abs((fpr(2:nthresh,:)-fpr(1:nthresh-1,:))) ...
    .* (tpr(1:nthresh-1,:)+tpr(2:nthresh,:))/2);

AUC_PPV_TPR = sum(abs((tpr(1:nthresh-1,:)-tpr(2:nthresh,:))) ...
    .* (ppv(1:nthresh-1,:)+ppv(2:nthresh,:))/2);

save(save_directory,...
    'registered_local_decorr', 'registered_local_motcorr_decorr', ...
    'registered_global_decorr', 'registered_global_motcorr_decorr', ...
    'registered_combined_decorr', 'registered_combined_motcorr_decorr', ...
    'registered_IBS', 'mask3D', 'ppv', 'tpr', 'fpr', 'thresh', ...
    'AUC_TPR_FPR', 'AUC_PPV_TPR');

end
