function func_output = process_US(trial_name, address, table_vals)
%%%%%

processed_US_file = [address,'US_alldatafor_',trial_name,'.mat'];

if exist(processed_US_file) == 0
    display('Processing ultrasound IQ data');
    set(0,'DefaultAxesFontSize', 15)
    set(0,'DefaultAxesLabelFontSize', 15)
    
    file_name         = 'bufApl0Out_0x0_0x0.data.dm.pmcr';
    az_angle          = table_vals(2);
    el_angle          = table_vals(3);
    frame_rate        = table_vals(3);
    
    r_min             = 0;
    azimuth_angle     = 2*pi*az_angle/360;
    elevation_angle   = 2*pi*el_angle/360;
    theta_min         = -azimuth_angle/2;
    theta_max         = azimuth_angle/2;
    phi_min           = -elevation_angle/2;
    phi_max           = elevation_angle/2;
    pre_abl_frames    = table_vals(5);
    
    sigma             = 3;
    max_decorr        = .1;
    max_IBS           = .1;
    dyn_range          = 3;
    max_Bmode         = 0.08;
    dB_range          = 65;
    
    rho_scale_fact    = 2;
    
    fclose('all');
    
    %%
    full_directory = [address,trial_name];
    data_directory    = dir([full_directory,'/IQDATA*']);
    full_path =strcat(full_directory,'/',{data_directory(1).name},'/',file_name);
    
    Dm_base = read_lbdump(full_path{1});
    which read_lbdump
    r_max = (Dm_base.H.rsz-1)*Dm_base.H.dr;
    inter_frame_time   = 1/frame_rate;
    
    for n = 1:length(data_directory)
        disp(['IQ file number: ',num2str(n)])
        
        full_path =strcat(full_directory,'/',{data_directory(n).name},'/',file_name);
        Dm = read_lbdump(full_path{1});
        
        usData(n) = USDataClass2(Dm.data,12,...
            Dm.Info,r_max,r_min,theta_min,theta_max,phi_min,phi_max,sigma,inter_frame_time);
        
        usData(n).scanConv_Frust();
        
        usData(n).compute3DDecorr();
        
        if n <= pre_abl_frames
            decorr_preablation(:,:,:,n) = usData(n).decorr;
        else
            rawData_cart(:,:,:,:,n-pre_abl_frames) = usData(n).rawData_cart;
            if pre_abl_frames ~= 0
                decorr(:,:,:,n-pre_abl_frames) = usData(n).decorr;
            end
        end
    end
    
    if pre_abl_frames ~= 0
        cum_decorr_preablation = max(decorr_preablation,[],4);
        outside_echo_ind = find(squeeze(rawData_cart(:,:,:,2,1))==0);
        cum_decorr = squeeze(max(decorr,[],4));
        cum_decorr(outside_echo_ind) = NaN;
        
        for i = 1:n-pre_abl_frames
            inst_decorr = squeeze(decorr(:,:,:,i));
            
            decorr_rho = rho_scale_fact*cum_decorr_preablation*...
                sum(inst_decorr.*cum_decorr_preablation,'all')/sum(cum_decorr_preablation.^2,'all');
            
            decorr_motion_corr(:,:,:,i) = (inst_decorr - decorr_rho)./(1 - decorr_rho);
        end
        cum_decorr_motion_corr = max(decorr_motion_corr,[],4);
        cum_decorr_motion_corr(outside_echo_ind) = NaN;
        cum_decorr_motion_corr(cum_decorr_motion_corr<0) = 0;
    end
    
    x_range = usData(1).x_range;
    y_range = usData(1).y_range;
    z_range = usData(1).z_range;
    r_max = usData(1).rmax;
    dr = usData(1).dr;
    
    xyz_range{1} = x_range;
    xyz_range{2} = y_range;
    xyz_range{3} = z_range;
    
    rawData_app = squeeze(rawData_cart(:,:,:,2,:));
    
    disp('Saving data...')
    
    if pre_abl_frames ~= 0
        save([address,'US_alldatafor_',trial_name,'.mat'],'rawData_cart',...
            'decorr','cum_decorr','cum_decorr_motion_corr','xyz_range',...
            'r_max','dr','rho_scale_fact','-v7.3')
    else
        save([address,'US_alldatafor_',trial_name,'.mat'],'rawData_cart',...
            'xyz_range','r_max','dr','rho_scale_fact','-v7.3')
    end
    %     clear rawData_cart ibs decorr decorr_simple decorr_alt cum_decorr cumIBS cum_decorr_simple cum_decorr_alt x_range y_range z_range r_max dr

else
    func_output = {};
    US_file = load(processed_US_file);
    
    rawData_app = squeeze(US_file.rawData_cart(:,:,:,2,:));
    func_output{1} = rawData_app;
    
    xyz_range = US_file.xyz_range;
    func_output{2} = xyz_range;
    
    %     decorr = US_file.decorr;
    if strcmp('LR',trial_name(end-1:end))
        cum_decorr = US_file.cum_decorr;
        func_output{3} = cum_decorr;
        
        cum_decorr_motion_corr = US_file.cum_decorr_motion_corr;
        func_output{4} = cum_decorr_motion_corr;
    end
end
end