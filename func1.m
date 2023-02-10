function [registered_vol,compute_registration_flag,load_registration_parameters_flag] = ...
    func1(apex_loc, apex_to_needle_dist, rotation_order, input_rotation_angles, output_path)
%%% Utility function for register.m script
%%% Author: Elmira Ghahramani Z
%%% Last update: 01/26/2022

btn1 = 'No';
btn2 = 'Yes';
defbtn = 'No';
dlgtitle_func1 = 'Re-registration warning';
quest_func1 = ...
    ['The registered file already exists. Do you want to apply registration again using apex loc [',...
    num2str(apex_loc(1:2)), '], needle track slice number ',num2str(apex_loc(3)),...
    ', apex to needle in Z ', num2str(apex_to_needle_dist), ...
    ', order of rotation ', rotation_order, ', and rotation angles [', ...
    num2str(input_rotation_angles), ']?'];

if exist(output_path, 'file')
    answer_func1 = questdlg(quest_func1,dlgtitle_func1,btn1,btn2,defbtn);
    
    if strcmp(answer_func1, 'No')
        registered_vol = load(output_path).registered_vol;
        load_registration_parameters_flag = true;
        compute_registration_flag = false;
    else
        registered_vol = [];
        load_registration_parameters_flag = false;
        compute_registration_flag = true;
    end
else
    registered_vol = [];
    load_registration_parameters_flag = false;
    compute_registration_flag = true;
end
end