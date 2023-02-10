function func2(registered_vol,output_path)
%%% Utility function for register.m script
%%% Author: Elmira Ghahramani Z
%%% Last update: 01/26/2022

btn1 = 'No';
btn2 = 'Yes';
defbtn = 'No';
quest_func2 = 'The file already exists. Do you want to overwrite?';
dlgtitle_func2 = 'Overwriting warning';

if exist(output_path, 'file')
    answer_func2 = questdlg(quest_func2,dlgtitle_func2,btn1,btn2,defbtn);
    if strcmp(answer_func2, 'Yes')
        save(output_path,'registered_vol');
    end
else
    save(output_path,'registered_vol');
end
end