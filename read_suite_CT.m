function vol = read_suite_CT(address)
%%%%%
flip_flag = true;
dicomdir = [address,'_dicom'];
dirinfo = dir(dicomdir);
info = dicominfo([dicomdir '/' dirinfo(end).name]);
info_obj = info;
info_obj.SeriesNumber = 1;
% 
% vol_all = single(squeeze(dicomreadVolume(dicomdir)));
% 
% if length(size(vol_all)) > 3
%     Nim = 5;
%     vol_all = reshape(vol_all,size(vol_all,1),size(vol_all,2), ...
%         Nim,size(vol_all,3)/Nim);
%     vol = squeeze(vol_all(:,:,1,:));
%     info_obj.SeriesNumber = Nim;
% else
%     vol = vol_all;
% end
files = dir(address);
vol = [];
for i = 3:length(files)
    directory = [address,'/',files(i).name];
    vol(:,:,i-2) = rgb2gray(imread(directory));
end

if flip_flag
    vol = flip(vol,3);
end

info_obj.vol = vol;
info_obj.VolumeSlices = size(vol,3);
end