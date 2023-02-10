function [image_vol, vol] = read_1day_post_ablation_image_vol(address)
%%%
% Author: Elmira Ghahramani Z.
% Last update: 01/31/2023
%%%

flip_flag = true;

files = dir(address);
vol = [];
n_junk = 3;
for i = n_junk+1:length(files)
    directory = [address,'/',files(i).name];
    vol(:,:,i-n_junk,:) = imread(directory);
end

if flip_flag
    vol = flip(vol,3);
end

image_vol = permute(vol,[1,2,4,3]);
vol = single(0.2989 * vol(:,:,:,1) + ...
             0.5870 * vol(:,:,:,2) + ...
             0.1140 * vol(:,:,:,3));

end