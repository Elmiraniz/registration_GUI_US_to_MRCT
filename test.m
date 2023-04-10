rot_trans_inv = inv(rotation_transform);

ZXYtranslated = [(Z(:)-Zc)'; (X(:)-Xc)'; (Y(:)-Yc)'];
ZXYtranslatedUnrotated = rot_trans_inv * ZXYtranslated;

Zunrotated = ZXYtranslatedUnrotated(1,:);
Zunrotated = reshape(Zunrotated,size(X));

Xunrotated = ZXYtranslatedUnrotated(2,:);
Xunrotated = reshape(Xunrotated,size(X));

Yunrotated = ZXYtranslatedUnrotated(3,:);
Yunrotated = reshape(Yunrotated,size(X));

tempCenter = rot_trans_inv * [Zc; Xc; Yc];

Xunrotated_US_half_ind = ceil(Xunrotated_US/2);
Zunrotated_US_half_ind = ceil(Zunrotated_US/2);

if Xunrotated_US_half_ind >= XApex_ind
    ind1_US = Xunrotated_US_half_ind - XApex_ind + 1;
else
    ind1_US = 1;
end
if Yunrotated_US_half_ind >= YApex_ind
    ind3_US = Yunrotated_US_half_ind - YApex_ind + 1;
else
    ind3_US = 1;
end
if Zunrotated_US_half_ind >= ZApex_ind
    ind5_US = Zunrotated_US_half_ind - ZApex_ind + 1;
else
    ind5_US = 1;
end 


maxIndX_MR = size(X,2);
maxIndY_MR = size(X,1);
maxIndZ_MR = size(X,3);

if Xunrotated_US_half_ind >= maxIndX_MR-XApex_ind
    ind2_US = Xunrotated_US_half_ind + maxIndX_MR-XApex_ind;
else
    ind2_US = 1;
end
if Yunrotated_US_half_ind >= YApex_ind
    ind4_US = Yunrotated_US_half_ind - YApex_ind + 1;
else
    ind5_US = 1;
end
if Zunrotated_US_half_ind >= ZApex_ind
    ind6_US = Zunrotated_US_half_ind - ZApex_ind + 1;
else
    ind6_US = 1;
end 
% % 
% % Zunrotated_US = tempCenter(1) + x;
% % Xunrotated_US = tempCenter(2) + y;
% % Yunrotated_US = tempCenter(3) + z;
% 
US_inMRsize = imresize3(imaging_vol,[size(imaging_vol,1)/Imdy, ...
    size(imaging_vol,2)/Imdx, size(imaging_vol,3)/Imdz]);
% 
% if min(Zunrotated_US(:)) <= min(Zunrotated(:))
%     minZunrotated = min(Zunrotated(:)); 
% else
%     minZunrotated = min(Zunrotated_US(:)); 
% end
% 
% if min(Xunrotated_US(:)) <= min(Xunrotated(:))
%     minXunrotated = min(Xunrotated(:)); 
% else
%     minXunrotated = min(Xunrotated_US(:)); 
% end
% 
% if min(Yunrotated_US(:)) <= min(Yunrotated(:))
%     minYunrotated = min(Yunrotated(:)); 
% else
%     minYunrotated = min(Yunrotated_US(:)); 
% end
% 
% 
% if max(Zunrotated_US(:)) >= max(Zunrotated(:))
%     maxZunrotated = max(Zunrotated(:)); 
% else
%     maxZunrotated = max(Zunrotated_US(:)); 
% end
% 
% if max(Xunrotated_US(:)) >= max(Xunrotated(:))
%     maxXunrotated = max(Xunrotated(:)); 
% else
%     maxXunrotated = max(Xunrotated_US(:)); 
% end
% 
% if max(Yunrotated_US(:)) >= max(Yunrotated(:))
%     maxYunrotated = max(Yunrotated(:)); 
% else
%     maxYunrotated = max(Yunrotated_US(:)); 
% end
% % finding index of unrotated US boundary
% X_US_Ind1 = find(abs(Xunrotated_US(1,:,1)-minXunrotated)<Imdx);
% X_US_Ind2 = find(abs(Xunrotated_US(1,:,1)-maxXunrotated)<Imdx);
% Y_US_Ind1 = find(abs(Yunrotated_US(:,1,1)-minYunrotated)<Imdy);
% Y_US_Ind2 = find(abs(Yunrotated_US(:,1,1)-maxYunrotated)<Imdy);
% Z_US_Ind1 = find(abs(Zunrotated_US(1,1,:)-minZunrotated)<Imdz);
% Z_US_Ind2 = find(abs(Zunrotated_US(1,1,:)-maxZunrotated)<Imdz);

%%
p = find(imaging_vol ~= -inf); % points for valid interpolation
    
    %image_initialization_time = toc
    
    %tic;
    % for each point to interpolate, find nearest previous neighbors, and
    %  distances from them along R, mu, nu directions
    imu = zeros(size(X));
    Lmu = zeros(size(X));
    imu(p) = floor((mu0(p) - mumin)/dmu) + 1;
    Lmu(p) = mu0(p) - muvec(imu(p))';
    
    inu = zeros(size(X));
    Lnu = zeros(size(X));
    inu(p) = floor((nu0(p) - numin)/dnu) + 1;
    Lnu(p) = nu0(p) - nuvec(inu(p))';
    
    iR = zeros(size(X));
    LR = zeros(size(X));
    iR(p) = floor((R0(p) - obj.rmin)/obj.dr) + 1;
    LR(p) = R0(p) - Rvec(iR(p))';

    %interpolation_initialization_time = toc
    %tic;
    Icart = zeros(size(R0));
    for ip = 1:length(p)
        q = p(ip);
        
        % differences to be used below, trying to save a few flops
        drmLR = obj.dr-LR(q);
        dmumLmu = dmu-Lmu(q);
        dnumLnu = dnu-Lnu(q);
        Icart(q) = Isph(iR(q),imu(q),inu(q)) ...
            * drmLR * dmumLmu * dnumLnu + ...
            Isph(iR(q)+1,imu(q),inu(q)) ...
            * LR(q) * dmumLmu * dnumLnu + ...
            Isph(iR(q),imu(q)+1,inu(q)) ...
            * drmLR * Lmu(q) * dnumLnu + ...
            Isph(iR(q),imu(q),inu(q)+1) ...
            * drmLR * dmumLmu * Lnu(q) + ...
            Isph(iR(q)+1,imu(q),inu(q)+1) ...
            * LR(q) * dmumLmu * Lnu(q) + ...
            Isph(iR(q),imu(q)+1,inu(q)+1) ...
            * drmLR * Lmu(q) * Lnu(q) + ...
            Isph(iR(q)+1,imu(q)+1,inu(q)) ...
            * LR(q) * Lmu(q) * dnumLnu + ...
            Isph(iR(q)+1,imu(q)+1,inu(q)+1) ...
            * LR(q) * Lmu(q) * Lnu(q);
    end
    obj.rawData_cart(:,:,:,volIndex) = Icart/(obj.dr*dmu*dnu);
    %interpolation_time = toc
    end   
%%
yaw = -pi/4;
pitch = 0;
roll = 0;
Rx = [1 0 0; 0 cos(yaw) -sin(yaw); 0 sin(yaw) cos(yaw)];
    Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
    Rz = [cos(roll) -sin(roll) 0; sin(roll) cos(roll) 0; 0 0 1];
switch rotation_order
        case rotation_order_cases{1}
            rotation_transform = Rz*Ry*Rx;
        case rotation_order_cases{2}
            rotation_transform = Ry*Rz*Rx;
        case rotation_order_cases{3}
            rotation_transform = Rz*Rx*Ry;
        case rotation_order_cases{4}
            rotation_transform = Rx*Rz*Ry;
        case rotation_order_cases{5}
            rotation_transform = Ry*Rx*Rz;
        case rotation_order_cases{6}
            rotation_transform = Rx*Ry*Rz;
end
    
p = find(imaging_vol ~= -inf);
rotatedxyz = rotation_transform*[x(p)';y(p)';z(p)'];
xr = rotatedxyz(1,:);
yr = rotatedxyz(2,:);
zr = rotatedxyz(3,:);

rot_US = zeros(size(imaging_vol));
for i = 1:length(p)
    temp = find(xr - x(p(i)) <= 1 & yr - y(p(i)) <= 1 & zr - z(p(i)) <= 1);
    if ~isempty(temp)
        rot_US(p(i)) = imaging_vol(temp(1));
    end
end
    
    