function compute3DDecorr( obj )
%COMPUTE3DDECORR Summary of this function goes here
%   Detailed explanation goes here

% COMPUTE3DDECORR Summary of this function goes here
%   Detailed explanation goes here

%
%%
%     % *Define Guassian Window*
x_range_length = size(obj.x_range,2);
y_range_length = size(obj.y_range,2);
z_range_length = size(obj.z_range,2);
sigx = obj.windowSigma/obj.dx;
sigy = obj.windowSigma/obj.dy;
sigz = obj.windowSigma/obj.dz;

x_mid = ceil(x_range_length/2+1);
y_mid = ceil(y_range_length/2+1);
z_mid = ceil(z_range_length/2+1);

sigfel = x_range_length/(2*pi*sigx);
sigfaz = y_range_length/(2*pi*sigy);
sigfra = z_range_length/(2*pi*sigz);
%     %coeffX =    1/sqrt(2*pi*sigx^2);
%     %coeffY =    1/sqrt(2*pi*sigy^2);
%     %coeffZ =    1/sqrt(2*pi*sigz^2);
% %     xmask = coeffX*exp(-(((1:vol_x_length)-x_mid).^2)/(2*sigx^2));
% %     ymask = coeffY*exp(-(((1:vol_y_length)-y_mid).^2)/(2*sigy^2));
% %     zmask = coeffZ*exp(-(((1:vol_z_length)-z_mid).^2)/(2*sigz^2));
% %     azMask   = filtFactAz .* exp(-((1:nSigPad)-azId).^2/2/sigFAz^2);
% %     raMask   = filtFactRa .* exp(-((1:nRowPad)-rangeId).^2/2/sigFRa^2);
% %     [am,rm]  = meshgrid(azMask,raMask);
% %     maskFilt = fftshift(am.*rm);
xmask = exp(-(((1:x_range_length)-x_mid).^2)/2/sigfel^2);
ymask = exp(-(((1:y_range_length)-y_mid).^2)/2/sigfaz^2);
zmask = exp(-(((1:z_range_length)-z_mid).^2)/2/sigfra^2);
[z_mask_mat,y_mask_mat,x_mask_mat] = ndgrid(zmask,ymask,xmask);
maskfilt = (fftshift(z_mask_mat.*y_mask_mat.*x_mask_mat));
% maskfilt = maskfilt/sum(maskfilt(:));
%
%     [x_mask_mat,y_mask_mat,z_mask_mat] = ndgrid(zmask,xmask,ymask);
%     %maskfilt = fftshift(x_mask_mat.*y_mask_mat.*z_mask_mat);
%     %maskfilt = maskfilt/sum(maskfilt(:));
%     maskfilt = (fftshift(x_mask_mat.*y_mask_mat.*z_mask_mat));
%     maskfilt = maskfilt/sum(maskfilt(:));
%     size(maskfilt)
%     rangee=exp(-2*((1:116)-58).^2*pi^2*9/(116^2));
%     azimuthh=exp(-2*((1:151)-76).^2*pi^2*9/(151^2));
%     elevationn=exp(-2*((1:151)-76).^2*pi^2*9/(151^2));
%     [rr,aa,ee]=ndgrid(rangee,azimuthh,elevationn);
%     maskfilt=fftshift(rr.*aa.*ee); maskfilt=maskfilt/sum(maskfilt(:));
% *compute windowed ibs and autocorr01*
%compute ibs and autocorr before windowing
%     rawData_cart_ROI = obj.rawData_cart.*cat(4,obj.ROILimits,obj.ROILimits);
rawData_cart = obj.rawData_cart;
obj.ibs = abs(rawData_cart).^2;
obj.autocorr01 = rawData_cart(:,:,:,1:(end-1)).*conj(rawData_cart(:,:,:,2:end));

% obj.ibs = abs(obj.rawData_cart).^2;
% obj.autocorr01 = obj.rawData_cart(:,:,:,1:(end-1)).*conj(obj.rawData_cart(:,:,:,2:end));
% set NaN values to small number
obj.autocorr01(find(isnan(obj.autocorr01))) = realmin('double');
obj.ibs(find(isnan(obj.ibs))) = realmin('double');
% compute windowed ibs
for currVolume = 1:size(obj.ibs,4)
    obj.ibs(:,:,:,currVolume) = abs(ifftn(fftn(obj.ibs(:,:,:,currVolume)).*maskfilt));  %IBS Term
end
%compute windowed autcorrelation and decorrelation
for currVolume = 1:(size(obj.ibs,4)-1)
    obj.autocorr01(:,:,:,currVolume) = abs(ifftn(fftn(obj.autocorr01(:,:,:,currVolume)).*maskfilt)); %R01
end
for currVolume = 1:(size(obj.ibs,4)-1)
    %obj.decorr(:,:,:,currVolume) = (1 - abs(obj.autocorr01(:,:,:,currVolume)).^2./(obj.ibs(:,:,:,currVolume).*obj.ibs(:,:,:,currVolume+1)))./obj.interFrameTime;
    R00 = obj.ibs(:,:,:,currVolume);
    R11 = obj.ibs(:,:,:,currVolume+1);
    B2 = R00.*R11; %beta^2
    R01 = abs(obj.autocorr01(:,:,:,currVolume)).^2; %|R01|^2
    tau = obj.interFrameTime*1000; %ms
    
    %B2mean = sum(B2(:).*obj.ROIMask(:))/sum(obj.ROIMask(:)); %%commented by Elmira 01262021
    indInsideEcho = find(obj(1).rawData_cart(:,:,:,1)~=0);
    B2meanInsideEcho = sum(B2(indInsideEcho))/nnz(obj(1).rawData_cart(:,:,:,1));
    
    B2insideEcho=zeros(size(obj(1).rawData_cart(:,:,:,1)));
    B2insideEcho(indInsideEcho)=B2(indInsideEcho);
   
    R01insideEcho=zeros(size(obj(1).rawData_cart(:,:,:,1)));
    R01insideEcho(indInsideEcho)=R01(indInsideEcho);
    % B2mean = sum(B2(:))/nnz(obj(1).rawData_cart(:,:,:,1)); %commented by Elmira 3162021
    
    obj.decorr(:,:,:,currVolume) = 2*(B2insideEcho-R01insideEcho)./...
        (B2insideEcho + B2meanInsideEcho)/tau;
    % obj.decorr(:,:,:,currVolume) = 2*(B2-R01)./(B2 + B2mean)/tau;%commented by Elmira 3162021
    %         obj.decorr(:,:,:,currVolume) = obj.decorr(:,:,:,currVolume).*obj.ROIMask;
    
    obj.decorr_simpleNormxIBS(:,:,:,currVolume) = (B2insideEcho-R01insideEcho)...
./B2insideEcho/tau;
%     obj.decorr_simpleNormxIBS(:,:,:,currVolume) = (B2-R01)./B2/tau;%commented by Elmira 3162021
    %         obj.decorr_simpleNormxIBS(:,:,:,currVolume) = obj.decorr_simpleNormxIBS(:,:,:,currVolume).*obj.ROIMask;
    
    obj.decorr_altNormxIBS(:,:,:,currVolume) = (B2insideEcho-R01insideEcho)...
./B2meanInsideEcho/tau;
%     obj.decorr_altNormxIBS(:,:,:,currVolume) = (B2-R01)./B2mean/tau;%commented by Elmira 3162021
    %         obj.decorr_altNormxIBS(:,:,:,currVolume) = obj.decorr_altNormxIBS(:,:,:,currVolume).*obj.ROIMask;
    
    %         obj.decorr_xIBS(:,:,:,currVolume) = 2*(B2-R01)./(1+B2mean./B2)/tau;
    %         obj.decorr_xIBS(:,:,:,currVolume) = obj.decorr_xIBS(:,:,:,currVolume).*obj.ROIMask;
    %obj.decorr(:,:,:,currVolume) = (1 - abs(obj.autocorr01(:,:,:,currVolume)).^2./(thisMean+obj.ibs(:,:,:,currVolume).*obj.ibs(:,:,:,currVolume+1)))./;
    %obj.decorr(:,:,:,currVolume) = 2*((B2 - (obj.autocorr01)).^2./(obj.interFrameTime*(mean(B2(:)) + B2)));
    
%     obj.BSquared = squeeze(B2);%.*obj.ROIMask; %commented by Elmira 3162021
    obj.BSquared = squeeze(B2insideEcho);%.*obj.ROIMask;
end
% set values outside of volume to small number
obj.autocorr01(find(isnan(obj.rawData_cart(:,:,:,1:(end-1))))) = realmin('double');
obj.ibs(find(isnan(obj.rawData_cart))) = realmin('double');
obj.decorr(find(isnan(obj.rawData_cart(:,:,:,1:(end-1))))) = realmin('double');
% max(obj.decorr(:))
end

