function [hF,hB] = imoverlay(B,F,climF,climB,cmap,alpha,x,y,ablated_region,axes)

fontSize = 15;

if isempty(axes)
    axes = gca;
end

% colormap(axes,cmap);

if isempty(climB)
    climB = [min(B(:)), max(B(:))];
end
% To have a grayscale background, replicate image to 3-channels
B = repmat(mat2gray(double(B), double(climB)),[1,1,3]);

% Display the back image
ticksY = round(min(y(:)):50:max(y(:)));
ticksX = round(min(x(:)):50:max(x(:)));

hB = imagesc(axes, x, y, B);axis image;
hold(axes,'on');
hF = imagesc(axes, x, y, F, climF);
set(gca,'YTick',ticksY);
set(gca,'XTick',ticksX);
ylabel('Y (mm)','FontSize', fontSize);
xlabel('X (mm)','FontSize', fontSize);

if ~isempty(ablated_region)
    hold on;
    ah = imagesc(axes, x, y, ablated_region);
    ah.AlphaData = ablated_region;
    axis image
end

% Make the foreground image transparent
alphadata = alpha .* (F >= climF(1));
set(hF,'AlphaData',alphadata);
hold(axes,'off');

% if exist('f')
%     set(f,'Visible','on');
% end
colormap(axes,cmap);
end