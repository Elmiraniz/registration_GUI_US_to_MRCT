rot_trans_inv = inv(rotation_transform);

temp = rot_trans_inv * [Z(:)' - Zc; X(:)' - Xc; Y(:)' - Yc];
x1 = temp(1,:); 
y1 = temp(2,:); 
z1 = temp(3,:);

x1 = reshape(x1,size(X));
y1 = reshape(y1,size(X));
z1 = reshape(z1,size(X));

ix = round((x1-min(x(:)))/1 + 1);
iy = round((y1-min(y(:)))/1 + 1);
iz = round((z1-min(z(:)))/1 + 1);

tempx = find(ix <= 0 | ix > size(imaging_vol,3));
tempy = find(iy <= 0 | iy > size(imaging_vol,2));
tempz = find(iz <= 0 | iz > size(imaging_vol,1));

tempyxz = [tempy; tempx; tempz];
iy(tempyxz) = 1;
ix(tempyxz) = 1;
iz(tempyxz) = 1;
iyixiz = sub2ind(size(imaging_vol),iz,iy,ix);

us_inMR_subsection = imaging_vol(iyixiz);
