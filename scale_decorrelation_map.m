function vol_scaled = scale_decorrelation_map(vol, max_decorr, dyn_range)

vol_scaled = vol / max_decorr;
vol_scaled(find(vol_scaled<=0)) = 1.e-16;
vol_scaled = (log10(vol_scaled) + dyn_range) / dyn_range;
vol_scaled(find(vol_scaled>1)) = 1;
vol_scaled(find(vol_scaled<0)) = 0;
vol_scaled(isnan(vol_scaled)) = 0;

end