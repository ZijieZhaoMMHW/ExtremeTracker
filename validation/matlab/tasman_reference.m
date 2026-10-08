function tasman_reference(root,data_root)
% Full original Tasman Sea case, using the preserved original MATLAB functions.
if nargin<2; data_root=fullfile(root,'data','tasman_sea'); end
addpath(fullfile(root,'validation','matlab'));
t=load(fullfile(data_root,'tracks_tas.mat'));
c=load(fullfile(data_root,'lonlat_used.mat'));
a=load(fullfile(data_root,'data_anom.mat'));
for norm_flag=0:1
 fprintf('MATLAB Tasman Sea: starting norm_flag=%d\n',norm_flag);
 timer=tic;
 evalc('[radius,snapshots,normalized,xunit,yunit]=spn(t.tracks,c.lon_used,c.lat_used,a.data_anom,norm_flag);');
 elapsed_seconds=toc(timer);
 save(fullfile(root,'validation','results',sprintf('matlab_tasman_spn_%d.mat',norm_flag)),...
  'radius','snapshots','normalized','xunit','yunit','norm_flag','elapsed_seconds','-v7');
 fprintf('MATLAB norm_flag=%d finished in %.6f seconds, radius %.9f km\n',norm_flag,elapsed_seconds,radius);
end
s=load(fullfile(data_root,'snaps.mat')); o=load(fullfile(data_root,'ocean_idx.mat'));
lon=1.25:2.5:360; lat=-88.75:2.5:90; res=2.5;
binary=double(~isnan(s.snaps)); binary(~o.ocean_idx)=NaN;
near=grid_nearby_index(lon,lat,res);
fraction=convzz(binary,near,o.ocean_idx);
smoothed=double(fraction>=0.5);
save(fullfile(root,'validation','results','matlab_tasman_smoothing.mat'),'binary','fraction','smoothed','near','lon','lat','-v7');
fprintf('Original global 5-degree smoothing example finished.\n');
end
