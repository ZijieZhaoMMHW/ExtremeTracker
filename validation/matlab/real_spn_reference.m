function real_spn_reference(root)
addpath(fullfile(root,'validation','matlab'));
data_root=getenv('MHW_DATA_DIR');
if isempty(data_root); data_root=root; end
a=load(fullfile(data_root,'mhw_ts.mat'));
c=load(fullfile(data_root,'lonlatoi.mat'));
lon=c.lon_used; lat=c.lat_used;
binary=isfinite(a.mhw_ts(:,:,1:30)) & a.mhw_ts(:,:,1:30)>0;
evalc('alltracks=mhwtrack(binary,lon,lat,4);');
chosen=[];
for k=1:numel(alltracks)
 t=alltracks(k);
 if numel(t.day)>=5 && numel(t.day)<=8 && all(cellfun(@numel,t.xloc)>=20)
  chosen(end+1)=k; %#ok<AGROW>
  if numel(chosen)==3; break; end
 end
end
assert(numel(chosen)>=1,'No suitable real normalization cases found');
tr=alltracks(chosen);
firstday=min(arrayfun(@(s)min(s.day),tr)); lastday=max(arrayfun(@(s)max(s.day),tr));
anomaly=a.mhw_ts(:,:,firstday:lastday);
% Numeric field as provided; no claim that its unknown variable is an SST anomaly.
for norm_flag=0:1
 [radius,snapshots,normalized,xunit,yunit]=spn(tr,lon,lat,anomaly,norm_flag);
 save(fullfile(root,'test','fixtures',sprintf('spn_real_%d.mat',norm_flag)),...
  'tr','lon','lat','anomaly','norm_flag','radius','snapshots','normalized','xunit','yunit','chosen','-v7');
end
fprintf('Actual coordinates, actual tracked shapes and supplied numeric field: %d events.\n',numel(chosen));
end
