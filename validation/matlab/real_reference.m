function real_reference(root,nums,binary_rule)
if nargin<3; binary_rule='raw'; end
addpath(fullfile(root,'validation','matlab'));
data_root=getenv('MHW_DATA_DIR');
if isempty(data_root); data_root=root; end
a=load(fullfile(data_root,'mhw_ts.mat'));
b=load(fullfile(data_root,'lonlatoi.mat'));
s=struct('size',size(a.mhw_ts),'class',class(a.mhw_ts),'nan_count',nnz(isnan(a.mhw_ts)),...
 'active_count',nnz(a.mhw_ts),'finite_count',nnz(isfinite(a.mhw_ts)),...
 'lon_range',[min(b.lon_used) max(b.lon_used)],'lat_range',[min(b.lat_used) max(b.lat_used)],'nums',nums);
if strcmp(binary_rule,'positive'); a.mhw_ts=isfinite(a.mhw_ts) & a.mhw_ts>0; end
s.binary_rule=binary_rule; s.binary_active_count=nnz(a.mhw_ts);
disp(s); tic;
diary(fullfile(root,'validation','results',['matlab_real_' binary_rule '.log']));
tracks=mhwtrack(a.mhw_ts,b.lon_used,b.lat_used,nums);
s.elapsed_seconds=toc; s.n_tracks=numel(tracks);
save(fullfile(root,'validation','results',['matlab_real_' binary_rule '.mat']),'tracks','s','-v7.3');
diary off; disp(s);
end
