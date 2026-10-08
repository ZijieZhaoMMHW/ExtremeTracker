function generate_reference(root)
% Read the original local source; save compact independent MATLAB fixtures.
addpath(fullfile(root,'validation','matlab'));
out=fullfile(root,'test','fixtures');
lon=(0:30:330)'; lat=(-3:3)';
inputs={}; names={};
z=false(12,7,8); z(3,3,:)=true; inputs{end+1}=z; names{end+1}='isolated_continuation';
z=false(12,7,8); z([1 12],3,:)=true; inputs{end+1}=z; names{end+1}='periodic';
z=false(12,7,8); z(3:7,3,1:4)=true; z([3 4 6 7],3,5:8)=true; inputs{end+1}=z; names{end+1}='split';
z=false(12,7,8); z([3 4 6 7],3,1:4)=true; z(3:7,3,5:8)=true; inputs{end+1}=z; names{end+1}='merge';
z=false(12,7,8); z(3:7,3,1:3)=true; z([3 4 6 7],3,4:5)=true; z(3:7,3,6:8)=true; inputs{end+1}=z; names{end+1}='split_merge';
z=false(12,7,8); z(2:4,1:2,:)=true; z(8:10,6:7,:)=true; inputs{end+1}=z; names{end+1}='two_events';
z=true(12,7,8); inputs{end+1}=z; names{end+1}='all_active';
rng(2718); for k=1:12; inputs{end+1}=rand(12,7,12)<0.23; names{end+1}=sprintf('random_%02d',k); end
outputs=cell(size(inputs)); errors=cell(size(inputs));
for k=1:numel(inputs)
 try; evalc('outputs{k}=mhwtrack(inputs{k},lon,lat,1);'); errors{k}='';
 catch e; errors{k}=e.message; end
end
save(fullfile(out,'overlap.mat'),'inputs','outputs','errors','names','lon','lat','-v7');

lon_aux=[0;1;4;90;180;359]; lat_aux=[-2;0;3]; res=2;
near=grid_nearby_index(lon_aux,lat_aux,res);
data=reshape(1:18,6,3); data(1)=NaN; data(13)=NaN; ocean=true(6,3); ocean(3,2)=false;
smooth=convzz(data,near,ocean);
save(fullfile(out,'auxiliary.mat'),'lon_aux','lat_aux','res','near','data','ocean','smooth','-v7');

rng(314); masks=cell(8,1); tracks=cell(8,1); infos=cell(8,1);
for k=1:8
 masks{k}=rand(6,5,3,7)<0.045;
 [tracks{k},infos{k}]=track_lonlat_4d_bwconncomp(masks{k},'MinDuration',1,'Verbose',false);
end
time=[8;2;2;10;5;11;6]; global_size=[6 5 3 12];
[mapped,mapped_info]=track_lonlat_4d_bwconncomp(masks{1},'MinDuration',1,'Verbose',false,'Time',time,'MaskSize',global_size,'AddIdx5Alias',true);
save(fullfile(out,'lonlat4d.mat'),'masks','tracks','infos','time','global_size','mapped','mapped_info','-v7');

seamEdges=[3 10;6 7]; masks=cell(5,1); tracks=cell(5,1); infos=cell(5,1);
for k=1:5
 masks{k}=rand(3,2,2,3,7)<0.04;
 [tracks{k},infos{k}]=track_llc_tilewise_4d_bwconncomp(masks{k},seamEdges,'MinDuration',1,'Verbose',false);
end
save(fullfile(out,'llc4d.mat'),'seamEdges','masks','tracks','infos','-v7');

% Normalization: irregular, concave shapes, split objects, seam arithmetic
% means, non-affine anomaly fields, land NaNs, and both norm flags.
lon=(10:0.25:13)'; lat=(-42:0.25:-39)';
tr=struct('day',(11:16)','xloc',{{}},'yloc',{{}});
for j=1:6
 img=false(13,13); img(4:9,4:9)=true; img(6:9,7:9)=false;
 if mod(j,2)==0; img(4,10)=true; end
 [x,y]=find(img); tr.xloc{j,1}=single(x); tr.yloc{j,1}=single(y);
end
[yy,xx]=meshgrid(lat,lon); anomaly=zeros(13,13,6);
for j=1:6; anomaly(:,:,j)=sin(xx)+cos(yy)+j^2+xx.*yy/100; end
anomaly(5,5,:)=NaN;
for norm_flag=0:1
 [radius,snapshots,normalized,xunit,yunit]=spn(tr,lon,lat,anomaly,norm_flag);
 save(fullfile(out,sprintf('spn_%d.mat',norm_flag)),'tr','lon','lat','anomaly','norm_flag','radius','snapshots','normalized','xunit','yunit','-v7');
end
fprintf('Synthetic MATLAB reference fixtures generated.\n');
end
