function tasman_boundary_fixtures(root)
% Regression references for holes touching an exterior contour at a vertex.
f=load(fullfile(root,'data','tasman_sea','tracks_tas.mat'));
c=load(fullfile(root,'data','tasman_sea','lonlat_used.mat'));
r=load(fullfile(root,'validation','results','matlab_tasman_spn_0.mat'),'radius');
radius=linspace(0,r.radius,100)'; angle=linspace(0,2*pi,100);
x=radius*sin(angle); y=radius*cos(angle);
days=[7 164 165 170 188 218]; geometries=cell(numel(days),1);
for q=1:numel(days)
    j=days(q); xs=f.tracks.xloc{j}; ys=f.tracks.yloc{j};
    lon=c.lon_used(xs); lat=c.lat_used(ys); cx=mean(lon); cy=mean(lat);
    px=6371.*(lon*pi/180-cx*pi/180).*cos((lat*pi/180+cy*pi/180)/2);
    py=6371.*(lat*pi/180-cy*pi/180);
    points=[px(:),py(:)]; k=boundary(points(:,1),points(:,2),1);
    [in,on]=inpolygon(x,y,points(k,1),points(k,2));
    geometries{q}=struct('local_day',j,'absolute_day',f.tracks.day(j), ...
        'points',points,'boundary_indices',k,'inside',in|on);
end
save(fullfile(root,'test','fixtures','boundary_tasman.mat'),'geometries','x','y','-v7');
fprintf('Saved six original MATLAB Tasman Sea boundary regressions.\n');
end
