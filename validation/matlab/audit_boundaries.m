function audit_boundaries(root)
f=load(fullfile(root,'test','fixtures','spn_real_0.mat'));
geometries={}; critical=[]; maxalpha=[]; spectra={};
for i=1:numel(f.tr)
 tr=f.tr(i);
 for j=1:numel(tr.day)
  x=double(tr.xloc{j}); y=double(tr.yloc{j});
  lon=f.lon(x); lat=f.lat(y); cx=mean(lon); cy=mean(lat);
  px=6371*(lon*pi/180-cx*pi/180).*cos((cy*pi/180+lat*pi/180)/2);
  py=6371*(lat*pi/180-cy*pi/180);
  p=[px(:),py(:)]; k=boundary(p,1); s=alphaShape(p,Inf);
  geometries{i,j}=struct('points',p,'indices',k,'day',tr.day(j));
  critical(i,j)=criticalAlpha(s,'one-region'); a=alphaSpectrum(s); maxalpha(i,j)=max(a); spectra{i,j}=a;
 end
end
save(fullfile(root,'test','fixtures','boundary_audit.mat'),'geometries','critical','maxalpha','spectra','-v7');
end
