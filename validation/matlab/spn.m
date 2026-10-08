%% radius
function [radius_s,sst_sp,sst_spn_full,x_full,y_full]=spn(tracks,lon_used,lat_used,data_anom,norm_flag)

tracks_new=squeeze(struct2cell(tracks));

radius_sp=cell(length(tracks),1);
[lat2,lon2]=meshgrid(lat_used,lon_used);

for i=1:length(tracks)
    tic
    t_here=tracks_new{1,i};
    radius_s=cell(length(t_here),1);
    xloc=tracks_new{2,i};
    yloc=tracks_new{3,i};
    for j=1:length(t_here)
        %         mhw_here=sst_anom(:,:,t_here(j));
        ind_here=sub2ind([length(lon_used) length(lat_used)],xloc{j},yloc{j});
        img_here=zeros(length(lon_used),length(lat_used));
        img_here(ind_here)=1;
        bw=bwconncomp(img_here,8);
        %sst_b=NaN(100,100,bw.NumObjects);
        radius_b=NaN(1,1);
        for b=1:1
            idx_here=ind_here;
            center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
            dist_here=geodist(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
            
            radius_b(b)=nanmax(dist_here);
        end
        %         sst_s(:,:,j)=nansum(sst_b.*reshape(area_b,1,1,length(area_b)),3)./nansum(area_b);
        %center_s(j,:)=nansum(center_b.*area_b,1)./nansum(area_b);
        %         x_s{j}=x_b;y_s{j}=y_b;
        %area_s{j}=area_b;
        radius_s{j}=radius_b;
    end
    radius_sp{i}=radius_s;
    toc
end
radius_s=NaN(length(radius_sp),1);
for i=1:length(radius_sp)
    radius_here=radius_sp{i};
    radius_here=cellfun(@nanmax,radius_here);
    radius_s(i)=nanmax(radius_here);
end

%% spn
tracks_new=squeeze(struct2cell(tracks));
t_used=tracks_new(1,:);
t_min=nanmin(cellfun(@nanmin,t_used));
t_max=nanmax(cellfun(@nanmax,t_used));

% mfile=matfile(['/Volumes/new_drive/lagrangian/sun/sst_anom']);
% data_anom=mfile.sst_anom(:,:,t_min:t_max);

sst_spn_full=NaN(100,100,5,length(tracks),1);

for m=1:1
    tic
    sst_anom=data_anom(:,:,:,m);
    sst_sp=cell(length(tracks),1);
    x_sp=cell(length(tracks),1);
    y_sp=cell(length(tracks),1);
    area_sp=cell(length(tracks),1);
    [lat2,lon2]=meshgrid(lat_used,lon_used);
    res=100;
    
    for i=1:length(tracks)
        t_here=tracks_new{1,i}-t_min+1;
        sst_s=NaN(res,res,length(t_here));
        x_s=cell(length(t_here),1);
        y_s=cell(length(t_here),1);
        center_s=NaN(length(t_here),2);
        area_s=cell(length(t_here),1);
        xloc=tracks_new{2,i};
        yloc=tracks_new{3,i};
        for j=1:length(t_here)
            mhw_here=sst_anom(:,:,t_here(j));
            ind_here=sub2ind([length(lon_used) length(lat_used)],xloc{j},yloc{j});
            img_here=zeros(length(lon_used),length(lat_used));
            img_here(ind_here)=1;
            %bw=bwconncomp(img_here,8);
            sst_b=NaN(res,res,1);
            x_b=NaN(res,res,1);
            y_b=NaN(res,res,1);
            %         area_b=NaN(1,1);
            center_b=NaN(1,2);
            for b=1:1
                idx_here=ind_here;
                center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
                center_b(b,:)=center_here;
                %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
                %             +(lat2(idx_here)-center_here(2)).^2);
                %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
                [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
                %             area_here=length(idx_here)*0.0625;
                %             area_b(b)=area_here;
                radius_here=radius_s(i);
                radius_range=linspace(radius_here,0,res);
                radius_range=radius_range(end:-1:1);
                x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
                for r=1:res
                    x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
                    y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
                end
                [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
                F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','none');
                sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
                k=boundary(lond(idx_here),latd(idx_here),1);
                if isempty(k)
                    k=1:length(lond(idx_here));
                end
                lon2k=lond(idx_here);lat2k=latd(idx_here);
                [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
                if norm_flag==0
                sst_b_here(~in & ~on)=0;
                end
                sst_b(:,:,b)=sst_b_here;
                x_b(:,:,b)=x_full;
                y_b(:,:,b)=y_full;
            end
            sst_s(:,:,j)=sst_b;
            center_s(j,:)=center_b;
            x_s{j}=x_b;y_s{j}=y_b;
        end
        sst_sp{i}=sst_s;
        x_sp{i}=x_s;
        y_sp{i}=y_s;
    end
    
    t_used=0:0.25:1;
    radius_range=linspace(1,0,res);
    radius_range=radius_range(end:-1:1);
    x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
    for r=1:res
        x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
        y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
    end
    
    sst_spn=NaN(res,res,5,length(sst_sp));
    ta=(0:0.25:1)';
    for i=1:length(sst_sp)
        tic
        sst_here=sst_sp{i};
        tb=(1:size(sst_here,3))./size(sst_here,3);
        sst_b=NaN(50,50,length(ta));
        for x=1:size(sst_here,1)
            for y=1:size(sst_here,2)
                ts_b=squeeze(sst_here(x,y,:));
                ts_a=interp1(tb,ts_b,ta,'linear','extrap');
                sst_b(x,y,:)=ts_a;
            end
        end
        sst_spn(:,:,:,i)=sst_b;
        toc
    end
    sst_spn_full(:,:,:,:,m)=sst_spn;
    
    toc
end
end
function [d2km,x,y]=geodist(lat1,lat2,lon1,lon2)
radius=6371;
lat1=lat1*pi/180;
lat2=lat2*pi/180;
lon1=lon1*pi/180;
lon2=lon2*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
% a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon/2).^2;
% c=2*atan2(sqrt(a),sqrt(1-a));
% d1km=radius*c;    %Haversine distance

x=radius.*deltaLon.*cos((lat1+lat2)/2);
y=radius.*deltaLat;
d2km=sqrt(x.*x + y.*y); %Pythagoran distance

end