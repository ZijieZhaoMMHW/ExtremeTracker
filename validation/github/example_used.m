load('snaps');
snaps=double(~isnan(snaps));
lon_used=1.25:2.5:360; % raw is 2.5 grid
lat_used=-88.75:2.5:90;
res=2.5; % 5 degree smoothing - radius 2.5

grids=grid_nearby_index(lon_used,lat_used,res); % find nearby grids
load('ocean_idx');
snaps(~ocean_idx)=nan;
data_re=convzz(snaps,grids,ocean_idx); % do the KNN smoothing
data_re=double(data_re>=0.5);

subplot(2,1,1);
m_proj('miller','lon',[0 360],'lat',[-66.5 75]);
snaps(isnan(snaps))=0;
m_pcolor(lon_used,lat_used,snaps');
shading flat
m_coast('patch',[0.7 0.7 0.7],'linewidth',2);
m_grid('xtick',[],'ytick',[],'linewidth',2);

ttl=title(['Raw'],'fontsize',18,'fontweight','bold');

subplot(2,1,2);
m_proj('miller','lon',[0 360],'lat',[-66.5 75]);
m_pcolor(lon_used,lat_used,data_re');
shading flat
m_coast('patch',[0.7 0.7 0.7],'linewidth',2);
m_grid('xtick',[],'ytick',[],'linewidth',2);

ttl=title(['KNN 5^{o}'],'fontsize',18,'fontweight','bold');
