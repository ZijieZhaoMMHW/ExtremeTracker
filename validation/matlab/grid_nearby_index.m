function opt=grid_nearby_index(lon_used,lat_used,res)

% lon_used=1.25:2.5:360;
% lat_used=-88.75:2.5:90;

grid_nearby=cell(length(lon_used),length(lat_used));
% res=2.5;

for x=1:size(grid_nearby,1)
    for y=1:size(grid_nearby,2);
        x_range=[lon_used(x)-res lon_used(x)+res];
        if x_range(1)>=0 && x_range(2)<360
            logicx=lon_used>=x_range(1) & ...
                lon_used<=x_range(2);
            %x_here=find(logicx);
        elseif x_range(1)<0 && x_range(2)<360
            logicx=(lon_used>=0 & lon_used<=x_range(2)) |...
                (lon_used>=(360+x_range(1)) & lon_used<360);
            %x_here=find(logicx);
        elseif x_range(1)>=0 && x_range(2)>=360
            logicx=(lon_used>=x_range(1) & lon_used<360) |...
                (lon_used>=0 & lon_used<=(x_range(2)-360));
        end
        x_here=find(logicx);
        y_here=find(lat_used>=(lat_used(y)-res) & ...
            lat_used<=(lat_used(y)+res));
        
        [x2,y2]=meshgrid(x_here,y_here);
        x2=x2(:);
        y2=y2(:);
        
        x2(x2>length(lon_used))=x2(x2>length(lon_used))-length(lon_used);
        x2(x2<1)=length(lon_used)+x2(x2<1);
        
        idx_rm=y2<1 | y2>length(lat_used) ;
        y2(idx_rm)=[];
        x2(idx_rm)=[];
        
        grid_nearby{x,y}=[x2 y2];
    end
    
end

opt=cellfun(@(x)sub2ind([length(lon_used) length(lat_used)],x(:,1),x(:,2)),grid_nearby,'UniformOutput',false);
