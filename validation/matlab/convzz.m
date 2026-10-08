function data_re=convzz(data_used,idx_full,ocean_idx)
data_re=NaN(size(data_used));
[x,y]=find(ocean_idx);
for i=1:length(x)
    idx_here=idx_full{x(i),y(i)};
    data_re(x(i),y(i))=nanmean(data_used(idx_here));
end
end