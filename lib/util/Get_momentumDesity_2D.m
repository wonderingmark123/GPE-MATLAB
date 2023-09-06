function D=Get_momentumDesity_2D(DensityCell,x,y)

if  isgpuarray( DensityCell)
    DensityCell = gather(DensityCell);

end

len = length(DensityCell);
D = [];
for i = 1:len
    DensityNow = DensityCell{i};
    DensityNow = DensityNow(x,y);
    DensityNow = reshape(DensityNow,1,[]);
    D = [D;DensityNow];
end
end