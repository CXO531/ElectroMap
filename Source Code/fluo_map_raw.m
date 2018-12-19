% function for map showing signal levels across tissue
function [map]=fluo_map_raw(imagestack)

[rows cols] = size(imagestack(:,:,1))
 


premap = zeros(rows,cols);


for row = 1:rows
    for col = 1:cols
        signal=imcomplement(squeeze(imagestack(row,col,:))); 
        signal = (double(signal));      
              signal=signal-min(signal);
              premap(row,col)=max(signal);
    end
end

%normalise 0.02 to 1.02 %0.02 to make sure doesn't come up as white 
maxsig=max(max(premap));
minsig=min(min(premap));
map=premap-minsig;
map=map/(maxsig-minsig);
            