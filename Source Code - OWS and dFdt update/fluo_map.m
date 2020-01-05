
function [map]=fluo_map(framerate,maskedimage,imagestack,tfilt,avbeat)
% function for map showing signal levels across tissue
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary

[rows cols] = size(avbeat(:,:,1))
 
if tfilt == 3 % filter for later if needed
            d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
end
 
order =3; %SG filter variables
framesize =11;

premap = zeros(rows,cols);

for row = 1:rows
    for col = 1:cols
        if maskedimage(row,col) ~= 0
        signal=imcomplement(squeeze(avbeat(row,col,:))); 
        signal = (double(signal));
      
            %Temporal Filter
            if tfilt == 2
            signal = sgolayfilt(signal, order,framesize);
            end
            
            if tfilt == 3
            signal=filtfilt(d,signal);
            end
            
              %set to zero
              signal=signal-min(signal);
%               if row == 30 && col == 18
%                   figure,
%                   plot(signal,'r')
%               end
%               if row == 16 && col == 15
%                   figure,
%                   plot(signal,'b')
%               end
              premap(row,col)=max(signal);
        else premap(row,col) = NaN;
        end
    end
end
%figure,
%normalise 0.02 to 1.02 %0.02 to make sure doesn't come up as white 
maxsig=max(max(premap))
minsig=min(min(premap))
map=premap
% map=premap-minsig;
% map=map/(maxsig-minsig);
% map=map+0.02;

            