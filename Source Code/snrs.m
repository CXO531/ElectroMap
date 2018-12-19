function [mapSNRr,mapSNRdb,allSNRr,allSNRdb]=SNRs(images,mask,before,after,tfilt);
% Function for calculating SNR from amplitude and s/d of noise . 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary


 if tfilt == 3
    d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
 end
order=3; framesize=11;
for r = 1:size(images,1)
    for c = 1:size(images,2)
        if mask(r,c) ~= 0
            signalav = imcomplement(squeeze((images(r,c,:))));
            signalav = (double(signalav));
            
            if tfilt == 2
                signalav = sgolayfilt(signalav, order,framesize);
            end
            if tfilt == 3
                signalav=filtfilt(d,signalav);
            end
            signalav=signalav-min(signalav);
            [maxval, maxInd] = max(signalav);
            noiseinds=zeros(numel(signalav),1);
            
            for m=-before:after
                if (maxInd+m) > 1 && (maxInd+m) <= numel(signalav)
                    noiseinds(maxInd+m)=1;
                end
            end
            for j=1:numel(signalav)
                if noiseinds(j) == 0
                    noise(j)=signalav(j);
                else
                    noise(j)=NaN;
                end
            end
            amp1(r,c)=(maxval-nanmean(noise));
noise1(r,c)=(nanstd(noise)-nanmean(noise));
mapSNRr(r,c)=(maxval-nanmean(noise))/(nanstd(noise)); %changed because pks amp also shifted by mean noise. 
mapSNRdb(r,c)=20*log10(mapSNRr(r,c));

        else
            noise1(r,c)=NaN;
            mapSNRr(r,c)=NaN;
            mapSNRdb(r,c)=NaN;
        end
    end
end

allSNRr=mapSNRr(isnan(mapSNRr)~=1);
allSNRdb=mapSNRdb(isnan(mapSNRdb)~=1);
