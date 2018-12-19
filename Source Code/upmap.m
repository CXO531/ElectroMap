function [map1] = upmap(images,mask)
% function for producing map of time at which maximum upstroke occurs
% (Based on - Zemlin et al 2008)
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary

order=3; framesize=11;      
[rows cols] = size(images(:,:,1));
rawmap=zeros(rows,cols);
usespline=0;
for row = 1: rows;
    for col = 1: cols;
        if mask(row,col) ~=0
            try
        signalav=[];
        maxInd=[];
        dpol=[];
        % diff signal, find upstroke
        %signalav= sgolayfilt(double(squeeze(images(row, col,:))), order,framesize);
        signalav= squeeze(images(row, col,:));
        signalav=imcomplement(signalav);
        time=1:length(signalav);
        fittime=min(time):0.2:max(time);
        if row == 32 && col == 32
        end
        %normalise signalav
        signalav=signalav-min(signalav);
        signalav=signalav/max(signalav);
        if usespline==1
        ss=spline(time,signalav,fittime);
        else 
        ss=signalav
        end
            
        %ads=diff(signalav);
        ads=diff(ss);
        [~,maxInd]=max(signalav);
        if row == 30 && col == 30
%            figure,
%            subplot(1,3,1)
%            hold on
%            plot(time,signalav,'kx')
%            plot(fittime,ss,'b-');
%            pss=pchip(time,signalav,fittime);
%            plot(fittime,pss,'r-');
%            dno=diff(signalav);
%            dno=dno-min(dno);
%            dno=dno./max(dno);
%            dyes=diff(ss);
%            dyes=dyes-min(dyes);
%            dyes=dyes./max(dyes);
%            dpyes=diff(pss);
%            dpyes=dpyes-min(dpyes);
%            dpyes=dpyes./max(dpyes);
%            subplot(1,3,2)
%            hold on
%            plot(time(2:end),dno,'kx')
%            plot(fittime(2:end),dyes,'b-')
%            plot(fittime(2:end),dpyes,'r-')
%            subplot(1,3,3)
%            hold on
%            figure
        end
        
%         if col == 31 && row == 6
%             figure,
%             plot(signalav)
%             title(num2str(col))
%         end
%         
%         if col == 30 && row == 37
%             figure,
%             plot(signalav)
%             title(num2str(col))
%         end
% 
%         if col == 25 && row == 61
%             figure,
%             plot(signalav)
%             title(num2str(col))
%             figure
%         end
%         
        if length(ads) == maxInd-1
            maxInd=maxInd -1;
        end
        ds_up=ads(1:maxInd);
        [~, upstroke] = max(ads);
        
        vf=ss(upstroke);
        rawmap(row,col)=vf;
        catch
            rawmap(row,col)=0;
        end
        
  end
end
end



map=rawmap;
map=medfilt2(map);
map1 = map.*double(mask);
