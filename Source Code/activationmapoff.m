function [actmap,aoff]=activationmapoff(pix,framerate,images,mask,velalgo,before,tfilt,usespline,splineN,repolap,inoff);
% Function for creating activation map from an image stack.
% This version also outputs the offset  
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary

exposure = 1/framerate; %as input in GUI is in kHz, so 1=1000Hz which also means exposure of 1ms
pix=pix/10000; %input in micrometers, this converts it to cms
scalelength=1;
before=round(before/exposure);
checkimage=images; %save an untouched imagestack for repol map
usespline=usespline-1;
splineN=1/splineN;
count=0;
order=3;
framesize=11;
if usespline==0
   splineN=1;
end
[rows cols]=size(images(:,:,1))
if velalgo == 1 || velalgo == 2 || velalgo == 4 || velalgo == 6
for row = 1: rows;
    for col = 1: cols;
        if mask(row,col) ~=0
        count = count +1;
        signalav=[];
        maxInd=[];
        dpol=[];
        % diff signal, find upstroke
        signalav=double(squeeze(images(row, col,:)));
        if tfilt == 2
        signalav= sgolayfilt(signalav, order,framesize);
        end
        signalav=imcomplement(signalav);
        time=1:length(signalav);
        fittime=min(time):splineN:max(time);
        if usespline == 1
        signalav=spline(time,signalav,fittime);
        end
        ads=diff(signalav);
        [~,maxInd]=max(signalav);
        if length(ads) == maxInd-1
            maxInd=maxInd -1
        end
        ds_up=ads(1:maxInd);
        [~, upstroke] = max(ads);
        ds=diff(signalav);
        [~,downstroke] = min(ds)
        
        %find dpol, calc max d2F/dt2 up
        for i =1:upstroke
            dpol(i) = signalav(i);
        end
        ds=smooth(diff(dpol));
        d2s=diff(ds);
        [~,sdstart] = max(d2s);
       
       if velalgo == 1
       rawmap(row, col)=maxInd;
       end 
        
       if velalgo == 2
       rawmap(row, col)=upstroke;
       end
       
       if velalgo == 4
       rawmap(row, col)=sdstart;
       end
       
       if velalgo == 6
       rawmap(row, col)=downstroke;
       end
       
        else
            rawmap(row,col)=0;
        end
    end
end
end

%% repolmap
if velalgo == 5
highVal =[];
lowVal=[];
checkimage=imcomplement(images);
for row = 1: rows;
    for col = 1: cols;
        highVal =[];
        lowVal=[];
        count = count +1;
        if tfilt == 2
        checkimage(row, col,:) = sgolayfilt(double(squeeze(checkimage(row, col,:))), order,framesize);
        end
        dsigav = diff(checkimage(row, col,:));
        

        [maxi, maxInd] = max(checkimage(row, col,:));
        [mini, minInd] = min(checkimage(row, col,:));
        amp=maxi-mini;
        [pks,locs]=findpeaks(squeeze(checkimage(row,col,:)),'MINPEAKHEIGHT',(amp/2)+mini);
        maxInd=[];
        maxi=[];
        if isempty(pks) == 0
        maxi=pks(1);
        maxInd=locs(1);
        else
        [maxi, maxInd] = max(checkimage(row, col,:));    
        end
        
        mini=[];minInd=[];

        [mini, minInd] = min(checkimage(row, col, maxInd:end));
        
        minInd=minInd+maxInd; % want minimum after peak
        [maxup, upstroke] = max(dsigav(1:maxInd-2));
        midi=(maxi-mini)*(1-(repolap/100));
        midi=midi+mini;
        if isempty(upstroke) == 1
            upstroke = 1;
        end
        
        if numel(pks) > 1
           stopp=locs(2);
        else
           stopp=size(images,3)-1;
        end
        
         countr=0;
        for i = maxInd:size(images,3)-1
           
            if checkimage(row, col,i) > midi && checkimage(row, col,i+1) < midi
                countr=countr+1;
                if countr==1
                highVal= checkimage(row, col,i);
                lowVal=checkimage(row, col,i+1);
                end
            end
        end


            % Determines points for line equations
            if isempty(highVal) == 0 && isempty(lowVal) == 0
            midi;
            y1 = highVal;
            y2 = lowVal;


            x1 = find(checkimage(row, col,maxInd:end)==highVal);
            x2 = find(checkimage(row, col,maxInd:end)==lowVal);
            x1=x1+maxInd;
            x2=x2+maxInd;

            if numel(x1) > 1
            x1=x1(numel(1));
            end
            if numel(x2)>1
            x2=x2(1);
            end
            m = (y2-y1)/(x2-x1);
            % Line constant, should be same for both c1 and c2
            c1 = y1-(m.*x1);
            c2 = y2-(m.*x2);

            % Time 
            Ti = (midi-c1)/m;
        npeaks(count) = Ti;
        rawmap(row, col)=Ti;
        
%         if row == 56 && col == 31
%             figure,
%             plot(squeeze(checkimage(row, col,:)),'b')
%             hold on
%             plot(locs,pks,'go')
%             plot(maxInd,(checkimage(row,col,maxInd)),'xb')
%                         line([0 120],[midi midi])
%             Ti
%             plot(round(Ti),(checkimage(row,col,round(Ti))),'or')
% 
%             figure,
%         end
%         
%          if row == 31 && col == 31
%             figure,
%             plot(squeeze(checkimage(row, col,:)),'r')
%             hold on
%             plot(locs,pks,'go')
%             plot(maxInd,(checkimage(row,col,maxInd)),'xb')
%                         line([0 120],[midi midi])
%             Ti
%             plot(round(Ti),(checkimage(row,col,round(Ti))),'or')
% 
%             figure,
%         end
        
            else 
                     
        npeaks(count) = 0;
        rawmap(row, col)=0;
            end
  end
    end
end


%% upstroke midpoint - 2/5/17 - CLEAN THIS UP A BIT NOW CHNAGED TO BEFORE ONLY (LIKE APD) 
if velalgo == 3
highVal =[];
lowVal=[];
for row = 1: rows;
    for col = 1: cols;
        if mask(row,col) ~=0
        count = count +1;
        signalav=[];
        maxInd=[];
        dpol=[];
        % diff signal, find maxi, mini and dpol start
%% come back to this         
% %         signalav= sgolayfilt(double(squeeze(images(row, col,:))), order,framesize);
% %         signalav=imcomplement(signalav);
% %         ads=diff(signalav);
% %         [maxi,maxInd]=max(signalav);
% %          ds_up=ads(1:before+round(20/exposure));            
% %          [~, upstroke] = max(ds_up);
% %         
% %         for i =1:upstroke
% %             dpol(i) = signalav(i);
% %         end
% %         
% %         ds=smooth(diff(dpol));
% %         d2s=diff(ds);
% %         [mini,sdstart] = max(d2s);
% %         
% %         % find midpoint
% %         Ti=maxInd-sdstart;
%%        
        
        
        highVal =[];
        lowVal=[];
        count = count +1;
        images(row, col,:) = sgolayfilt(double(squeeze(images(row, col,:))), order,framesize);
        dsigav = diff(imcomplement(images(row, col,:)));
%find depol point        

s_dsigav = smooth(diff((imcomplement(images(row, col,:)))));
dsigav_up=s_dsigav(1:before+round(40/exposure));
[~, s_upstroke] = max(dsigav_up); 
%[s_maxval, s_maxInd] = max(imcomplement(images(row, col,:)));
%[~, s_upstroke] = max(s_dsigav(1:s_maxInd-2));

%[maxi, maxInd] = max(imcomplement(images(row, col,1:(s_upstroke+(round(20/exposure))))));
[maxi, maxInd] = max(imcomplement(images(row, col,:)));
if isempty(s_upstroke) == 0
for i =1:s_upstroke
s_dpol(i) = imcomplement(images(row,col,i));
end

% sD = find(diff(s_dpol)>0);
% st = diff([0,round(diff(diff(sD)))==0,0]);
% sp = find(st==1);
% sq = find(st==-1);
% [smaxlen, sind] = max(sq-sp);
% sfirst = sp(sind);
% sdstart = sD(sfirst);% depolarisation start point
% 
% % ds=smooth(diff(s_dpol));
% %             d2s=diff(ds);
% %             [~,sdstart] = max(d2s);
% 
% 
% mini=imcomplement(images(row, col, sdstart));  
%if isempty(mini) == 1
    mini=min(imcomplement(images(row, col, 1:maxInd)));  
%end
        
        [maxup, upstroke] = max(dsigav(1:maxInd-2));
        midi=(maxi-mini)*0.5;
        midi=midi+mini;
        if isempty(upstroke) == 1
            upstroke = 0;
        end
        count50= 0; %switch to find first time 50% amplitude reached (i.e. in upstroke, not repol)
        for i = 1:before+round(50/exposure)
            if imcomplement(images(row, col,i)) < midi && imcomplement(images(row, col,i+1)) > midi && count50 == 0
                lowVal= imcomplement(images(row, col,i));
                highVal=imcomplement(images(row, col,i+1));
                count50=count50+1;
            end
        end
        count50=0;
            % Determines points for line equations
            if isempty(highVal) == 0 && isempty(lowVal) == 0
            midi;
            y1 = highVal;
            y2 = lowVal;


            x1 = find(imcomplement(images(row, col,:))==highVal);
            x2 = find(imcomplement(images(row, col,:))==lowVal);
            x1=x1+1;
            x2=x2+1;
            
            if numel(x2) > 1
            x2=x2(numel(x2));
            end
            if numel(x1)>1  
            x1=x1(1);
            end
            m = (y2-y1)/(x2-x1);
            % Line constant, should be same for both c1 and c2
            c1 = y1-(m.*x1);
            c2 = y2-(m.*x2);
           
            % Time 
            Ti = (midi-c1)/m;
            if isempty(Ti) == 1
                disp('fuck up')
                Ti=NaN;
            end
        npeaks(count) = Ti;
        rawmap(row, col)= Ti;
        
            else 
        npeaks(count) = 0;
        rawmap(row, col)=0;
            end
end

if isempty(s_upstroke) ~= 0
        npeaks(count) = 0;
        rawmap(row, col)= 0;
end
        else 
            rawmap(row,col) = 0;
        end
    end
    end
end

A = unique(rawmap);
A;
difference = diff(A);
dif = diff(difference); % Second difference
if rows == 1 || cols == 1
    train = diff([false, round(dif)==0, false]); %work around for line stacks 
else
pretrain=[false; round(dif)==0; false];
train = diff(pretrain);% train of sequencial values
end
p = find(train==1);
q = find(train==-1);
[maxlen, ind] = max(q-p);
first = p(ind);
last = q(ind)-1;
if A(first) == 0
    disp('HI!')
    A(first)=1;
end
offset = A(first)-1;
if isempty(offset) == 1
    %errordlg('Unable to compute conduction velocity. Possibly due to difference in activation time being too short at this framerate');
end
%map = rawmap-offset;
map=rawmap;
offset;
if isempty(offset) == 0;
    inoff
    if isempty(inoff) == 1
map=map-offset;
aoff=offset;
    else
        map=map-inoff;
    end
end

% map = map.*double(mask);
% miniT=min(min(map(map~=0)));
% map=(map-miniT)+1;
map = map.*double(mask);
isomap  = medfilt2(map, 'symmetric');
actmap=isomap.*(exposure);
       