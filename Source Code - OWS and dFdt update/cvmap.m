function [map1,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CVmap,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedev,vari,SE]...
    = cvmap(pix,framerate,images,mask,outs,velminin,velmaxin,velalgo,MINt,MAXt,winsize,before,wint,repolopt,repolap,tfilt,usespline,splineN)

% Function for creating activation map and calculating CV using Bayly multi-vector method. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary


%% DEFINE VARIABLES FOR LATER USE
% pixel size in cm
% image = 1 beat average
% exposure ms
exposure = 1/framerate; %as input in GUI is in kHz, so 1=1000Hz which also means exposure of 1ms
pix=pix/10000; %input in micrometers, this converts it to cms


usespline=usespline-1;
splineN=1/splineN;
if usespline==0
    splineN=1;
end

%% ACTIVATION MAP
% detects point of activation and stores the time into a 2D matrix
% activation time determined by peak


[rows,cols] = size(images(:,:,1));
rawmap=zeros(rows,cols);
count = 0;
% SG filter variables
order = 3;
framesize = 11;
npeaks = rows*cols;

%% get vf* map if required
if outs == 9
    vfmap=upmap(images,mask);
end

%% peak activation map
if velalgo ==1
    for row = 1: rows
        for col = 1: cols
            count = count +1;
            signalav = sgolayfilt(double(squeeze(images(row, col,:))), order,framesize);
            signalav=imcomplement(signalav);
            time=1:length(signalav);
            fittime=min(time):splineN:max(time);
            if usespline == 1
                signalav=spline(time,signalav,fittime);
            end
            [~, peak] = max(signalav); %peak = peak indice, not value
            rawmap(row, col)=peak*splineN;
        end
    end
end
%% max upstroke, max downstroke, start and end maps
if velalgo == 2 || velalgo == 4 ||velalgo == 6 || velalgo == 7
    for row = 1: rows
        for col = 1: cols
            if mask(row,col) ~=0
                count = count +1;
                % diff signal, find upstroke
                signalav=double(squeeze(images(row, col,:)));
                if tfilt == 2
                    signalav= sgolayfilt(signalav, order,framesize);
                end
                signalav=imcomplement(signalav);
                
                rd=diff(signalav);
                [~,minind]=min(signalav);
                if minind < rd
                    [~,downstroke] = min(rd(1:minind));
                else
                    [~,downstroke] = min(rd);
                end
                time=1:length(signalav);
                fittime=min(time):splineN:max(time);
                if usespline == 1
                    signalav=spline(time,signalav,fittime);
                end
                ads=diff(signalav);
                [~,maxInd]=max(signalav);
                %         if row == 36 && col == 47
                %             figure,
                %             plot(signalav)
                %             drawnow()
                %             pause(10)
                %         end
                if length(ads) == maxInd-1
                    maxInd=maxInd -1;
                end
                ds_up=ads(1:maxInd);
                [~, upstroke] = max(ds_up);
                
                %find dpol, calc max d2F/dt2 up

                dpol = signalav(1:upstroke);
                ds=smooth(diff(dpol));
                d2s=diff(ds);
                [~,sdstart] = max(d2s);
                %
                %         if row == 32 && col == 7
                %             figure,
                %             hold on
                %             plot(signalav-min(signalav),'k')
                %             plot(rd-min(rd),'r')
                %             plot(downstroke,signalav(downstroke)-min(signalav),'ok')
                %             %plot(d2s-min(d2s),'b')
                %             figure,
                %         end
                %find repol, calc downstoke and max d2F/dt2 end
                repol=zeros(1,numel(signalav));
                for j =1:maxInd
                    repol(j) = signalav(maxInd);
                end
                for j =maxInd:numel(signalav)
                    repol(j) = signalav(j);
                end
                rds=smooth(diff(repol));
                %[~,downstroke] = min(rds);
                rd2s=diff(rds);
                if downstroke+round(25/exposure) < numel(rd2s)
                    [~,rsdstart] = max(rd2s(1:downstroke+round(25/exposure))); %the 25 exposure term is to avoid pickind up depol of next point if window too big
                else
                    [~,rsdstart] = max(rd2s);
                end
                
                %rd=diff(signalav);
                %[~,downstroke] = min(rd);
                
                
                if isempty(upstroke) == 1
                    upstroke = 0;
                end
                if isempty(sdstart) == 1
                    sdstart = 0;
                end
                if isempty(downstroke) == 1
                    disp('hi')
                    downstroke = 0;
                end
                if isempty(rsdstart) == 1
                    rsdstart = 0;
                end
                
                if velalgo == 2
                    npeaks(count) = upstroke;
                    rawmap(row, col)=upstroke;
                end
                
                if velalgo == 4
                    npeaks(count) = sdstart;
                    rawmap(row, col)=sdstart;
                end
                
                if velalgo == 6
                    npeaks(count) = downstroke;
                    rawmap(row, col)=downstroke;
                end
                
                if velalgo == 7
                    npeaks(count) = rsdstart;
                    rawmap(row, col)=rsdstart;
                end
                
                rawmap(row, col)=rawmap(row,col)*splineN;
            else
                rawmap(row,col)=0;
            end
        end
    end
end
%% repol point

%22/11/17 - Changed to be general repol map for any APD
if velalgo == 5
    checkimage=imcomplement(images);
    for row = 1: rows
        for col = 1: cols
            highVal =[];
            lowVal=[];
            count = count +1;
            if tfilt == 2
                checkimage(row, col,:) = sgolayfilt(double(squeeze(checkimage(row, col,:))), order,framesize);
            end
            dsigav = diff(checkimage(row, col,:));
            
            
            [maxi, ~] = max(checkimage(row, col,:));
            [mini, ~] = min(checkimage(row, col,:));
            amp=maxi-mini;
            [pks,locs]=findpeaks(squeeze(checkimage(row,col,:)),'MINPEAKHEIGHT',(amp/2)+mini);
            if isempty(pks) == 0
                maxi=pks(1);
                maxInd=locs(1);
            else
                [maxi, maxInd] = max(checkimage(row, col,:));
            end
            
            
            [mini, minInd] = min(checkimage(row, col, maxInd:end));
            
            minInd=minInd+maxInd; % want minimum after peak
            [~, upstroke] = max(dsigav(1:maxInd-2));
            midi=(maxi-mini)*(1-(repolap/100));
            midi=midi+mini;

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
                s_upstroke=[];
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
                %images(row, col,:) = sgolayfilt(double(squeeze(images(row, col,:))), order,framesize);
                dsigav = diff(imcomplement(images(row, col,:)));
                
                %find depol point
                
                s_dsigav = smooth(diff((imcomplement(images(row, col,:)))));
                dsigav_up=s_dsigav;
                [~, s_upstroke] = max(dsigav_up);
                %[s_maxval, s_maxInd] = max(imcomplement(images(row, col,:)));
                %[~, s_upstroke] = max(s_dsigav(1:s_maxInd-2));
                
                %[maxi, maxInd] = max(imcomplement(images(row, col,1:(s_upstroke+(round(20/exposure))))));
                [maxi, maxInd] = max(imcomplement(images(row, col,:)));
                if isempty(s_upstroke) == 0
                        s_dpol = imcomplement(squeeze(images(row,col,1:s_upstroke)));
                    
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
                    for i = 1:length(images(1,1,:))-1
                        if imcomplement(images(row, col,i)) < midi && imcomplement(images(row, col,i+1)) > midi && count50 == 0
                            lowVal= imcomplement(images(row, col,i));
                            highVal=imcomplement(images(row, col,i+1));
                            count50=count50+1;
                        end
                    end

                    % Determines points for line equations
                    if isempty(highVal) == 0 && isempty(lowVal) == 0
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
                            disp('mess up')
                            Ti=NaN;
                        end
                        npeaks(count) = Ti;
                        rawmap(row, col)= Ti;
                    end
                else
                    npeaks(count) = 0;
                    rawmap(row, col)=0;
                end
            else
                rawmap(row,col)=0;
            end
        end
    end
end


%% CALCULATE THE OFFSET TO THE FIRST PART OF THE ACTIVATION WAVE IS 1ms
% sort the unique values and find difference between them - this will find
% the timeframe where the activations begins as the difference should be
% close to the exposure time

A = unique(rawmap);
difference = diff(A);
dif = diff(difference); % Second difference
% [~, minInd] = min(floor(difference)); % finds index of smallest difference
[~, minInd] = min(abs(floor(dif))); % finds index of smallest difference


% figure('name', 'time distribution'), plot(npeaks, '.'), title('Activation time distribution');
% xlabel('pixel');
% ylabel('activation time (ms)');

A = unique(rawmap);
A
difference = diff(A)
dif = diff(difference) % Second difference
if rows == 1 || cols == 1
    train = diff([false, round(dif)==0, false]); %work around for line stacks
else
    pretrain=[false; round(dif)==0; false]
    train = diff(pretrain)% train of sequencial values
end
p = find(train==1);
q = find(train==-1);
[maxlen, ind] = max(q-p)
first = p(ind)
last = q(ind)-1
if A(first) == 0
    disp('HI!')
    A(first)=1
end
offset = A(first)-1
if isempty(offset) == 1
    %errordlg('Unable to compute conduction velocity. Possibly due to difference in activation time being too short at this framerate');
end
%map = rawmap-offset;
map=rawmap;

if isempty(offset) == 0
    map=map-offset;
end
% map = map.*double(mask);
% miniT=min(min(map(map~=0)));
% map=(map-miniT)+1;
size(map)
size(mask)
map = map.*double(mask);
% if min(min(map(mask~=0))) <= 0
%     %min(map(mask~=0))
%     %map=map-min(map(mask~=0));
% end

%% median filter that shit
isomap  = medfilt2(map, 'symmetric');
map1=isomap.*(exposure);
isomap=isomap.*exposure;

%% find edge point to get rid of quivers later
se=strel('square',winsize)
centmap=imerode(mask,se);
%% DETERMINING THE ACTIVATION POINTS
k=0;
[M,N]=size(isomap);
edgemap=zeros(M,N);
isomap(isomap<=0)=NaN;
x=zeros(1,M*N);
y=zeros(1,M*N);
t=zeros(1,M*N);
for m=1:M,
    for n=1:N,
        k=k+1;
        x(k)=n;
        y(k)=m;
        t(k)=isomap(m,n);
        if centmap(m,n) == 0 && mask(m,n) == 1
            t(k)=NaN; %get rid of edge points
        end
    end;
end;
n0=find(~isnan(t));
t=t(n0);
x=x(n0);
y=y(n0);

act_t=t;
act_x=x;
act_y=y;

%% figure,  again for refrences
%plot3(x,y,t,'.k')
%title('activation points');
%zlabel('time (ms)', 'FontSize', 20);
%xlabel('x', 'FontSize', 20);
%ylabel('y', 'FontSize', 20);

%% ESTIMATE VELOCITY VECTORS


XYT=[x',y',t'];
if velmaxin > 1
    [quivers_X,quivers_Y,~,quivers_vx,quivers_vy,v,fitxyt]=master_vel(XYT,MINt,MAXt,velminin,velmaxin,winsize,wint);
else
    [quivers_X,quivers_Y,~,quivers_vx,quivers_vy,v,fitxyt]=master_vel(XYT,MINt,MAXt,0,100,winsize,wint);
end
factor = pix*1000;
v=v.*factor;
%% outlier removal
vout=v;
quivers_vxout=quivers_vx;
quivers_vyout=quivers_vy;
quivers_Xout=quivers_X;
quivers_Yout=quivers_Y;
if outs == 2
    for i=1:length(v)
        if v(i) > velmaxin || v(i) < velminin
            vout(i) = NaN;
            quivers_vxout(i)=NaN;
            quivers_vyout(i)=NaN;
            quivers_Xout(i)=NaN;
            quivers_Yout(i)=NaN;
        end
    end
end
if outs == 9
    for i=1:length(quivers_X)
        if vfmap(quivers_Y(i),quivers_X(i)) < velminin || vfmap(quivers_Y(i),quivers_X(i)) > velmaxin
            vout(i) = NaN;
            quivers_vxout(i)=NaN;
            quivers_vyout(i)=NaN;
            quivers_Xout(i)=NaN;
            quivers_Yout(i)=NaN;
        end
    end
end

if outs == 3 || outs == 4 || outs == 5 || outs == 6 || outs == 7 || outs == 9
    if outs == 9
        onedev=std(vout);
        velmax=mean(vout)+(outs-2)*onedev;
        velmin=mean(vout)-(outs-2)*onedev;
        for i=1:length(vout)
            if vout(i) > velmax || vout(i) < velmin
                vout(i) = NaN;
                quivers_vxout(i)=NaN;
                quivers_vyout(i)=NaN;
                quivers_Xout(i)=NaN;
                quivers_Yout(i)=NaN;
            end
        end
    end
    onedev=std(v);
    velmax=mean(v)+(outs-2)*onedev;
    velmin=mean(v)-(outs-2)*onedev;
    for i=1:length(v)
        if v(i) > velmax || v(i) < velmin
            vout(i) = NaN;
            quivers_vxout(i)=NaN;
            quivers_vyout(i)=NaN;
            quivers_Xout(i)=NaN;
            quivers_Yout(i)=NaN;
        end
    end
end

if outs == 8
    if isempty(v) == 0
        vout=deleteoutliers(v);
    else
        vout=NaN
    end
end
vout=vout(~isnan(vout));
onedev=std(vout);
SE=onedev/sqrt(numel(vout));
vari=var(vout);
quivers_vxout=quivers_vxout(~isnan(quivers_vxout));
quivers_vyout=quivers_vyout(~isnan(quivers_vyout));
quivers_Xout=quivers_Xout(~isnan(quivers_Xout));
quivers_Yout=quivers_Yout(~isnan(quivers_Yout));
CV = mean(v);
CVmap=zeros(rows,cols);
CVXmap=zeros(rows,cols);
% size(CVXmap)
% quivers_Xout
% quivers_Yout
CVXmap(sub2ind(size(CVXmap),quivers_Yout,quivers_Xout)) = quivers_vxout;
sCVXmap=CVXmap.*CVXmap;
CVYmap=zeros(rows,cols);
CVYmap(sub2ind(size(CVYmap),quivers_Yout,quivers_Xout)) = quivers_vyout;
sCVYmap=CVYmap.*CVYmap;
%construct cv map
CVmap=sqrt(sCVXmap+sCVYmap);
CVmap=CVmap.*factor;
