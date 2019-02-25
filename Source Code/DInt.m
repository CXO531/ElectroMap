function [map,meann,alll] = DInt(framerate,maskedimage,imagestack,outs,cmin,cmax,tfilt,before,after,peakdist,frame_1,frame_last,t);
% Function for creating diastolic interval map from an image stack. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary


[rows cols num_images] = size(imagestack(:,:,:));
counter = 0;

 if tfilt == 3
            d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
 end
tic;
mask=maskedimage;
baseline=zeros(size(maskedimage));
DI = zeros(rows,cols);

exposure = 1/framerate; %in milliseconds
before=round(before/exposure);
after=round(after/exposure);

imagestack=imagestack(:,:,(frame_1-before):(frame_last+after));
[rows cols num_images] = size(imagestack(:,:,:));
count = 0;
order =3;
framesize =11;
figcount = 0;

for row=1:rows;
    for col=1:cols;
     try   
sdstart=[]; Tim = [];
if mask(row,col) ~= 0  
%try
%find pixel signal
    sig=imcomplement((squeeze(imagestack(row,col,:))));


sig=double(sig);
sig=sig-min(sig);

 if tfilt == 2
 sig = sgolayfilt(sig, order,framesize);
 end
            
 if tfilt == 3
 sig=filtfilt(d,sig);
 end
 peakdist;
 [pks locs]=findpeaks(sig,'MINPEAKHEIGHT', (max(sig)-min(sig))/2, 'MINPEAKDISTANCE', peakdist);
 
 if (locs(1)-before) < 1 %if first peak too close too start or last peak too close to end 
    start = 2;
 else start = 1;
 end
 
if (locs(end)+after) > num_images
finish =numel(locs) -1;
else finish = numel(locs);
end

%% make double peak average
for i = start
    tstart=locs(i)-before;
    tend=locs(i+1)+after;
    overalldoublepeak=sig(tstart:tend);
end

if  (locs(finish) +after ) <= numel(sig)
for i=start+1:finish-1
    tstart=locs(i)-before;
    tend=locs(i+1)+after;
    doublepeak=sig(tstart:tend);
    if numel(doublepeak) > numel(overalldoublepeak)
    for j=1:numel(overalldoublepeak);
    overalldoublepeak(j)=overalldoublepeak(j)+doublepeak(j);
    end
    end
     if numel(doublepeak) <= numel(overalldoublepeak)
    for j=1:numel(doublepeak);
    overalldoublepeak(j)=overalldoublepeak(j)+doublepeak(j);
    end
     end   
end
end
overalldoublepeak=overalldoublepeak-min(overalldoublepeak);
[opks olocs]=findpeaks(overalldoublepeak,'MINPEAKHEIGHT', (max(overalldoublepeak)-min(overalldoublepeak))/2, 'MINPEAKDISTANCE', peakdist);
if numel(opks) < 2
[opks olocs]=findpeaks(overalldoublepeak,'MINPEAKHEIGHT', (max(overalldoublepeak)-min(overalldoublepeak))/5, 'MINPEAKDISTANCE', peakdist);
end
peak=overalldoublepeak;
%% find apd of first peak

            maxval=opks(1);
            maxInd=olocs(1);
            baseline=min(peak);
            APD = (maxval-baseline)*(1-(t/100))+baseline;
            
 
             maxInd
             olocs(2)
             before
            checkSignal = peak(maxInd:olocs(2)-before);
            
            minInd = find(abs(checkSignal-APD)==min(abs(checkSignal-APD)));
            if minInd == numel(checkSignal)
                minInd=minInd - 1;
            end
            % Locates points above and below APD90
            if(checkSignal(minInd) > APD)
                highVal = checkSignal(minInd);
                if minInd < numel(checkSignal)
                lowVal = checkSignal(minInd+1);
                else
                highVal = checkSignal(minInd-1);   
                 lowVal = checkSignal(minInd);
                end
            else
                highVal = checkSignal(minInd-1);
                lowVal = checkSignal(minInd);
            end

            % Solution for detecting a slight positive gradient next to APD region
%             if (lowVal - highVal)>0
%                 highVal = checkSignal(minInd-2);
%                 lowVal = checkSignal(minInd-1);
% %                 disp('**There was a positive gradient found APD70**');
%             end

            % Determines points for line equations
            y1 = highVal(1);
            y2 = lowVal(1);
            
            x1 = find(peak==highVal(1));
            x2 = find(peak==lowVal(1));
            % If encountering a flat region gradient becomes inf so only takes the
            % section of the flat part closest to the APD point
            x1=x1(end);
            x2=x2(1);
            % Gradient of line y=mx+c
            m = (y2-y1)/(x2-x1); 
            % Line constant, should be same for both c1 and c2
            c1 = y1-(m.*x1);
            c2 = y2-(m.*x2);

            % Time and APD70
            Ti = (APD-c1)/m;
            Tim=Ti;
%% find dpol point of 2nd peak
%peak2=peak(olocs(2)-before:end);
% dpeak2=diff(peak2);
% [p2maxval, p2maxInd] = max(peak2);
% [~, p2upstroke] = max(dpeak2(1:p2maxInd-2));
% dpol = peak2(1:p2upstroke);
% sD = find(diff(dpol)>0);
% st = diff([0,round(diff(diff(sD)))==0,0]);
% sp = find(st==1);
% sq = find(st==-1);
% [smaxlen, sind] = max(sq-sp);
% sfirst = sp(sind);
% sdstart = sD(sfirst); % depolarisation start point
% sdstart=sdstart+(olocs(2)-before);

% use d2f/dt2 instead
peak2=peak(olocs(2)-(before-1):end);
dpol=peak2(1:before);
ds=diff(dpol);
d2s=diff(ds);
[~,sdstart] = max(d2s);
sdstart=sdstart+(olocs(2)-(before-1));


%% DI calc
if isempty(sdstart) == 0 && isempty(Tim) == 0
DI(row,col)=(sdstart-Tim);
else DI(row,col)=NaN;
end

             end
catch error
    DI(row,col)=NaN;
% if row == 20 && col == 50
% figure,
% plot(overalldoublepeak)
% hold on
% plot(olocs(2),overalldoublepeak(olocs(2)),'sr')
% plot(round(Tim),overalldoublepeak(round(Tim)),'bo')
% plot(sdstart,overalldoublepeak(sdstart),'go')
% figure,
% olocs
% opks
% end

end % if not bg end
%       catch
%          DI(row,col) = NaN;
    %end
end %col end
    end %row end

premap = DI;

%% stuff for output

alll=premap(~isnan(premap)&~isinf(premap));
alll=alll(alll>0); % get rid of -ve APD values
alll=alll*exposure;
if outs == 2
   alll=alll(alll>cmin);
   alll=alll(alll<cmax);
end

if outs == 3 || outs == 4 || outs == 5 || outs == 6 || outs == 7
   onedev=std(alll);
   maxout=mean(alll)+(outs-2)*onedev;
   minout=mean(alll)-(outs-2)*onedev;
   alll=alll(alll>minout);
   alll=alll(alll<maxout);
end

if outs == 8
    if isempty(alll) ~= 1
    alll=deleteoutliers(alll);
    end
end
onedev=std(alll);
SE=onedev/sqrt(numel(alll));
vari=var(alll);
meann = mean(alll);

map=medfilt2(premap);
for i=1:rows
    for j=1:cols
        if map(i,j) <= 0
            map(i,j) = NaN;
        end
    end
end
map=map*exposure;


