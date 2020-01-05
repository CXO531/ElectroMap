
function [map,meann,alll,onedev,vari,SE,gofmap] = tautest3(tstar,tend,framerate,maskedimage,imagestack,avbeat,outs,cmin,cmax,tfilt,rcut)
% Function for calculating relaxation constant 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary

[rows cols] = size(imagestack(:,:,1))
counter = 0;

            if tfilt == 3
            d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
            end
tic;
baseline=zeros(size(maskedimage));
premap = zeros(rows,cols);
gofmap = zeros(rows,cols);
exposure = 1/framerate; %in milliseconds
time=[1:1:length((avbeat(1,1,:)))];
time=time.*exposure;

%Check Inputted fit paramaters make sense
if tstar > tend || tend > 100 || tstar > 100 || tend < 0  || tstar < 0 || rcut > 1 || rcut < 0
    errordlg('Please check inputted tau start and end values (both should be between 0 and 100) and r2 cut off (0 and 1)')
end


count = 0;
order =3;
framesize =11;
avbeat=imcomplement(avbeat);
col_maskedimage=maskedimage(:);
col_avbeat=reshape(avbeat,rows*cols,length(avbeat(1,1,:)));
opts = optimset('Display','off');
p=[];
expmod=@(p,t) p(1)*exp(-p(2)*(t));
check = 0;
for row = 1:rows*cols
        if col_maskedimage(row) ~= 0
            ei=[];
        try
            % APD
            signalav = (squeeze((col_avbeat(row,:))));
            signalav = (double(signalav));
            signalav=signalav-min(signalav);
            % can normalise here but taken out, explained below
            if tfilt == 2
            signalav = sgolayfilt(signalav, order,framesize);
            end
            
            if tfilt == 3
            signalav=filtfilt(d,signalav);
            end

            % Calculates the baseline from first 10 points
            %BLval = signalav(1:(round(10/exposure))); % want to use 10ms not just 10 frames
            BLval=min(signalav);
            baseline(row) = mean(BLval);
            
            
            % Calculate APD
            
            [maxi, maxInd] = max(signalav(1:(end-5)));
            [lpks,llocs]=findpeaks(signalav,'MinPeakHeight',maxi/2);
            maxInd=llocs(1);
            checkSignal = signalav(maxInd:end); %checkSignal(26)=190;
            taustart=(maxi-baseline(row))*((1-(tstar/100)));
            taustart=taustart+baseline(row);
            tauend=(maxi-baseline(row))*((1-(tend/100)));
            tauend=tauend+baseline(row);
            %minInd = find(abs(checkSignal-taustart)==min(abs(checkSignal-taustart))
            minInd = find(checkSignal<taustart,1);
            
            %check new for 'beat'
            [~,newbeatInd]=min(checkSignal);

            if tend ~= 100
            min90 = find(checkSignal<tauend,1);
            if isempty(min90) == 1
                min90 = length(checkSignal)
            end
            else min90= newbeatInd;
            end
            if minInd == length(checkSignal)
                minInd = minInd-1;
            end
            if min90 == length(checkSignal) 
                min90 = min90-1;
            end
            
            if newbeatInd < minInd
                minInd = newbeatInd;
            end
              fittimeind=find(signalav==checkSignal(minInd),1);
              
              fitend = find(signalav==checkSignal(min90),1);
              
              if isempty(fitend) == 1
                  fitend=length(signalav);
              end
              if length(time(fittimeind:fitend)) <= 2
                  fitend=length(signalav);
              end
            %Find tau
            if isrow(signalav) == 1
                signalav=signalav';
            end
            
            fitsignalav=signalav(fittimeind:fitend);
            fittime=time(fittimeind:fitend)-time(fittimeind);
            if length(fittime) > 2
            if isempty(p) == 1
            p0=[fitsignalav(1),0];
            else p0=[fitsignalav(1),p(2)]; %use answer from last pixel to make close guess
            end
            count=count+1;
            p=nlinfit(fittime,fitsignalav',expmod,p0);
            tau=(1/p(2));
            %work out r2
            ei=fitsignalav'-expmod(p,fittime); %residuals
            
            %residual sum of squares
            ssres=(ei.*ei);
            ssres=sum(ssres) ;
            
            % total sum of squares
            ybar=mean(fitsignalav);
            adjy=fitsignalav-ybar;
            sstot=adjy.*adjy;
            sstot=sum(sstot);
            
            goff=1-(ssres/sstot);
            
            
            if goff < rcut %goodness of fit criteria
               tau=NaN;
               goff=NaN;
            end
            else 
                 tau = NaN;
                 goff=NaN;
            end
        catch error
            
            tau = NaN;
            goff=NaN;
            
       end
        
        premap(row) = tau;
        gofmap(row) = goff;
       end
    end
premap=reshape(premap,rows,cols)
gofmap=reshape(gofmap,rows,cols)


%% Infomation 

% disp('Some information:');
% 
% maxxxx = max(premap(premap<inf));
% minnnn = min(premap(premap>0));
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
    alll=deleteoutliers(alll);
end
onedev=std(alll);
SE=onedev/sqrt(numel(alll));
vari=var(alll);
meann = mean(alll);
% 
% disp(['APD',num2str(t)]);
% disp(['max: ', num2str(maxxxx),'ms'])
% disp(['min: ', num2str(minnnn),'ms']);
% disp(['mean :' , num2str(meannnn),'ms']);
% disp(' ');

map=medfilt2(premap);
gofmap=medfilt2(gofmap);
% get rid of -ve APD values in the map as well
for i=1:rows
    for j=1:cols
        if map(i,j) <= 0
            map(i,j) = NaN;
            gofmap(i,j) = NaN;
        end
    end
end
%gofmap
map=map*exposure;
%imshow(gofmap,[])
%colormap('jet')

