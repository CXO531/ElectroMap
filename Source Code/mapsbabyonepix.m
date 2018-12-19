
function [signalav] = mapsbabyonepix(startopt,framerate,t,maskedimage,imagestack,avbeat,row,col,colopt,before,after,apdblopt,apdblnum,tstar,tend,normalise,tfilt)
% function for taking a single beat image file and extracting and displaying 
% signal from chosen pixel
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary
[rows cols num] = size(imagestack(:,:,:))
[~,~,num] = size(avbeat(:,:,:))

tic;
 if tfilt == 3
            d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
            end

premap30 = zeros(rows,cols);

exposure = 1/framerate; %in milliseconds
before=round(before/exposure);
timelength = round((before+after)/exposure)
if num >= timelength
avbeat=avbeat(:,:,1:timelength);
time=(1:timelength)*exposure;
else time =(1:num)*exposure;
end
count = 0;
order =3;
framesize =11;
tend=tend/100;
            % APD
            APD = t/100;
            signalav = imcomplement(squeeze(avbeat(row,col,:)));
            signalav = (double(signalav)); % can normalise here but taken ouT, explained below
%           plot(signalav), title(['row: ', num2str(row), ' col:', num2str(col)]);
            if tfilt == 2
            signalav = sgolayfilt(signalav, order,framesize);
            end
            if tfilt == 3
              signalav=filtfilt(d,signalav);  
            end
            signalav=signalav-min(signalav);
            if normalise == 1
                signalav=signalav./max(signalav);
            end
%             hold on
%             plot(sgolayfilt(signalav, order, framesize),'r');
%             hold off
%             pause(0.005);
            
            % Obtain Initial Criteria for APDs

        

            % Calculates the baseline from first 10 points
%             BLval = signalav(1:(round(10/exposure/1000))); % want to use 10ms not just 10 frames

            
            %baseline= min(min(signalav));
           
            
                        %Possible Refrence points
            
            dsigav = diff(signalav);
            [maxval, maxInd] = max(signalav);
            maxInd
            dsigav_up=dsigav(1:maxInd);            
            % 1, max upstroke (upstroke)
            [~, upstroke] = max(dsigav_up);
            
            % 2, Peak (maxInd)
            %[maxval, maxInd] = max(signalav(1:(upstroke+round(25/exposure)))); %peak assotiead with first upstroke shortley after maxupstroke (set here to 20ms)
            [maxval, maxInd] = max(signalav);
            % 3, Depol point (sdstart) (should be changed to max d2f/dt2???) 
            for i =1:upstroke
                dpol(i) = signalav(i);
            end
            % Uncomment for figure of derivatives
%              
%                 hold on
%                 plot(time,signalav,'kx');
%                 ds=diff(signalav);
%                 d2s=diff(ds);
%                 plot(time(1:end-1),ds,'bx-')
%                 plot(time(2:end-1),d2s,'rx-')
%                 axis tight
%                 hy = plot(xlim, [0 0], 'Color', [.8 .8 .8]);
            
            %New dpol find code (!! Do comparsion of these for thesis!)
            ds=smooth(diff(dpol));
            d2s=diff(ds);
            [~,sdstart] = max(d2s);
            
            % Old dpol point find code (Ting)   
%             sD = find(diff(dpol)>0)
%             st = diff([0,round(diff(diff(sD)))==0,0])
%             sp = find(st==1)
%             sq = find(st==-1)
%             [smaxlen, sind] = max(sq-sp)
%             sfirst = sp(sind)
%             sdstart = sD(sfirst)
            %pause(10)
            
            % 4, Depol midpont
            mini=signalav(sdstart);  
            if isempty(mini) == 1
                mini=min(signalav(1:upstroke));  
            end
            maxi=maxval;
            midi=(maxi-mini)*0.5;
            midi=midi+mini;
%             if before+round(50/exposure) > length(signalav)
%                 midicheck=length(signalav)-1
%             else midicheck = before+round(50/exposure) 
%             end
%             for i = 1:midicheck-1 
%             if signalav(i) < midi && signalav(i+1) > midi
%                 lowVal= signalav(i);
%                 highVal=signalav(i+1);
%             end
%             end

            ind1=find(signalav > midi);
            if isempty(sdstart) == 1
                sdstart=1;
            end
            ind1=ind1(ind1>sdstart);

            if isempty(ind1) == 1
                ind1 = 2
            end
            ind1=ind1(1)
            ind2=ind1-1;
            lowVal=signalav(ind2);
            highVal=signalav(ind1);

            
            
            % repol points 
            repol=zeros(1,numel(signalav))
            for j =1:maxInd
                repol(j) = signalav(maxInd);
            end
            for j =maxInd:numel(signalav)
                repol(j) = signalav(j);
            end
            rds=smooth(diff(repol));
            [~,downstroke] = min(rds);
            rd2s=diff(rds);
            if downstroke+round(25/exposure) < numel(rd2s)
            [~,rsdstart] = max(rd2s(1:downstroke+round(25/exposure))); %the 25 exposure term is to avoid pickind up depol of next point if window too big
            else
            [~,rsdstart] = max(rd2s);
            end
            % Determines points for line equations
            if isempty(highVal) == 0 && isempty(lowVal) == 0
            midi
            y1 = highVal;
            y2 = lowVal;


            x1 = find(signalav==highVal);
            x2 = find(signalav==lowVal);
            x1=x1
            x2=x2

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
            depol_mid = (midi-c1)/m;
            if isempty(depol_mid) == 1
                disp('fuck up')
                depol_mid=NaN;
            end
            end
            
            %plot(dsigav)
            % Calculate APD
            %figure, plot(signalav)
            if apdblopt == 1
            blsec=(round(apdblnum/exposure));
            BLval = signalav(1:blsec); 
            baseline = mean(BLval)
            end
            
            if apdblopt == 2
            blsec=(round(apdblnum/exposure));
            BLval = signalav((length(signalav)-blsec):length(signalav)); 
            baseline = mean(BLval);
            end
            
            if apdblopt == 3
            baseline= min(signalav(1:upstroke))
            end
            
            if apdblopt == 4
                aftsig=signalav(maxInd:length(signalav));
                baseline = min(aftsig);
            end
            
           maxval
            APD = (maxval-baseline)*(1-APD)+baseline
            
            checkSignal = signalav(maxInd:end);%checkSignal(26)=190;
            [~,min2]=min(checkSignal);
            checkSignal=checkSignal(1:min2); %ignore 2nd beat if present
            % [~, minInd] = min(abs(checkSignal - APD70));

            %minInd = find(abs(checkSignal-APD)==min(abs(checkSignal-APD)));;
            minInd = find(checkSignal<APD,1);
            % Locates points above and below APD70
            if(checkSignal(minInd) > APD)
                highVal = checkSignal(minInd);
                lowVal = checkSignal(minInd+1);
            else
                highVal = checkSignal(minInd-1);
                lowVal = checkSignal(minInd);
            end

            % Solution for detecting a slight positive gradient next to APD region
            if (lowVal - highVal)>0
                highVal = checkSignal(minInd-2);
                lowVal = checkSignal(minInd-1);
%                 disp('**There was a positive gradient found APD70**');
            end

            % Determines points for line equations
            y1 = highVal;
            y2 = lowVal;

%             x1 = time(find(signalav==highVal)); % 'time' no value here
%             x2 = time(find(signalav==lowVal));
            if isempty(highVal) == 1 && isempty(lowVal) == 0
                highVal=lowVal;
            end
            highVal
            lowVal
            if isempty(highVal) == 1
                highVal = 10
                lowVal = 9
            end
            x1 = find(signalav==highVal);
            x2 = find(signalav==lowVal);
            % If encountering a flat region gradient becomes inf so only takes the
            % section of the flat part closest to the APD point
            try
            x1=x1(end);
            x2=x2(1);
            % Gradient of line y=mx+c
            m = (y2-y1)/(x2-x1); 
            % Line constant, should be same for both c1 and c2
            c1 = y1-(m.*x1);
            c2 = y2-(m.*x2);

            % Time and APD
            Ti = (APD-c1)/m;
            
            if startopt == 1
            apd = Ti-(upstroke);
            elseif startopt == 2
            apd = Ti-(maxInd);
            elseif startopt == 3
            apd = Ti-(sdstart);
            elseif startopt == 4
            apd = Ti-(depol_mid);
            end
            catch
                Ti=0;
                apd=0;
            end
            Ti=Ti*exposure;
            apd=apd*exposure;
if colopt == 1;
hold on
if length(time) > length(signalav)
    time=time(1:length(signalav));
elseif length(signalav) > length(time)
    signalav=signalav(1:length(time));
end

%% Find me man tau bruv
exposure = 1/framerate; %in milliseconds
time=[1:1:length((avbeat(1,1,:)))];
time=time.*exposure;
count = 0;
order =3;
framesize =11;


signalav = imcomplement(squeeze((avbeat(row,col,:))));
            signalav = (double(signalav));
            
           %signalav = sgolayfilt(signalav, order,framesize);
           signalav=signalav-min(signalav);
            if normalise == 1
                signalav=signalav./max(signalav);
            end
            % Calculate APD
            
            [maxi, maxInd] = max(signalav);
            checkSignal = signalav(maxInd:end); %checkSignal(26)=190;
            taustart=(maxi)*((1-(tstar/100)));
            taustart=taustart+baseline;
            minInd = find(abs(checkSignal-taustart)==min(abs(checkSignal-taustart)));
            length(checkSignal);
            if minInd == length(checkSignal)
                minInd = minInd-1;
            end
            checkSignal
            if minInd == 0
                minInd = 1
            end
            if(checkSignal(minInd) < taustart)  %only care about lowval for tau
                minInd-1
                lowVal = checkSignal(minInd);
            else
                lowVal = checkSignal(minInd);
            end
            minInd=minInd(1);
            fittimeind = find(signalav==lowVal(1));
            
            % Find apd90 and stop exp fit from this point onwards
            
            APD90 = (maxval-baseline)*(1-tend)+baseline;
            min90 = find(abs(checkSignal-APD90)==min(abs(checkSignal-APD90)));
            min90 = min90(1);
            if min90 == length(checkSignal)
                min90 = minInd-1;
            end
           if min90 == 0
               min90=1
           end
           min90
           APD90
            if(checkSignal(min90) > APD90) && min90 > 1
                lowVal90 = checkSignal(min90+1);
            else
                lowVal90 = checkSignal(min90);
              end
            
           fitend = find(signalav==lowVal90);
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
            if length(time(fittimeind:fitend)) > 2
            [f,gof] = fit(time(fittimeind:fitend)',signalav(fittimeind:fitend),'exp1');
            fitdecaybi=@(a,b,c,d,x) a*exp(b*x)+c*exp(d*x);
            fitdecaymono=@(e,g,x) e*exp(g*x);
            go=gof.rsquare;
            tau_early=(-1/f.b);
            tau_early=tau_early.*exposure;
            %tau_late=(-1/f.d);
            end

          
plot(time, signalav,'b', 'LineWidth', 1);
plot(time, signalav,'kx', 'LineWidth', 2);
plot(time(fittimeind:fitend), signalav((fittimeind:fitend)),'rx', 'LineWidth', 2);
plot(time(sdstart),signalav(sdstart),'go')
plot(time(rsdstart),signalav(rsdstart),'ro')
plot(time(maxInd),signalav(maxInd),'go')
plot(time(upstroke),signalav(upstroke),'go')
plot(time(downstroke),signalav(downstroke),'bo')
xlim([0 max(time)])
ax=gca;
axis tight
xx=get(ax,'Xlim');
yy=get(ax,'Ylim');

c='r';


baseline
hy = plot(xlim, [baseline baseline], 'Color', [.8 .8 .8]);

if startopt == 1
hx = plot([time(upstroke) time(upstroke)], ylim, 'Color', [.8 .8 .8]);
elseif startopt == 2
hx = plot([time(maxInd) time(maxInd)], ylim, 'Color', [.8 .8 .8]);
elseif startopt == 3
hx = plot([time(sdstart) time(sdstart)], ylim, 'Color', [.8 .8 .8]);
elseif startopt == 4
depol_mid=depol_mid*exposure;
hx = plot([depol_mid depol_mid], ylim, 'Color', [.8 .8 .8]);
end
ha = plot([Ti Ti], [yy(1) APD] ,'LineStyle', '--', 'Color', c); 
hb = plot([xx(1) Ti], [APD APD],'LineStyle', '--', 'Color', c);
plot(time(fittimeind:fitend),fitdecaymono(f.a,f.b,time(fittimeind:fitend)),'r-','LineWidth',2);
%plot(time(fittimeind:fitend),fitdecaymono(f.c,f.d,time(fittimeind:fitend)),'g-','LineWidth',2);
text(0.5,(0.95),['APD', num2str(t),': ',num2str(apd),'ms'],'Units','normalized');
text(0.5,(0.85),['Tau: ',num2str(tau_early),'ms'],'Units','normalized');
end

if ischar(colopt) == 1
    hold on
signalav=signalav-min(signalav);
if length(time) > length(signalav)
    time=time(1:length(signalav));
elseif length(signalav) > length(time)
    signalav=signalav(1:length(time));
end

plot(time, signalav, colopt, 'LineWidth', 2);
xlim([0 max(time)])
ax=gca;
axis tight
end

