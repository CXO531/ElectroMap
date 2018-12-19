
function [apd,tau,DI,ttp,maxvelup] = segapd(apdopt,tauopt,DIopt,ttpopt,maxupvelopt,Fluo,framerate,f,before,after,tfilt,t,startopt,apdblopt,apdblnum,tstar,tend);
% Function for doing analysis of paramaters measured in single file analysis GUI
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary


%% Normalise option
Fluo=imcomplement(Fluo);
Fluo=Fluo-min(Fluo);
Fluo=Fluo./max(Fluo);
order=3;framesize=11;
exposure = 1/framerate;
time=[1:(before+after)]*exposure%in milliseconds
before=round(before/exposure);
after=round(after/exposure);
baseline=[];

if f(1) < before
    startloc =2;
else startloc =1
end
if f(end) + after > numel(Fluo) 
    endloc=1;
else endloc = 0;
end
locRange = startloc:(numel(f)-endloc);
%% Overlay

for x = -before:after
    overlay(x+before+1) = sum(Fluo(f(locRange)+x));
end

%% Filter
 if tfilt == 3
            d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
            signalav=filtfilt(d,overlay);
 end
 
 if tfilt == 2
            signalav = sgolayfilt(overlay, order,framesize);
 end
            
 %% Find apd, ttp, and maxupvel
 % APD
 if apdopt == 1 || ttpopt == 1 || maxupvelopt == 1
        try
            % APD
            APD = t/100;
            signalav = (double(signalav));
            % can normalise here but taken out, explained below

            dsigav = diff(signalav);
            dsigav_up=dsigav(1:before+round(20/exposure)); %if two beats present (e.g alternan) will find first dpol (20ms buffer)
            
          
            %Possible Refrence points
            % 1, max upstroke (upstroke)
            [maxvelup, upstroke] = max(dsigav_up);
            
            % 2, Peak (maxInd)
            [maxval, maxInd] = max(signalav(1:(upstroke+round(25/exposure)))); %peak assotiead with first upstroke shortley after maxupstroke (set here to 20ms)
            
            % 3, Depol point (sdstart) (should be changed to max d2f/dt2???) 
            for i =1:upstroke
                dpol(i) = signalav(i);
            end
            %New dpol find code (!! Do comparsion of these for thesis!)
            ds=diff(dpol);
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
                      
            % 4, Depol midpont
            mini=signalav(sdstart);  
            if isempty(mini) == 1
                mini=min(signalav(1:upstroke));  
            end
            maxi=maxval;
            midi=(maxi-mini)*0.5;
            midi=midi+mini;
            for i = 1:before+round(50/exposure)  
            if signalav(i) < midi && signalav(i+1) > midi
                lowVal= signalav(i);
                highVal=signalav(i+1);
            end
            end
            
                       % Determines points for line equations
            if isempty(highVal) == 0 && isempty(lowVal) == 0
            midi;
            y1 = highVal;
            y2 = lowVal;


            x1 = find(signalav==highVal);
            x2 = find(signalav==lowVal);
            x1=x1;
            x2=x2;

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
            
            
            %[maxval, maxInd] = max(signalav);
            % Basline calc

            if apdblopt == 1
            blsec=(round(apdblnum/exposure));
            BLval = signalav(1:blsec); 
            baseline = mean(BLval);
            end
            
            if apdblopt == 2
            blsec=(round(apdblnum/exposure));
            BLval = signalav((length(signalav)-blsec):length(signalav)); 
            baseline = mean(BLval);
            end
            
            if apdblopt == 3
            baseline = min(signalav(1:upstroke));
            end
            
               
            if apdblopt == 4
                aftsig=signalav(maxInd:length(signalav));
                baseline = min(aftsig);
            end
            % Calculate APD
            
            APD = (maxval-baseline)*(1-APD)+baseline;
            
            checkSignal = signalav(maxInd:end); %checkSignal(26)=190;
            [~,min2]=min(checkSignal);
            checkSignal=checkSignal(1:min2); %ignore 2nd beat if present
            % [~, minInd] = min(abs(checkSignal - APD70));
            minInd = find(abs(checkSignal-APD)==min(abs(checkSignal-APD)));

            % Locates points above and below APD
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
            x1 = find(signalav==highVal);
            x2 = find(signalav==lowVal);
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
            
                     
            if startopt == 1
            apd = Ti-(upstroke);
            elseif startopt == 2
            apd = Ti-(maxInd);
            elseif startopt == 3
            apd = Ti-(sdstart);
            elseif startopt == 4
            apd = Ti-(depol_mid); 
            end
            
            apd=apd*exposure;
            ttp=(maxInd-sdstart)*exposure;
            maxvelup=(maxvelup)/exposure; %divide as changes units from Fluo/frame to Fluo/ms
            
        catch error
           apd = NaN;
           ttp = NaN;
           maxupvel = NaN;
        end
 else  apd = NaN;
       ttp = NaN;
       maxupvel = NaN;
 end

%% Tau
if tauopt == 1
               
      if isempty(baseline) == 1
            if apdblopt == 1
            blsec=(round(apdblnum/exposure));
            BLval = signalav(1:blsec); 
            baseline = mean(BLval);
            end
            
            if apdblopt == 2
            blsec=(round(apdblnum/exposure));
            BLval = signalav((length(signalav)-blsec):length(signalav)); 
            baseline = mean(BLval);
            end
            
            if apdblopt == 3
            baseline = min(signalav(1:upstroke));
            end
            
               
            if apdblopt == 4
                aftsig=signalav(maxInd:length(signalav));
                baseline = min(aftsig);
            end
      end
   [maxi, maxInd] = max(signalav);
   checkSignal = signalav(maxInd:end); %checkSignal(26)=190;
   taustart=(maxi-baseline)*((1-(tstar/100)));
   taustart=taustart+baseline;
   minInd = find(abs(checkSignal-taustart)==min(abs(checkSignal-taustart)));
   if minInd == length(checkSignal)
                minInd = minInd-1;
            end
            if(checkSignal(minInd) < taustart)  %only care about lowval for tau
                lowVal = checkSignal(minInd-1);
            else
                lowVal = checkSignal(minInd);
            end

            fittimeind = find(signalav==lowVal);
            
            % Find apd90 and stop exp fit from this point onwards
            
            APD90 = (maxi-baseline)*(1-(tend/100))+baseline;
            min90 = find(abs(checkSignal-APD90)==min(abs(checkSignal-APD90)));
            if min90 == length(checkSignal)
                min90 = minInd-1;
            end
           
              if(checkSignal(min90) > APD90)  
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
            fitdecay=@(a,b,x) a*exp(b*x);
            tau=(-1/f.b);
            end
      
      
else tau =NaN;
end


%% Diastolic Interval
if DIopt == 1
else DI=NaN;
end


if maxupvelopt == 1
else maxvelup=NaN;
end
% figure,
% hold on
% plot(time, signalav,'b', 'LineWidth', 1);
% plot(time, signalav,'kx', 'LineWidth', 2);
% plot(time(sdstart),signalav(sdstart),'go')
% plot(time(maxInd),signalav(maxInd),'go')
% plot(time(upstroke),signalav(upstroke),'go')
% plot(time(fittimeind:fitend),fitdecay(f.a,f.b,time(fittimeind:fitend)),'r-','LineWidth',2);
% plot(time(fittimeind:fitend), signalav((fittimeind:fitend)),'rx', 'LineWidth', 2);
% xlim([0 max(time)])
% ax=gca;
% axis tight
% xx=get(ax,'Xlim');
% yy=get(ax,'Ylim');
% if startopt == 1
% hx = plot([time(upstroke) time(upstroke)], ylim, 'Color', [.8 .8 .8]);
% elseif startopt == 2
% hx = plot([time(maxInd) time(maxInd)], ylim, 'Color', [.8 .8 .8]);
% elseif startopt == 3
% hx = plot([time(sdstart) time(sdstart)], ylim, 'Color', [.8 .8 .8]);
% elseif startopt == 4
% %depol_mid=depol_mid*exposure;
% hx = plot([depol_mid depol_mid], ylim, 'Color', [.8 .8 .8]);
% end
% ha = plot([Ti Ti], [yy(1) APD] ,'LineStyle', '--', 'Color', 'r'); 
% hb = plot([xx(1) Ti], [APD APD],'LineStyle', '--', 'Color', 'r');
% hy = plot([baseline baseline], ylim, 'Color', [.8 .8 .8]);
 