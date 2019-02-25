
function [map,meann,alll,onedev,vari,SE] = mapsbaby(startopt,framerate,t,maskedimage,imagestack,avbeat,outs,cmin,cmax,tfilt,before,apdblopt,apdblnum,medianfilter)
% function for taking a single beat image file and making apd/cad maps
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary



[rows, cols] = size(imagestack(:,:,1));

 if tfilt == 3
    d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
 end

baseline=zeros(size(maskedimage));
premap = zeros(rows,cols);
singalall=fluo_map(framerate,maskedimage,imagestack,tfilt,avbeat);
singalall=(singalall-min(min(singalall)))/(max(max(singalall))-min(min(singalall)));
exposure = 1/framerate; %in milliseconds
before=round(before/exposure);

order =3;
framesize =11;
for row = 1:rows
    for col = 1:cols
        if maskedimage(row,col) ~= 0
        try
            APD = t/100;
            signalav = imcomplement(squeeze((avbeat(row,col,:))));
            signalav = (double(signalav));
            % can normalise here but taken out, explained below
%             
%             if row==40 && col==40
%                 figure,
%                 plot(signalav)
%             end
            
            if tfilt == 2
            signalav = sgolayfilt(signalav, order,framesize);
            end
            
            if tfilt == 3
            signalav=filtfilt(d,signalav);
            end
%             
%             if row==40 && col==40
%                 hold on
%                 plot(signalav)
%                 figure,
%             end
            
            dsigav = diff(signalav);
            dsigav_up=dsigav(1:before+round(20/exposure)); %if two beats present (e.g alternan) will find first dpol (20ms buffer)
            
          
            %Possible Refrence points
            % 1, max upstroke (upstroke)
            [~, upstroke] = max(dsigav_up);
            
            % 2, Peak (maxInd)
            [maxval, maxInd] = max(signalav(1:(upstroke+round(25/exposure)))); %peak assotiead with first upstroke shortley after maxupstroke (set here to 20ms)
            
            % 3, Depol point (sdstart) (should be changed to max d2f/dt2???) 
            dpol=signalav(1:upstroke);
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
%             for i = 1:before+round(50/exposure)  
%             if signalav(i) < midi && signalav(i+1) > midi
%                 dlowVal= signalav(i);
%                 dhighVal=signalav(i+1);
%             end
%             end
            ind1=find(signalav > midi);
%             ind1
%             if isempty(sdstart) == 0
%             ind1=ind1(ind1>sdstart);
%             end
            if isempty(ind1) == 1
                ind1 = 1;
            end
            ind1=ind1(1);
            ind2=ind1-1;
            if ind2 == 0
                ind1 = 2;
                ind2 = 1;
            end
            dlowVal=signalav(ind2);
            dhighVal=signalav(ind1);
            

            % Determines points for line equations
            if isempty(dhighVal) == 0 && isempty(dlowVal) == 0

            dy1 = dhighVal;
            dy2 = dlowVal;


            dx1 = find(signalav==dhighVal);
            dx2 = find(signalav==dlowVal);

            if numel(dx2) > 1
            dx2=dx2(numel(dx2));
            end
            if numel(dx1)>1  
            dx1=dx1(1);
            end
            dm = (dy2-dy1)/(dx2-dx1);
            % Line constant, should be same for both c1 and c2
            dc1 = dy1-(dm.*dx1);
            dc2 = dy2-(dm.*dx2);

            % Time 
            depol_mid = (midi-dc1)/dm;
            if isempty(depol_mid) == 1 || isnan(depol_mid) == 1

            end
            end
            
            
            %[maxval, maxInd] = max(signalav);
            % Basline calc

            if apdblopt == 1
            blsec=(round(apdblnum/exposure));
            BLval = signalav(1:blsec); 
            baseline(row,col) = mean(BLval);
            end
            
            if apdblopt == 2
            blsec=(round(apdblnum/exposure));
            BLval = signalav((length(signalav)-blsec):length(signalav)); 
            baseline(row,col) = mean(BLval);
            end
            
            if apdblopt == 3
            baseline(row,col) = min(signalav(1:upstroke));
            end
            
               
            if apdblopt == 4
                aftsig=signalav(maxInd:length(signalav));
                baseline(row,col) = min(aftsig);
            end
            % Calculate APD
            
            APD = (maxval-baseline(row,col))*(1-APD)+baseline(row,col);
            
            checkSignal = signalav(maxInd:end); %checkSignal(26)=190;
            [~,min2]=min(checkSignal);
            checkSignal=checkSignal(1:min2); %ignore 2nd beat if present
            % [~, minInd] = min(abs(checkSignal - APD70));
            % minInd = find(abs(checkSignal-APD)==min(abs(checkSignal-APD)));
            minInd = find(checkSignal<APD,1);
            % Locates points above and below APD
                highVal = checkSignal(minInd-1);
                lowVal = checkSignal(minInd);


            % Solution for detecting a slight positive gradient next to APD region
            if (highVal - lowVal)<0
                %highVal = checkSignal(minInd-2);
                lowVal = find(checkSignal<APD,2);
%                 disp('**There was a positive gradient found APD70**');
            end

            % Determines points for line equations
            y1 = highVal;
            y2 = lowVal;

%             x1 = time(find(signalav==highVal)); % 'time' no value here
%             x2 = time(find(signalav==lowVal));
            x1 = maxInd+minInd;
            x2 = maxInd+minInd-1;
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
            if isnan(apd) == 1
            end
            end
            
            
            
            if row == 30 && col == 31
%                 figure,
%                 hold on
%                 plot(signalav,'x')
%                 plot(x1,APD,'ks')
%                 plot(x1,signalav(x1),'ob')
%                 plot(x2,signalav(x2),'og')
%                 plot(upstroke,signalav(upstroke),'or')
%                 figure,
            end
            
        catch error
           apd = 0;
        end
        if isempty(apd) == 1
            apd=0;
        end
        premap(row,col) = apd;
       end
    end
end
% figure,

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
    if isempty(alll) ~= 1
    alll=deleteoutliers(alll);
    end
end
onedev=std(alll);
SE=onedev/sqrt(numel(alll));
vari=var(alll);
% 
% min(alll)
% max(alll)
% pause(10)

meann = mean(alll);


premap=premap*exposure;
if exist('medianfilter') == 0
    medianfilter=1;
end
if medianfilter == 1
map=medfilt2(premap);
else
map=premap;
end
%map(map<=0)=0;




