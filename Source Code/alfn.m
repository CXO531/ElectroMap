function [normsig,sigPks,siglocs,bPks,blocs,APDs,...
          ampalt,loadalt,apdalt,...
          smalllocs,biglocs,smallPks,bigPks,smallloads_t,bigloads_t,smallloads,bigloads] = alfn(Fluo,framerate,Pks,locs,before,after,tfilt,apd);
% Function for calcualting alternan paramaters from single pixel trace 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary
      
      smalllocs=[];
biglocs=[];
smallPks=[];
bigPks=[];
smallloads_t=[];
bigloads_t=[];
smallloads=[];
bigloads=[];

%% Normalise option
Fluo=double(Fluo);
Fluo=imcomplement(Fluo);
Fluo=Fluo-min(Fluo);
Fluo=Fluo./max(Fluo);
order=3;framesize=11;
exposure = 1/framerate; %in milliseconds
before=round(before/exposure);
after=round(after/exposure);

%% Find minimum amplitude before peaks
bPks=zeros(length(Pks)-1,1);
blocs=zeros(length(locs)-1,1);
bPks(1)=0; blocs(1)=1;

for i = 1: length(Pks)
    length(Fluo);
    if locs(i)-before > 0 && locs(i)+after < numel(Fluo)
    signal=Fluo((locs(i)-before):(locs(i)+after));
    end
    [sigPks(i), siglocs(i)] = max(signal);
    before_sig=signal(1:siglocs(i));
    [bPks(i), blocs(i)]=min(before_sig);
    siglocs(i)=siglocs(i)+(locs(i)-before); %shift times to match signal
    blocs(i)=blocs(i)+(locs(i)-before);
end
hold on
Load=bPks;
Peak=sigPks;
sigPks=double(sigPks);
Amplitude=sigPks-bPks;
%plot(Time(blocs),bPks,'og','MarkerSize',3,'LineWidth',10);

%% Finding Amplitude Alternans and the effect of load

%Find the small and large beats
s_count=0;b_count=0;
for i = 1:length(Pks)-2
   if Peak(i+1) < Peak (i+2)
      s_count=s_count+1;
      smalllocs(s_count)=siglocs(i+1);
      smallPks(s_count)=sigPks(i+1);
      smallloads(s_count)=bPks(i+1);
      smallloads_t(s_count)=blocs(i+1);
   elseif Peak(i+1) >= Peak (i+2)
      b_count=b_count+1;
      biglocs(b_count)=siglocs(i+1);
      bigPks(b_count)=sigPks(i+1);
      bigloads(b_count)=bPks(i+1);
      bigloads_t(b_count)=blocs(i+1);
   end
end

for i = 2:length(sigPks)-1
    if sigPks(i) < sigPks(i+1)
    %Alternan_Ratio(i)=1-(sigPks(i)/sigPks(i+1));
    sigPks;
    Alternan_Ratio(i)=(sigPks(i+1)-sigPks(i))/sigPks(i);
    ampalt(i)=((sigPks(i+1))/sigPks(i))-1;
    elseif sigPks(i) >= sigPks(i+1)
    %Alternan_Ratio(i)=1-(sigPks(i+1)/sigPks(i));
    Alternan_Ratio(i)=(sigPks(i)-sigPks(i+1))/sigPks(i+1);
    ampalt(i)=(sigPks(i)/sigPks(i+1))-1;
        else
            ampalt(i)=0;
    end
    ampalt(i)=ampalt(i)*100;
end

for i = 2:length(sigPks)-1
    if (sigPks(i)-bPks(i)) < (sigPks(i+1)-bPks(i+1))
    Load_Alternan_Ratio(i)=((sigPks(i+1)-bPks(i+1))-(sigPks(i)-bPks(i)))/(sigPks(i)-bPks(i));
    Load_Index(i)=Load_Alternan_Ratio(i)-Alternan_Ratio(i);
    end
    if (sigPks(i)-bPks(i)) >= (sigPks(i+1)-bPks(i+1))
    Load_Alternan_Ratio(i)=((sigPks(i)-bPks(i))-(sigPks(i+1)-bPks(i+1)))/(sigPks(i+1)-bPks(i+1));
    Load_Index(i)=abs(Load_Alternan_Ratio(i)-Alternan_Ratio(i)); 
    end
end

%% APD alternans
q=apd/100;

for i=2:length(sigPks)
    try
    signalav=Fluo(locs(i)-before:locs(i)+after);
            if tfilt == 2
            signalav = sgolayfilt(signalav, order,framesize);
            end
    [maxval,maxInd] = max(signalav);
    amp=sigPks(i)-bPks(i);
    APD = (amp*(1-q))+bPks(i);
    
            checkSignal = signalav(maxInd:end); 
            [~,min2]=min(checkSignal);
            checkSignal=checkSignal(1:min2); %ignore 2nd beat if present
            minInd = find(abs(checkSignal-APD)==min(abs(checkSignal-APD)));
            
            if(checkSignal(minInd) > APD)
                highVal = checkSignal(minInd);
                lowVal = checkSignal(minInd+1);
            else
                highVal = checkSignal(minInd-1);
                lowVal = checkSignal(minInd);
            end
            y1 = highVal;
            y2 = lowVal;
            x1 = find(signalav==highVal);
            x2 = find(signalav==lowVal);
            x1=x1(end);
            x2=x2(1);
            m = (y2-y1)/(x2-x1); 
            c1 = y1-(m.*x1);
            Ti = (APD-c1)/m;
            APDs(i)=(Ti-maxInd)/exposure; 
            
    catch
        APDs(i)=NaN;
    end
end    

for i = 2:length(sigPks)-1
      if isnan(APDs(i)) == 0 && isnan(APDs(i+1)) == 0
      APD_Alternan_Ratio(i)=abs(APDs(i)-APDs(i+1));
      else APD_Alternan_Ratio(i) = 0;
      end
%     if isnan(APDs(i)) == 0 && isnan(APDs(i+1)) == 0
%     if APDs(i) < APDs(i+1)
%     APD_Alternan_Ratio(i)=1-(APDs(i)/APDs(i+1));
%     end
%     if APDs(i) >= APDs(i+1)
%     APD_Alternan_Ratio(i)=1-(APDs(i+1)/APDs(i));
%     end
%     else APD_Alternan_Ratio(i) = 0;
%     end
end

%% Shift method
% 
% for i = 2:length(Pks)-5
% peak1=Fluo(blocs(i):blocs(i)+after);
% peak2=Fluo(blocs(i+1):blocs(i+1)+after);
% di=peak1(:)-peak2(:);
% %Int(i) = trapz(1:after+1, abs(di));
% end
% 

%outputs
normsig=Fluo';
bPks=bPks(2:end)';
blocs=blocs(2:end)';
sigPks=sigPks(2:end)';
siglocs=siglocs(2:end)';
APDs=APDs(2:end)';
ampalt=ampalt(2:end)';
loadalt=Load_Alternan_Ratio(2:end)'*100; %make them %
apdalt=APD_Alternan_Ratio(2:end)';
% shiftalt=45;


