
function [map,meann,alll,onedev,vari,SE] = ttpnew(t1,t2,framerate,t,maskedimage,imagestack,avbeat,outs,cmin,cmax,tfilt,before,apdblopt,apdblnum)
% Function for calculating ttp
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary


%Check Inputted fit paramaters make sense
if t1 < 0 || t1 > 100 || t2 < 0 || t2 > 100
    errordlg('Please check inputted ttp start and end values (both should be between 0 and 100)')
end

[rows cols] = size(imagestack(:,:,1))
counter = 0;

 if tfilt == 3
            d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
            end
tic;
maskedimage;
baseline=zeros(size(maskedimage));
premap = zeros(rows,cols);
singalall=fluo_map(framerate,maskedimage,imagestack,tfilt,avbeat);
singalall=(singalall-min(min(singalall)))/(max(max(singalall))-min(min(singalall)))
exposure = 1/framerate; %in milliseconds
before=round(before/exposure);
count = 0;
order =3;
framesize =11;
for row = 1:rows
    for col = 1:cols
        dpol=[];
        if maskedimage(row,col) ~= 0
        try
            % APD
            signalav = imcomplement(squeeze((avbeat(row,col,:))));
            signalav = (double(signalav));
            % can normalise here but taken out, explained below
            if tfilt == 2
            signalav = sgolayfilt(signalav, order,framesize);
            end
            
            if tfilt == 3
            signalav=filtfilt(d,signalav);
            end
            
            signalav=signalav-min(signalav);
            [maxi maxInd] = max(signalav);
            [mini minind] = min(signalav(1:maxInd));
            upstroke = signalav(minind:maxInd);
            upstroke = upstroke-min(upstroke);
            maxi=max(upstroke);
           %starttime
           a1=(t1*maxi)/100;
           ind_above_1=find(upstroke>a1);
           du=diff(ind_above_1);
           duind=find(du>1,1,'last');
           if isempty(duind) == 1
               duind = 0
           end
           ind_above_1=(ind_above_1(duind+1));
           length(upstroke);
           if ind_above_1 >1
           ind_below_1=ind_above_1-1;
           m1=(upstroke(ind_above_1)-upstroke(ind_below_1)); %don't need /x2-x1 as this is always 1
           c1=upstroke(ind_above_1)-(m1*ind_above_1);
           c2=upstroke(ind_below_1)-(m1*ind_below_1);
           a1_ind=((a1-c1)/m1);
           
%            linfit = polyfit([ind_below_1, ind_above_1], [upstroke(ind_below_1), upstroke(ind_above_1)], 1);
%            a=linfit(1);
%            b=linfit(2);
%            a1_ind=fzero(@(x) (a*x+b)-a1, [ind_below_1, ind_above_1]);
           else 
               a1_ind=1
           end
           
            %endttime
            
           b1=(t2*maxi)/100;
           ind_above_2=find(upstroke>b1,1);
           ind_below_2=ind_above_2-1;
           m2=(upstroke(ind_above_2)-upstroke(ind_below_2)); %don't need /x2-x1 as this is always 1
           c1=upstroke(ind_above_2)-(m2*ind_above_2);
           c2=upstroke(ind_below_2)-(m2*ind_below_2);
           b1_ind=((b1-c1)/m2);
           
           
           
%            
%            linfit = polyfit([ind_below_2, ind_above_2], [upstroke(ind_below_2), upstroke(ind_above_2)], 1);
%            a=linfit(1);b=linfit(2);
%            
%            b1_ind=fzero(@(x) (a*x+b)-b1, [ind_below_2, ind_above_2]);
           
           apd=b1_ind-a1_ind;
           
           
            
         catch error
             apd = NaN;
         end
        if isempty(apd) == 1
            apd=NaN;
        end
        premap(row,col) = apd;
       end
    end
end


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
sigcount=0
s=[0.25,0.5,0.75,0.1];
if outs == 9 || outs == 10 || outs == 11 || outs == 12
    siglevel=s(outs-8);
    alll=[];
    for r=1:rows
        for c=1:cols
            if singalall(r,c) > siglevel && isnan(premap(r,c)) ~= 1 && isinf(premap(r,c)) ~= 1 && premap(r,c) > 0 
                sigcount=sigcount+1;
                alll(sigcount)=premap(r,c)*exposure;
            else premap(r,c) = NaN;
            end
        end
    end
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
% get rid of -ve APD values in the map as well
for i=1:rows
    for j=1:cols
        if map(i,j) <= 0
            map(i,j) = 0;
        end
    end
end
map=map*exposure;


