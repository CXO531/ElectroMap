function [locs,minimas,q2locs,avgCL,numofpeaksoverall,newpeakheight]=Omseg2(signal,peakheight,peakdist,minpeakheight,minpeakdist,minboundary,segchoice,minmumofpeaks,numimages,div,before,after)
% Function for segmenting a signal based on CL according to user settings
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary

%% Detect Peaks and minimums
normsig=signal-min(signal);
normsig=normsig./max(normsig);
insig=imcomplement(normsig);
insig=insig-min(insig);
insig=smooth(insig);
%find peaks
maxfluo=max(signal);
newpeakheight=maxfluo*peakheight;

%make start up to before and end up to after zero so no peaks can be found here
signal(1:before)=zeros(1,before);
signal(end-after+1:end)=zeros(1,after);

[pks, locs] = findpeaks(signal, 'MINPEAKHEIGHT',newpeakheight, 'MINPEAKDISTANCE', peakdist);
remove1=0;
removeend=0;

minimas=[];
%% Split up 'averages'
CL=zeros((length(locs)),3);
CL(:,1)=locs(:);
CLdiff=zeros(1,(length(CL)-1));
%Calculate CLs
if length(locs) >= 2
CL(1,2)=locs(2)-locs(1); %first peak as no refrence before, so use next peak
end
for i=2:(length(locs))
    CL(i,2)=locs(i)-locs(i-1);
end

for i=2:(length(locs)-2)
    CLdiff(i)=(locs(i+1)-locs(i))-(locs(i)-locs(i-1));
    if abs(CLdiff(i)) > minboundary
        CL(i,3)=1;
        if CL(i-1,3) == 1
            CL(i,3) = 0;
        end
    end
end

%Find diffrenet regions of constant CL , and position in signal.
q1=[];k=1;j=0; q2=[]; q2locs=[];
if isempty(locs) == 1 %no peaks
    h=errordlg('No action potential detected')
    waitfor(h);
end

if length(locs) > 1
    for i=1:(length(locs))
        if CL(i,3) == 0
            q1(k,i)=[CL(i,2)];
        else q1(k,i)=[CL(i,2)];
            k=k+1;
        end
    end
     
    for k=1:length(q1(:,2))
        numofpeaks=sum(q1(k,:)~=0);
        if numofpeaks >= minmumofpeaks
           j=j+1;
           q2(j,:)=q1(k,:);
        end
        numofpeaksoverall = numofpeaks;
    end
    
    if isempty(q2) ==1
        errordlg('Did not find enough peaks, please reduce minimum number required in Signal Processing Options. If signal possible arrhythmia suggest min number and min boundary both reduced to 1');
    end
     avgCL=[];
    for j=1:length(q2(1,:))
        for i=1:length(q2(:,1))
            if q2(i,j) ~=0
                q2locs(i,j)=locs(j);
            else
                q2locs(i,j)=0;
            end
        end
    end
    avgCL=zeros(2,length(q2(:,1)));
    for i=1:length(q2(:,1))
        pos=mean(find(q2(i,:)));
        %avgCL(1,i)=round(pos);
        avgCL(1,i)=(pos);
        avgCL(2,i) = mean(nonzeros(q2(i,:)));
    end
end    

%% check if first or last peak too close to end
% if q2locs(1,1) < bframes
%    q2locs=q2locs(:,2:end);
% end
% if q2locs(length(q2locs(:,1)),length(q2locs(1,:)))+aframes > num_images
%     q2locs=q2locs(:,1:end-1);
% end

%% Sub-Segmentation
if segchoice == 2 || segchoice == 3
    q2locsnew=[];
    k=0;col=0;row=1;
    for i=1:length(q2locs(:,1))
        for j=1:length(q2locs(1,:))
            if q2locs(i,j) ~= 0;
                k=k+1;
                colum=mod(k,div);
                
                if colum == 0
                    col=div;
                else col=colum;
                end
                
                q2locsnew(row,col)=q2locs(i,j);
                if j == length(q2(1,:)) %if loop to stop martix dim error
                    j=length(q2(1,:))-1;
                end
                
                
                if col == div
                    row=row+1;
                    k=0;
                end
                if j+1 < length(q2locs(1,:))
                    if q2locs(i,j+1)==0
                        row = row+2;
                    end
                end
            end
        end
    end
    
    
    if div > size(q2locsnew,2)
        div=size(q2locsnew,2)
    end
        
    q2locsnew=[q2locsnew;zeros(2,div)]; %added to allow 'buffer' at end to allow idneftication of last constat CL region.
    firstrow=1; A=[]; B=[];
    
    for i=1:length(q2locsnew(:,1))
        check=any((q2locsnew(i,:)));
        if check == 0
            A=q2locsnew(firstrow:i,:);
            A=sort(A,1);
            A=sort(A,2);
            B=[B;A];
            firstrow=i+1;
        end
    end
    
    q2locs=B(any(B,2),:);
    
     if segchoice == 3
        row =1; C=[];
        for i = 1:(length(B(:,1))-1)
            if any(B(i+1,:)) == 0
                C(row,:)=B(i,:);
                row = row+1;
            end
        end
        q2locs= C(any(C,2),:);
     end
     cld=zeros(length(q2locs(:,1)),length(q2locs(1,:))-1)
     %update avgcl for subseq
     for i=1:length(q2locs(:,1));
        if any((q2locs(i,:))) == 1
            pos=mean(find(q2locs(i,:)));
            %avgCL(1,i)=round(pos);
            avgCL(1,i)=pos;
            for j=1:(length(q2locs(1,:))-1) %%does nothing when 1 peak?
                disp('f')
                cld(i,j)=q2locs(i,j+1)-q2locs(i,j);
                if q2locs(i,j+1) == 0 || q2locs(i,j) == 0
                    cld(i,j)=0;
                end
            end
            
            if length(q2locs(1,:)) ~= 1
                avgCL(2,i) = mean(nonzeros(cld(i,:)));
            else
                if i < length(q2locs)
                    avgCL(2,i)= locs(i+1)-locs(i);
                else avgCL(2,i)=avgCL(2,i-1);
                end %to deal with last peak
            end
        end
    end
    if length((q2locs(1,:))) == 1 %%only one peak
        q=1;
        for i=1:length(CL(:,1));
            if q < length(q2locs)
                if CL(i,1) == q2locs(q); %%check this peak in settings
                    avgCL(2,q) = CL(i,2);
                    q=q+1;
                end
            end
        end
    end
    
    avgCL=avgCL(:,any(avgCL,1));  %columns
    %avloc=locs(avgCL(1,:))
end
