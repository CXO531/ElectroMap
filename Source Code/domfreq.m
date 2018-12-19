function[freqmap] = domfreq(mask,imagestack,framerate,minf,maxf,fbin,winopt,tfilt)

% Function for creating DF map from an image stack. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary


wb=waitbar(0.1,'Calculating Dominant Frequncies');
[rows,cols,num] = size(imagestack(:,:,:));
freqmap=zeros(rows,cols);
imagestack=double(imagestack);
framerate=framerate*1000;
T=(1/framerate);
Fs=framerate;
L=num;
t=(0:L-1)*T;
order = 3;
framesize = 11;
if winopt == 0
winvec=ones(num,1)
end
if winopt == 1
winvec=hann(num)';
end

padf=(Fs)/(fbin*num)
lp = padf*(num);
lpa = nextpow2(lp);
lpad=2.^lpa;

 if tfilt == 3
    d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
 end

for row = 1:rows
    waitbar(0.1+(0.9*(row/rows)),wb)
    for col = 1:cols
        if mask(row,col)~=0
        pixelsignal=squeeze(imagestack(row,col,:))';
        %pixelsignal=imcomplement(pixelsignal);
        if tfilt == 2
        pixelsignal = sgolayfilt(pixelsignal, order,framesize);
        end
        
        if tfilt == 3
            pixelsignal=filtfilt(d,pixelsignal);
        end
        pixelsignal = pixelsignal - mean(pixelsignal); 
        pixelsignal=pixelsignal.*winvec;
        pixelsignal=pixelsignal(1,:);

        Y=fft(pixelsignal,lpad);
        Y=Y(1:lpad/2+1);
        P2=abs(Y/L);
        P1 = P2;
        P1(2:end-1) = 2*P1(2:end-1);
        f = 0:(Fs/lpad):(Fs/2);
        [o,oi]=find(f<maxf);
        [u,ui]=find(f>minf);
        [maxpower,maxind]=max(P1(min(ui):max(oi)));
        [~,xx]=find(P1==maxpower);
        xx=xx(1);
        freqmap(row,col)=f(xx);
        else
           freqmap(row,col)=NaN;
        end
    end
end
delete(wb)

freqmap=double(freqmap);
mask=double(mask);
freqmap=freqmap.*mask;


