function [preimages,images,averages,mask] = OMimprocess(fname,im,rect,num_images,cropchoice,mask,sfilt,sfiltsize,inversion,tfilt,removef,camopt,sfiltsigma,pbefore,pafter)
%funtion for reading in and processimg all images in the tif stack/ rhs
[rows cols] = size(im);
[token,remain] = strtok(fname,'.');
fileisrsh=0;
if strcmp(remain, '.mat') == 1
    fileisrsh = 2;
    S=whos('-file',fname);
    imvar=S.name;
    X=load(fname);
    v=struct2cell(X);
    images=cell2mat(v);
    images=double(images);
    images=images-min(min(images));
    images=images./max(max(images));
    figure,
    imshow(images(:,:,1),[])
    images=images*((2^16)-1);
    images=uint16(images);
    dbs=16;
end

%% Load images and store data in array
wb=waitbar(0.2,'Processing Images');
averages = zeros(1, num_images);
sigma=sfiltsigma;
%spatial filters
if sfilt == 1
h = fspecial('average', 1);
end
if sfilt == 2
h = fspecial('gaussian', sfiltsize,sigma); 
end
if sfilt == 3
h = fspecial('average', [sfiltsize sfiltsize]);
end

%image load and spatial filter
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')

wbfactor=floor(num_images/6);
wbcount=0;



if fileisrsh == 0
    images=zeros(rows,cols,num_images);
preimages=zeros(rows,cols,num_images);
    TifLink = Tiff(fname, 'r');
for j  = 1:num_images;  
    TifLink.setDirectory(j);
    A=TifLink.read();
    if cropchoice == 1
    A = imcrop(A, rect);
    end
    if inversion == 0
       A=imcomplement(A);
    end
    preimages(:,:,j)=A;
    A = imfilter(A, h);
    A = uint16(A);
    mask=uint16(mask);
    B = A.*mask;
    images(:,:,j) = B;
    averages(:,j) = mean2(B);
    if mod(j,wbfactor) == 0
       wbcount=wbcount+1;
    waitbar(0.3+0.1*wbcount,wb,'Processing Images');
    end
end
TifLink.close();
end
%averages(:,1:20)=ones(20,1)*averages(:,21);

if fileisrsh == 2;
for j  = 1:num_images;
    A=images(:,:,j);
    if cropchoice == 1
    A = imcrop(A, rect);
    end
    if inversion == 0
       A=imcomplement(A);
       A=A-min(A);
    end
    preimages(:,:,j)=(A);
    A = imfilter(A, h);
    mask=uint16(mask);
    B = A.*mask;
    images(:,:,j) = B;
    averages(:,j) = mean2(B);
    if mod(j,wbfactor) == 0
       wbcount=wbcount+1;
    waitbar(0.3+0.1*wbcount,wb,'Processing Images');
    end
end
end

delete(wb)

if fileisrsh == 1
averages  = imcomplement(averages);
end

if fileisrsh == 0;
averages  = imcomplement(averages);
end
tav=averages-min(averages);


%% Frame removal
if removef == 1
    disp('hi')
    pfind=figure
    plot(tav)
    [opks,olocs]=findpeaks(tav,'MINPEAKHEIGHT',max(tav)*0.75,'MINPEAKDISTANCE',3);
    hold on
    plot(olocs,opks,'or')
    waitfor(pfind)
    before=pbefore;
    after=pafter;
    pfcount=0;
    pframes=0;
    % find pulse frames
    for i=1:length(olocs)
        for k=-before:after
            pframe=[olocs(i)+k];
            pcheck=ismember(pframes,pframe);
            if any(pcheck) == 0
                pfcount=pfcount+1;
                pframes(pfcount)=pframe;
            end
            
        end
    end
    
    %remove pulse frames
    for j=1:length(pframes)
        pframes(j)
        if pframes(j) > 1
        C=images(:,:,pframes(j)-1);
        images(:,:,pframes(j))=C;
        end
    end
    images(:,:,1)=images(:,:,2);
    % reclac avs
    averages=[];
    for i=1:length(images)
    D=images(:,:,i);
    averages(:,i) = mean2(D);
    end
    averages=imcomplement(averages);
end

if tfilt ==3 %filter averages
d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
averages=filtfilt(d,averages);
end
                                                        
if tfilt == 2 || tfilt == 1
averages = sgolayfilt(averages, 3,11);
end



