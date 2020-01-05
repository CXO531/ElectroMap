function [num_images,rect,mask,im,I,boundaries,camopt,frame1,fluoim,rois,n,tstack] = OMimload(fname,cropchoice,quinnieopt,threshop,threshman,rect,inversion,camopt,fluo_opt,roinum,roisum,tstack)
% Function for taking raw image file and applying thersholding, final
% output - regradless of input, will be uint16 mask. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 tstack=[];
% Last Updated -
 tic
% Update Summary
id='MATLAB:imagesci:tifftagsread:expectedTagDataFormat';
warning('off',id)
id='MATLAB:imagesci:tiffmexutils:libtiffWarning';
warning('off',id)
%% Work out filetype. 0=tiff 2=mat
[qq,token,remain] = fileparts(fname);
dbs=0;
fileisrsh=0; %file is tiff stack 
n=[];
if strcmp(remain, '.mat') == 1 %file is .mat file 
    fileisrsh = 2;
    S=whos('-file',fname);
    imvar=S.name;
    X=load(fname);
    v=struct2cell(X);
    sv = size(v,1);
    if sv == 2
     qn = questdlg('Which Dataset?', ...
	'Which Dataset?', ...
	'1','2','1');
% Handle response
switch qn
    case '1'
        n = 1;
    case '2'
        n = 2;
end
       images=v{n};
       images=double(images);
       % errordlg('Multiple variables found in .MAT file!. Please edit file to contain only images variable')
    else
        v
    images=cell2mat(v);
    images=double(images);
    end
    frame1=images(:,:,1);
    frame1=uint16(frame1);
    images=images-min(images,[],3);
    images=images./max(max(max(images)));
    %plot(squeeze(images(25,40,:)),'b')
    images=images*((2^16)-1);
    images=uint16(images);
    %plot(squeeze(images(25,40,:)),'r')
    dbs=16;
end

%% Get info and normlaise to help thersholding
wb = waitbar(0,'Loading Images');
if fileisrsh == 0
info = imfinfo(fname);
%% work out file BitDepth
dbs=info.BitDepth;
num_images = numel(info);
num_images=num_images-1;
frame1=imread(fname,2);
%  TifLink = Tiff(fname, 'r');
%tsStack = TIFFStack(fname)



InfoImage=imfinfo(fname);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
FileID = tifflib('open',fname,'r');
rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
 
for i=1:NumberImages
   waitbar(i/NumberImages,wb,'Loading Images');
   drawnow()
   tifflib('setDirectory',FileID,i-1);
   % Go through each strip of data.
   rps = min(rps,nImage);
   for r = 1:rps:nImage
      row_inds = r:min(nImage,r+rps-1);
      stripNum = tifflib('computeStrip',FileID,r);
      tstack(row_inds,:,i) = tifflib('readEncodedStrip',FileID,stripNum-1);
   end
end
tifflib('close',FileID);







if num_images > 2000 %For speed, only do raw sig level from first 2000 frames - should be setting? 
    fluoims = 2000;
else
    fluoims = num_images;
end
for j  = 1:fluoims
%     TifLink.setDirectory(j);
%     A=TifLink.read();
 A=tstack(:,:,j);
    rawimages(:,:,j)=A;
    
end
fluoim=fluo_map_raw(rawimages);
if fluo_opt == 1
    im = frame1;
elseif fluo_opt == 2
    im= fluoim;
end
end


if fileisrsh == 2 %.mat
num_images=size(images,3);
if num_images > 2000
    fluoims = 2000;
else
    fluoims = num_images;
end
  for j  = 1:fluoims;  
    rawimages(:,:,j)=images(:,:,j);
  end
fluoim=fluo_map_raw(rawimages);
if fluo_opt == 1
    %check is any frames NaNs
    TF = isnan(frame1(1,1)); 
    if TF == 0  
    im = frame1;
    else 
        newTF=isnan(rawimages(1,1,:));
        k=find(newTF);
        im=(rawimages(:,:,k(1)));
    end
elseif fluo_opt == 2
    im=fluoim;
end
end



imvec=reshape(im,1,(size(im,1)*size(im,2)));
i3=im-min(imvec);
i3=double(i3);
i3=i3/double(max(imvec)-min(imvec)); 
%% image cropping
% square crop from user
if cropchoice == 1 && isempty(rect) == 1
figure,
imshow(im, [],'InitialMagnification', 800) 
cropfigure=gcf;
title('Make your selection and press enter');
[~,rect]=imcrop;% manually select ROI
hold on
rectangle('Position', rect, 'EdgeColor', 'r');
hold off
i3 = imcrop(i3, rect);
im = imcrop(im, rect);
fluoim=imcrop(fluoim,rect);
frame1=imcrop(frame1,rect);
close(cropfigure)
else
    rect=[];
end

%square crop saved from before
isempty(rect)
if cropchoice == 0 && isempty(rect) == 0
        disp('why are you here?')
i3 = imcrop(i3, rect);
im = imcrop(im, rect);
fluoim=imcrop(fluoim,rect);
frame1=imcrop(frame1,rect);
end
rois=zeros(size(im));
% manual roi crop
if quinnieopt == 1
figure,
imshow(i3, [],'InitialMagnification', 800) 
roifigure=gcf;
title('Select ROI and press enter'); 
if roinum == 1 
BWpoly=roipoly;
rois(:,:,1)=BWpoly;
end
if roinum == 2 || roinum == 3
   
BWpoly1=roipoly;
     polyfig=gcf;
 close(polyfig);
 figure,
 imshow(i3, [],'InitialMagnification', 800)
 hold on
 title('Select ROI and press enter'); 
   [B1,L,N,A] = bwboundaries(BWpoly1);
B1=B1{1};
   plot(B1(:,2), B1(:,1),...
       'r','LineWidth',2);
    BWpoly2=roipoly;
    BWpoly2=BWpoly2;
 polyfig=gcf;
 close(polyfig);
 rois(:,:,1)=BWpoly1;
 rois(:,:,2)=BWpoly2;
 BWpoly=BWpoly2+BWpoly1;
if roinum == 3
    figure,
 imshow(i3, [],'InitialMagnification', 800) 
 hold on
 [B2,L,N,A] = bwboundaries(BWpoly2);
 B2=B2{1};
  plot(B1(:,2), B1(:,1),...
       'r','LineWidth',2);
  plot(B2(:,2), B2(:,1),...
       'g','LineWidth',2);
 title('Select ROI and press enter'); 
    BWpoly3=roipoly;
    BWpoly3=BWpoly3;
    rois(:,:,3)=BWpoly3;
 polyfig=gcf;
 close(polyfig);
 BWpoly=(BWpoly)+BWpoly3;
end
end
if roisum == 0
    BWpoly(BWpoly>1)=0;
else
    BWpoly(BWpoly>1)=1;
end
close(roifigure)
end


%% Threshold image
if threshop == 1 %auto threshold 
    threshold=(graythresh(i3)*0.8)*max(imvec);
end

if threshop == 2 %manual threshold (percent change from auto)
    threshold=(graythresh(i3)*0.8)*(1-(-threshman/100))*max(imvec);
end

mask=[];
if dbs == 32
mask = uint32(im>=threshold);    
elseif dbs == 16
mask = uint16(im>=threshold);
elseif dbs == 8
mask = uint8(im>=threshold);
elseif dbs == 0
mask=double(im>=threshold);
end
% end

if quinnieopt == 1 %apply custom ROI
    BWpoly=uint16(BWpoly);
    mask=mask.*BWpoly;
end

close(wb)
class(im)
class(mask)
I=im.*mask;
[boundaries] = bwboundaries(mask);
class(tstack)
toc
