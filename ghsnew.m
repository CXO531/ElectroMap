clear all
clc

%% NEW GSH FILE TESTER

[fname,pathname]=uigetfile('*.gsh')
filename=[pathname,fname]
fid=fopen(filename,'r','b')
fstr=fread(fid,'int8=>char'); %string of chars held in .rsh file
fstr=fstr';
fclose(fid);

fstr;

%% get info from header file
k = strfind(fstr,'Data Size =')
xs=str2num(fstr(k+13:k+15)); %x
ys=str2num(fstr(k+17:k+19)); %y
k = strfind(fstr,'sample_time='); 
exposure=str2num(fstr(k+13:k+16)); framerate=1/exposure;
k = strfind(fstr,'frame_number='); 
num=str2num(fstr(k+13:k+16));
k = strfind(fstr,'dual_cam'); dual = str2num(fstr(k+9));
[token,remain] = strtok(fname,'.');
filelist{1}=[token,'.gsd'];
 
num_of_files = 1;
k=1;
    fid=fopen([pathname,filelist{1}],'r','l');
    fseek(fid,256,-1); %move to start of deets
    ghsdata=fread(fid, 6, 'int16');
    xskip=ghsdata(3);
    yskip=ghsdata(4);
    xpixs=ghsdata(5);
    ypixs=ghsdata(6);
    
    % Background Image
    fseek(fid,972,-1); %move to background image
    bgimage = fread(fid, xs * ys, 'int16');
    bgimage = reshape(bgimage, [xs, ys]);
    bgimage = bgimage(xskip+1:xskip+xpixs,yskip+1:yskip+ypixs); %trim
    
    %Data
    fdata=fread(fid,xs*ys*num,'int16');
    fdata=reshape(fdata,[xs,ys,num]);
    fdata=fdata(xskip+1:xskip+xpixs,yskip+1:yskip+ypixs,:);

    
    images=fdata;
    images1 = images - repmat(bgimage, [1, 1, size(images, 3)]);
    images1=imcomplement(images1);
    save([token,'.mat'],'images1');
    
    fclose(fid)