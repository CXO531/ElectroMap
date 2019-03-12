clear all

%% NEW RHS FILE CONVERTER
% 
[fname,pathname]=uigetfile('*.rsh')
filename=[pathname,fname]
fid=fopen(filename,'r','b')
fstr=fread(fid,'int8=>char'); %string of chars held in .rsh file
fstr=fstr';
fclose(fid);

fstr;

%% get info
k = strfind(fstr,'DataXsize='); 
xs=str2num(fstr(k+11:k+13)); %x
k = strfind(fstr,'DataYsize='); 
ys=str2num(fstr(k+11:k+13)); %y
k = strfind(fstr,'sample_time='); 
exposure=str2num(fstr(k+13:k+16)); 
framerate=1/exposure;
k = strfind(fstr,'frame_number='); 
num=str2num(fstr(k+13:k+16));
k = strfind(fstr,'dual_cam'); dual = str2num(fstr(k+9));
k = strfind(fstr,'Data-File-List'); filestart=k+16;


%% find data files

file_string=fstr(filestart:end);
filelist={};
[token,remain] = strtok(file_string,'.');

for i=1:100000
    str1=file_string;
    str2=[token,'(',num2str(i-1),').rsd'];
    k = findstr(str1,str2);
    if isempty(k) == 0
    filelist{i}=[token,'(',num2str(i-1),').rsd'];
    else
         num_of_files =i-1;
    break
    end
end


%% load data
k=1;
for i=1:num_of_files
    fid=fopen([pathname,filelist{i}],'r','l');
    fdata=fread(fid,inf,'int16=>int32');
    fclose(fid);
    fdata = reshape(fdata,12800,[]);
    step = 1;
    if dual ~= 0
        step=2;
    end
    for j = 1:step:size(fdata,2);
                   if dual == 0 
                    framedata = fdata(:,j);
                    framedata = reshape(framedata,128,100);
                    cmosData(:,:,k) = framedata(:,:);
                    k=k+1;
                   end
                    if dual ~= 0
                    oneframe = fdata(:,j);
                    oneframe = reshape(oneframe,128,100);
                    cmosData(:,:,k) = oneframe(:,:);
                    oneframe2 = fdata(:,j+1);
                    oneframe2 = reshape(oneframe2,128,100);
                    cmosData2(:,:,k) = oneframe2(:,:);
                    k=k+1;   
                    end
    end
            clear fdata;
end


%% Get bg frame if new camera
    if strcmp(fstr(3), 'M')
        fid = fopen([pathname,token,'.rsm'], 'r', 'l');
        fdata = fread(fid, 'int16=>32')';
        fclose(fid);
        fdata = int32(reshape(fdata, 12800, []));
        fdata = reshape(fdata, 128, 100);
        bgimage = fdata(21:120, :);
        bgimage2 = fdata(21:120, :);
    else 
        bgimage = cmosData(21:120, :, 1);
         if dual ~= 0
        bgimage2 = cmosData2(21:120, :, 1);
         end
    end
    
    
    
    
    

%% trim and save data

images=cmosData(21:120,:,:);
images1 = images - repmat(bgimage, [1, 1, size(images, 3)]);
images1=imcomplement(images1);
images1=images1-min(min(min(images1)));
 if strcmp(fstr(3), 'U')
     images1(:,:,1)=bgimage;
 end
if dual > 0
images2=cmosData2(21:120,:,:);
images2 = images2 - repmat(bgimage2, [1, 1, size(images, 3)]);
images2=imcomplement(images2);
images2=images2-min(min(min(images2)));
 if strcmp(fstr(3), 'U')
     images2(:,:,1)=bgimage2;
 end
save([token,'.mat'],'images1','images2');
else
save([token,'.mat'],'images1');
end
