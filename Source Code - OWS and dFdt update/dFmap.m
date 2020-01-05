function [map,vf] = dFmap(images,mask,dFlevel,framerate,i_framerate,tfilt,invert,vfup,vfdown)

exposure=1/framerate; %in ms
[rows,cols,num] = size(images);
images=double(images);
order=3; framesize=5;

%% Normalise APs
Imax=max(images,[],3); %Find maximum value of each pixel
Imin=min(images,[],3); %Find minimum value of each pixel
interval=Imax-Imin;
images=(images-Imin)./interval; %% Normalise each pixel signal 0 to 1
if invert == 1
    images=imcomplement(images); %%Invert signal if required (i.e. when using a vm dye)
end

%% Interpolatation Settings
frames=[1:exposure:(num*exposure)]; %orignal frames
i_frames=[1:(1/i_framerate):(num*exposure)]; %new frames

%% calculate dFdt value at desired dFlevel (0-1)

for r=1:rows
    for c=1:cols
        if mask(r,c) ~= 0
            signal=squeeze(images(r,c,:));
%             if r == 25 && c == 25
%                 figure,
%                 plot(squeeze(signal'))
%                 %title(num2str(signal))
%             end
            if tfilt == 2
                signal = sgolayfilt(signal,order,framesize); %Temporal Filtering
            end
%            if r == 25 && c == 25
%                 hold on
%                 plot(squeeze(signal'))
%             end
% frames
% i_frames
% size(signal)
            signal=spline(frames,signal,i_frames); %Interpolation
            d_signal=diff(signal); %Differentiate signal
            if r == 15 && c == 15
            assignin('base','d_signal',d_signal.*(i_framerate));
            assignin('base','signal',signal);
            end
              if r == 18 && c == 15
            assignin('base','d_signal2',d_signal.*(i_framerate));
            assignin('base','signal2',signal);
            end
            if dFlevel ~= 0
                index=find(signal>dFlevel,1); %find first time amplitude is greater than dFlevel
                map(r,c)=d_signal(index); %find df/dt value at this point
                vf(r,c)=dFlevel;
            end
            if dFlevel == 0
                [maxval maxindex]= max(signal);
                d_signal.*(i_framerate);
                [val index]= max(d_signal(1:maxindex-1));
                map(r,c)=val;
                              if r == 18 && c == 18
            assignin('base','d_signal2',d_signal.*(i_framerate))
            assignin('base','signal2',signal)
            assignin('base','val',val)
            end
                if index + 81 < numel(signal)
                vf(r,c)=signal(index);
                else
                    vf(r,c)=signal(end);
                end
                
%                 if vf(r,c) < vfdown || vf(r,c) > vfup
%                     map(r,c)=NaN
%                 end
                
            end
        else
            map(r,c)=NaN;
            
        end
    end
end
framerate;
i_framerate;
map=map.*(i_framerate);
% map=medfilt2(map);
% vf=medfilt2(vf);
