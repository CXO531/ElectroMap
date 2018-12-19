function varargout = segEP(varargin)

% Main function for running single file analysis GUI. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary


% SEGEP MATLAB code for segEP.fig
%      SEGEP, by itself, creates a new SEGEP or raises the existing
%      singleton*.
%
%      H = SEGEP returns the handle to a new SEGEP or the handle to
%      the existing singleton*.
%
%      SEGEP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGEP.M with the given input arguments.
%
%      SEGEP('Property','Value',...) creates a new SEGEP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segEP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segEP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segEP

% Last Modified by GUIDE v2.5 13-Nov-2017 10:38:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segEP_OpeningFcn, ...
                   'gui_OutputFcn',  @segEP_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before segEP is made visible.
function segEP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segEP (see VARARGIN)

% Choose default command line output for segEP
handles.output = hObject;
h = findobj('Tag','ElectroMap');
g1data = guidata(h);


%stick activation map in there for lols
axes(handles.axes1);
wb=waitbar(0.5,'Transfering Activation Map');
[handles.amap,~,~,handles.act_t,~,~,~,~,~,~,handles.vout,handles.quivers_Xout,handles.quivers_Yout,handles.quivers_vxout,handles.quivers_vyout]...
    =cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),g1data.cvimages,g1data.mask,1,str2double(get(g1data.minvel,'String')),str2double(get(g1data.maxvel,'String')),get(g1data.velalgo,'Value'),...
    str2double(get(g1data.MINt,'String')),str2double(get(g1data.MAXt,'String')),str2double(get(g1data.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String'))); 
delete(wb)
jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [amap] = handles.amap;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(amap));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(amap, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);

imshow(g1data.im,[])
hold on
for i=1:size(g1data.boundaries,1) 
plot(g1data.boundaries{i}(:,2),g1data.boundaries{i}(:,1),'r','LineWidth',2);
end
%Set values from electromap
set(handles.APDs,'String',get(g1data.t,'String'));
set(handles.actims,'String',get(g1data.actmax,'String'));
set(handles.winsize,'String',get(g1data.winsize,'String'));
set(handles.minvel,'String',get(g1data.minvel,'String'));
set(handles.maxvel,'String',get(g1data.maxvel,'String'));
set(handles.mint,'String',get(g1data.MINt,'String'));
handles.apdmap=g1data.apdmap;
handles.apalll=g1data.apalll;
handles.actmap=g1data.actmap;

%set(handles.maxt,'String',get(g1data.MAXt,'String'));

%allow button press on image
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);

%points for single vector stuff
handles.pointA=[];
handles.pointB=[];
handles.pointC=[];

%points for point of interest
handles.apdpoints=[];
handles.apdpointcount=0;
handles.roipointcount=0;
handles.roiareas=zeros(1,1,1);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segEP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = segEP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function APDs_Callback(hObject, eventdata, handles)
% hObject    handle to APDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of APDs as text
%        str2double(get(hObject,'String')) returns contents of APDs as a double


% --- Executes during object creation, after setting all properties.
function APDs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to APDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cvopt.
function cvopt_Callback(hObject, eventdata, handles)
% hObject    handle to cvopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
handles.pointA=[];
handles.pointB=[];
handles.pointC=[];

%reset image
    axes(handles.axes1)
    jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [amap] = handles.amap;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(amap));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(amap, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
% Update handles structure

imshow(g1data.im,[])
hold on
for i=1:size(g1data.boundaries,1) 
plot(g1data.boundaries{i}(:,2),g1data.boundaries{i}(:,1),'r','LineWidth',2);
end

cvopt=get(handles.cvopt,'Value')
if cvopt == 1
    set(handles.winsize,'Enable','on')
    set(handles.binsize,'Enable','off')
    set(handles.veldistance,'Enable','off')
    set(handles.minvel,'Enable','on')
    set(handles.maxvel,'Enable','on')
    set(handles.mint,'Enable','on')
    set(handles.maxt,'Enable','on')
end
if cvopt == 2
    set(handles.winsize,'Enable','on')
    set(handles.binsize,'Enable','on')
    set(handles.veldistance,'Enable','off')
    set(handles.minvel,'Enable','on')
    set(handles.maxvel,'Enable','on')
    set(handles.mint,'Enable','on')
    set(handles.maxt,'Enable','on')
end
if cvopt == 3
    set(handles.winsize,'Enable','on')
    set(handles.binsize,'Enable','on')
    set(handles.veldistance,'Enable','off')
    set(handles.minvel,'Enable','on')
    set(handles.maxvel,'Enable','on')
    set(handles.mint,'Enable','on')
    set(handles.maxt,'Enable','on')
end
if cvopt == 4
    set(handles.winsize,'Enable','off')
    set(handles.binsize,'Enable','off')
    set(handles.veldistance,'Enable','on')
    set(handles.minvel,'Enable','off')
    set(handles.maxvel,'Enable','off')
    set(handles.mint,'Enable','off')
    set(handles.maxt,'Enable','off')
end
if cvopt == 5
    set(handles.winsize,'Enable','off')
    set(handles.binsize,'Enable','off')
    set(handles.veldistance,'Enable','on')
    set(handles.minvel,'Enable','off')
    set(handles.maxvel,'Enable','off')
    set(handles.mint,'Enable','off')
    set(handles.maxt,'Enable','off')
end
if cvopt == 6
    set(handles.winsize,'Enable','off')
    set(handles.binsize,'Enable','on')
    set(handles.veldistance,'Enable','on')
    set(handles.minvel,'Enable','off')
    set(handles.maxvel,'Enable','off')
    set(handles.mint,'Enable','off')
    set(handles.maxt,'Enable','off')
end
if cvopt == 7
    set(handles.winsize,'Enable','off')
    set(handles.binsize,'Enable','off')
    set(handles.veldistance,'Enable','off')
    set(handles.minvel,'Enable','off')
    set(handles.maxvel,'Enable','off')
    set(handles.mint,'Enable','off')
    set(handles.maxt,'Enable','off')
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns cvopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cvopt


% --- Executes during object creation, after setting all properties.
function cvopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cvopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%points for single vector stuff


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function winsize_Callback(hObject, eventdata, handles)
% hObject    handle to winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winsize as text
%        str2double(get(hObject,'String')) returns contents of winsize as a double


% --- Executes during object creation, after setting all properties.
function winsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binsize_Callback(hObject, eventdata, handles)
% hObject    handle to binsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binsize as text
%        str2double(get(hObject,'String')) returns contents of binsize as a double


% --- Executes during object creation, after setting all properties.
function binsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function veldistance_Callback(hObject, eventdata, handles)
% hObject    handle to veldistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of veldistance as text
%        str2double(get(hObject,'String')) returns contents of veldistance as a double


% --- Executes during object creation, after setting all properties.
function veldistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to veldistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minvel_Callback(hObject, eventdata, handles)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minvel as text
%        str2double(get(hObject,'String')) returns contents of minvel as a double


% --- Executes during object creation, after setting all properties.
function minvel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxvel_Callback(hObject, eventdata, handles)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxvel as text
%        str2double(get(hObject,'String')) returns contents of maxvel as a double


% --- Executes during object creation, after setting all properties.
function maxvel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mint_Callback(hObject, eventdata, handles)
% hObject    handle to mint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mint as text
%        str2double(get(hObject,'String')) returns contents of mint as a double


% --- Executes during object creation, after setting all properties.
function mint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxt_Callback(hObject, eventdata, handles)
% hObject    handle to maxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxt as text
%        str2double(get(hObject,'String')) returns contents of maxt as a double


% --- Executes during object creation, after setting all properties.
function maxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function actims_Callback(hObject, eventdata, handles)
% hObject    handle to actims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of actims as text
%        str2double(get(hObject,'String')) returns contents of actims as a double


% --- Executes during object creation, after setting all properties.
function actims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
saveopt=get(handles.saveas,'Value');
if saveopt == 1 || saveopt == 2 || saveopt == 3 || saveopt == 6 || saveopt == 7
[filename,pathname] = uiputfile({'*.csv';'*.mat'},'Hi');
end
if saveopt == 4 || saveopt == 5
[filename,pathname] = uiputfile({'*.mat'},'Hi');
end
% analysis options 
opt1=get(handles.apdopt,'Value');opt2=get(handles.Tauopt,'Value');opt3=0;;opt4=get(handles.ttpopt,'Value');opt5=get(handles.maxuopt,'Value');

[~,~,ext] = fileparts(filename);
file=[pathname,filename];
numsec=length(g1data.section);
%values from segEP GUI
    ts=str2double(get(handles.APDs,'String'));
    apdnames=cell(numel(ts)-1,1);
    tcount=0;
    for k = 1:2:(2*numel(ts))
        tcount=tcount+1;
        apdnames{k}=['APD',num2str(ts(tcount))];
        apdnames{k+1}=['APD',num2str(ts(tcount)),'_Stdev'];
    end
    acttimss=str2double(get(handles.actims,'String'));

%some values from elecromap
numsec=numsec;
Section=[1:numsec]';
Cycle_Length=[g1data.avgCL(2,1:numsec)]';

apdvals=zeros(numsec,length(ts)*2);
vouts=zeros(numsec,1);

%wholetable for angular dists
wholetable=zeros(721,numsec+1);
wholetable(:,1)=(-360:1:360);
binsize=str2double(get(handles.binsize,'String'))
360/binsize;
length(binsize/2:binsize:360-binsize/2);
wholetablemulti=zeros(360/binsize,numsec+1);
wholetablemulti(:,1)=(binsize/2:binsize:360-binsize/2);
cvopt=get(handles.cvopt,'Value')
if cvopt ~= 7
cvonedevs=zeros(numsec,1);
Amax=zeros(numsec,1);
Amaxdiff=zeros(numsec,1);
Conduction_Velocity1=zeros(numsec,1);
Conduction_Velocity2=zeros(numsec,1);
maxp=zeros(numsec,1);
time_at_max=zeros(numsec,1);
time_in_act=zeros(numsec,1);
time_in_repol=zeros(numsec,1);
FWHM=zeros(numsec,1);
ovAUC=zeros(numsec,1);
actAUC=zeros(numsec,1);
maxAUC=zeros(numsec,1);
repolAUC=zeros(numsec,1);
   b2bapmean = zeros(numsec,1);
   b2bapstd =  zeros(numsec,1);
   b2bcvmean = zeros(numsec,1);
   b2bcvstd =  zeros(numsec,1);
else
cvonedevs=cell(numsec,1)
Conduction_Velocity1=cell(numsec,1)
end
%%start going through sections 
wb=waitbar(0,'Producing whole file section analysis: section 1')
APDROI=zeros(numsec,handles.apdpointcount);tauROI=zeros(numsec,handles.apdpointcount);DIROI=zeros(numsec,handles.apdpointcount);ttpROI=zeros(numsec,handles.apdpointcount);maxuROI=zeros(numsec,handles.apdpointcount);
for j=1:numsec
    handles = guidata(hObject);
waitbar(j/numsec,wb,['Producing whole file section analysis: section ',num2str(j)]);
%% Store each peak into array

before=str2double(get(g1data.beforeGUI,'String'));
after=str2double(get(g1data.afterGUI,'String'));
% 
% exposure = time(2); % looks at the time to the first frame this is roughly the exposure
exposure=1/str2double(get(g1data.framerate,'String')); 
before = round(before/exposure) %1000 because we are dealing with ms
 after = round(after/exposure);
section_choice = j;
m=g1data.q2locs(section_choice,:) %peak positons of section
f=m(m~=0);
peaks = (length(f)); % ignores last peak as the signal may cut out 
                            % beyond "after" hence matrix dim error
numpeaks = peaks;
timeframe = 1:before+after+1;

%% OVERLAYING ALL BEATS TO MAKE AN AVERAGE BEAT
% total action potential duration
APtime = before+after;

% create empty matrix to fill later on
overlay = zeros(size(g1data.im,1), size(g1data.im,2), APtime);

% skip the first and last AP to forgo any possible errors exceeding matrix
% dimensions
if f(1) < before
    startloc =2;
else startloc =1
end
locRange = startloc:numel(f);

% fill matrix
% figure('name', 'overlay of all beats'),
if numel(locRange) == 1
if f(locRange)+after > size(g1data.images,3)
   after=size(g1data.images,3)-f(locRange);
end
elseif numel(locRange) > 1
if f(numel(locRange))+after > size(g1data.images,3)
   locRange=locRange(1:end-1);
end    
end
for x = -before:after
%    f
%    locRange
%    x+before+1
%    f(locRange)
   overlay(:,:,x+before+1) = sum(g1data.images(:,:,f(locRange)+x),3)./numel(f);
   overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(g1data.mask);
   %imshow(overlay(:,:,x+before+1),[], 'InitialMagnification', 500);
   % pause(0.01);
end

% title('overlay of all beats');
handles.cvimages=overlay;
%% WRITE TO TIFF STACK
% % normalise
minI = min(overlay(:));
maxI = max(overlay(:));
% 
averageBeat = overlay - minI;
averageBeat = (2^16-1)*averageBeat./(maxI);
% 
%make 16 bit
handles.averageBeat = uint16(averageBeat);

   %get APDs %apdvals( APD1;ST_DEV1;APD2;ST_DEV2 etc.. with rows being sections) 
    for i=1:length(ts)
        [apdmap,apdvals(j,(2*i-1)),~,apdvals(j,2*i)]=mapsbaby(get(g1data.aptime1,'Value'),str2double(get(g1data.framerate,'String')),ts(i),g1data.I,g1data.images,handles.averageBeat,g1data.outlier,str2double(get(g1data.cmin,'String')),str2double(get(g1data.cmax,'String')),get(g1data.tfilt,'Value'),str2double(get(g1data.beforeGUI,'String')),get(g1data.apdbl,'Value'),str2double(get(g1data.apdblnum,'String')));
        APDmap(:,:,j)=apdmap;
    end
    
    
    if saveopt == 2
       for i = 1:handles.apdpointcount
           if j == 1    
           APDROI(j,i)=APDmap(handles.apdpoints(i,2),handles.apdpoints(i,1));
           else
           APDROI(j,i)=APDmap(handles.apdpoints(i,2),handles.apdpoints(i,1),j);
           end
       end
       
    end
     if saveopt == 3
           for i = 1:handles.roipointcount
           roiapds=[];
           anopt=2;
           if anopt == 1
           [apd1,tau1,di1,ttp1,maxu1]=segapd(opt1,opt2,opt3,opt4,opt5,handles.av_sig(i,:),str2double(get(g1data.framerate,'String')),f,str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.afterGUI,'String')),get(g1data.tfilt,'Value'),ts(1),get(g1data.aptime1,'Value'),get(g1data.apdbl,'Value'),str2double(get(g1data.apdblnum,'String')),str2double(get(g1data.taustart,'String')),str2double(get(g1data.taufinish,'String')));    
           APDROI(j,i)=apd1;tauROI(j,i)=tau1;DIROI(j,i)=di1;ttpROI(j,i)=ttp1;maxuROI(j,i)=maxu1;
           end  
           if anopt == 2
           for k=1:size(handles.roiareas,3)
             roiapds(k)=APDmap(handles.roiareas(k,2,i),handles.roiareas(k,1,i),j);
           end    
           if j == 1
           APDROI(j,i)=mean(roiapds);
           else
           APDROI(j,i)=mean(roiapds);
           end
           end
           end
     end

    if saveopt == 6
       [actmap]=activationmapoff(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velalgo,'Value'),str2double(get(g1data.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')),str2double(get(g1data.t,'String')),[]); 
       tim=actmap(actmap>0);
       tim=tim-min(tim);
       allpts=numel(tim);
       xbins=0:0.01:max(tim);
       tissueact=100*cumsum(hist(tim,xbins))/allpts;
       actmaxes=str2double(get(handles.actims,'String'));
       actmin=0;

for k=1:numel(actmaxes)
    actmax=actmaxes(k);
Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if isempty(Imax) == 1 || isempty(Imin) == 1 
Imax=10000
Imin = 10000
end
if actmax < 100
Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);

if isempty(Imax) == 1 || isempty(Imin) == 1 
Imax=10000
Imin = 10000
end

Amax(j,k)=Imax*0.01
end
if k == 2
   Amax(j,k+1)=Amax(j,2)-Amax(j,1);
   if j>1
   Amaxdiff(j)=abs(Amax(j,k+1)-Amax(j-1,k+1));
   else
   Amaxdiff(j)=NaN;
   end
end
    
T=table(Section,Cycle_Length,Amax,Amaxdiff);
    end 
    
    if saveopt == 7
        try
%        [handles.map,~,~,handles.act_t,~,~,~,~,~,~,handles.vout,handles.quivers_Xout,handles.quivers_Yout,handles.quivers_vxout,handles.quivers_vyout,~,~,~,handles.repolmap,handles.repolact_t,~,~]...
%     =cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,1,str2double(get(g1data.minvel,'String')),str2double(get(g1data.maxvel,'String')),get(g1data.velalgo,'Value'),...
%     str2double(get(g1data.MINt,'String')),str2double(get(g1data.MAXt,'String')),str2double(get(g1data.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),1,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String'))); 
       
[actmap,aoff]=activationmapoff(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velalgo,'Value'),str2double(get(g1data.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')),str2double(get(g1data.t,'String')),[]); 
       tim=actmap(actmap>0);
[repolmap]=activationmapoff(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,5,str2double(get(g1data.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')),str2double(get(g1data.t,'String')),aoff); 
       rtim=repolmap(repolmap>0);
tim;
repolmin=min(tim); %save for repol calcs
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;
handles.xbins=xbins;
handles.tissueact=tissueact;
    
repolcurve=1
if repolcurve ==1
%rtim=handles.repolact_t;
rtim=rtim-repolmin;
allpts=numel(rtim);
rxbins=0:0.01:max(rtim);
rtissueact=100*cumsum(hist(rtim,rxbins))/allpts;
handles.rxbins=rxbins;
handles.rtissueact=rtissueact;


% calculate joint activation-repolarisation curve
for x=1:length(rxbins)
    if x>length(xbins)
       tissueact(x) = 100;
    end
end

jcurve=tissueact-rtissueact;
%jcurve stats
% figure,
% plot(jcurve,'x')
%max precetnage
maxp(j)=max(jcurve);
if j == 20
    figure,
    plot(jcurve)
end
%time at maxpercetnage
ipt=(find(jcurve==maxp(j))); %time at maxp
mpt=rxbins(ipt);
time_at_max(j)=mpt(end)-mpt(1);
time_in_act(j)=mpt(1);
time_in_repol(j)=rxbins(end)-mpt(end);

%FWHM
over50t=rxbins(find(jcurve>50));
FWHM(j)=over50t(end)-over50t(1);

%overall AUC

ovAUC(j)=trapz(rxbins,jcurve);

%split AUC
actAUC(j)=trapz(rxbins(1:ipt(1)),jcurve(1:ipt(1)));
maxAUC(j)=trapz(rxbins(ipt(1):ipt(end)),jcurve(ipt(1):ipt(end)));
repolAUC(j)=trapz(rxbins(ipt(end):end),jcurve(ipt(end):end));

T=table(Section,Cycle_Length,maxp,time_at_max,time_in_act,time_in_repol,FWHM,ovAUC,actAUC,maxAUC,repolAUC)

end
        catch
           
T=table(Section,Cycle_Length,maxp,time_at_max,time_in_act,time_in_repol,FWHM,ovAUC,actAUC,maxAUC,repolAUC);

end
    end
    
    
    if saveopt == 1
    %% Multi overall
    if cvopt == 7
       T=table(Section,Cycle_Length);
       for m=1:numel(apdnames)
       aptable=[];
       aptable=table(apdvals(:,m));
       aptable.Properties.VariableNames = {['APD',num2str(m)]};
       T=[T,aptable];
       end
    end
    
    if cvopt == 1 && saveopt ~= 6
        T=table(Section,Cycle_Length);
       for m=1:numel(apdnames)
       aptable=[];
       aptable=table(apdvals(:,m));
       aptable.Properties.VariableNames = {['APD',num2str(m)]};
       T=[T,aptable];
       end
    [actmap,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
    cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(g1data.velalgo,'Value'),...
          str2double(get(handles.mint,'String')),str2double(get(handles.maxt,'String')),str2double(get(handles.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')));
    Conduction_Velocity1(j)=mean(vout);
    CV1=table(Conduction_Velocity1);
    cv_onedev=table(cvonedevs)
    
    %% activation times 
vouts=mean(vout);
tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2double(get(handles.actims,'String'));
actmin=0;

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if isempty(Imax) == 1 || isempty(Imin) == 1 
Imax=10000
Imin = 10000
end
if actmax < 100
Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);

if isempty(Imax) == 1 || isempty(Imin) == 1 
Imax=10000
Imin = 10000
end

Amax(j)=Imax*0.01;
%timdiff=timmax-timmin;
if actmin == 0
 %   timdiff=timmax;
end
CV1
cv_onedev
Amax1=table(Amax)
%calculate B2B diffrences
if j>1
   b2bapmean(j) = abs(apdvals(j,1)-apdvals(j-1,1));
   b2bapstd(j) =  abs(apdvals(j,2)-apdvals(j-1,2));
   b2bcvmean(j) = abs(Conduction_Velocity1(j)-Conduction_Velocity1(j-1));
   b2bcvstd(j) = abs(cvonedevs(j)-cvonedevs(j-1));
   
end
   b2bapmeanT=table(b2bapmean); b2bapstdT=table(b2bapstd); b2bcvmeanT=table(b2bcvmean); b2bcvstdT=table(b2bcvstd);
    T=[T,CV1,cv_onedev,Amax1,b2bapmeanT,b2bapstdT,b2bcvmeanT,b2bcvstdT];
    %writetable(T,file,'Delimiter',',','WriteVariableNames',true);
    end
    
    %% Multi long + tran
    if cvopt == 2
        quiverv=[]
        vangle=[]
    [actmap,~,~,~,~,~,~,~,~,~,handles.vout,handles.quivers_Xout,handles.quivers_Yout,handles.quivers_vxout,handles.quivers_vyout]...
    =cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(g1data.velalgo,'Value'),...
     str2double(get(handles.mint,'String')),str2double(get(handles.maxt,'String')),str2double(get(handles.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String'))); 

    quiverv=[handles.quivers_vxout, handles.quivers_vyout,handles.vout,handles.quivers_Xout,handles.quivers_Yout];
for i =1:length(handles.quivers_vxout)
    vangle(i)=(atan(quiverv(i,1)/quiverv(i,2))*(180/pi));
    if quiverv(i,1) < 0
        vangle(i)=vangle(i)+180;
    end
    vangle(i)=vangle(i)+90;
end
vangle=vangle'
quiverv=[quiverv,vangle];
binsize=str2double(get(handles.binsize,'String'));
hva=hist(vangle,[binsize/2:binsize:360-binsize/2]);
anglevs=zeros(1,length(hva));
for i =1:length(handles.quivers_vxout)
       k=ceil(quiverv(i,6)/binsize);
       anglevs(k)=anglevs(k)+quiverv(i,3);
end
for i=1:length(hva)
    anglevs(i)=anglevs(i)/hva(i);
end
anglevs=anglevs'
min(anglevs)
Conduction_Velocity1(j)=min(anglevs)
Conduction_Velocity2(j)=max(anglevs)
CV_Slow=Conduction_Velocity1
CV_Fast=Conduction_Velocity2
APD=apdvals
    T=table(Section,Cycle_Length,APD,CV_Slow,CV_Fast);
    writetable(T,file,'Delimiter',',','WriteVariableNames',true);
    end
    %% Multi Dist
    if cvopt == 3
        quiverv=[]
        vangle=[]
    [actmap,~,~,~,~,~,~,~,~,~,handles.vout,handles.quivers_Xout,handles.quivers_Yout,handles.quivers_vxout,handles.quivers_vyout]...
    =cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(g1data.velalgo,'Value'),...
     str2double(get(handles.mint,'String')),str2double(get(handles.maxt,'String')),str2double(get(handles.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String'))); 

    quiverv=[handles.quivers_vxout, handles.quivers_vyout,handles.vout,handles.quivers_Xout,handles.quivers_Yout];
for i =1:length(handles.quivers_vxout)
    vangle(i)=(atan(quiverv(i,1)/quiverv(i,2))*(180/pi));
    if quiverv(i,1) < 0
        vangle(i)=vangle(i)+180;
    end
    vangle(i)=vangle(i)+90;
end
vangle=vangle'
quiverv=[quiverv,vangle];
binsize=str2double(get(handles.binsize,'String'));
hva=hist(vangle,[binsize/2:binsize:360-binsize/2]);
anglevs=zeros(1,length(hva));
for i =1:length(handles.quivers_vxout)
       k=ceil(quiverv(i,6)/binsize);
       anglevs(k)=anglevs(k)+quiverv(i,3);
end
for i=1:length(hva)
    anglevs(i)=anglevs(i)/hva(i);
    anglevs(i)
    wholetablemulti(i,j+1)=anglevs(i)
end
T=table(wholetablemulti)
writetable(T,file,'Delimiter',',','WriteVariableNames',false);
    end
    %% Single Manual
if cvopt == 4
    [map,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
    cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(g1data.velalgo,'Value'),...
          str2double(get(handles.mint,'String')),str2double(get(handles.maxt,'String')),str2double(get(handles.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')));

time_A=map(handles.pointA(2),handles.pointA(1));
time_B=map(handles.pointB(2),handles.pointB(1)); %times in ms
dist1=sqrt((handles.pointA(2)-handles.pointB(2))^2+(handles.pointA(1)-handles.pointB(1))^2)*str2double(get(g1data.pixelsize,'String')) %distance in um
speed1=dist1*0.0001/(time_B-time_A)*1000; %speed in cm/s
speed1=sqrt(speed1*speed1)
time1=sqrt((time_B-time_A)^2) 

Conduction_Velocity1(j)=speed1
time_C=map(handles.pointC(2),handles.pointC(1)); %times in ms
dist2=sqrt((handles.pointA(2)-handles.pointC(2))^2+(handles.pointC(1)-handles.pointC(1))^2)*str2double(get(g1data.pixelsize,'String')) %distance in um
speed2=dist2*0.0001/(time_C-time_A)*1000; %speed in cm/s
speed2=sqrt(speed2*speed2);
time2=sqrt((time_C-time_A)^2); 
Conduction_Velocity2(j)=speed2;
APD=apdvals
actmap=map
    T=table(Section,Cycle_Length,APD,Conduction_Velocity1,Conduction_Velocity2);
    writetable(T,file,'Delimiter',',','WriteVariableNames',true);
end

%% single auto
if cvopt == 5
    veldistance=str2double(get(handles.veldistance,'String'));
veldistancepix=round(veldistance*10000/(str2double(get(g1data.pixelsize,'String'))));

rrow=handles.rrow;
rcol=handles.rcol;
[map,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
    cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(g1data.velalgo,'Value'),...
          str2double(get(handles.mint,'String')),str2double(get(handles.maxt,'String')),str2double(get(handles.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')));
th = 0:pi/180:2*pi;
xunit = veldistancepix * cos(th) + rcol;
yunit = veldistancepix * sin(th) + rrow;
xco=100000000000;
yco=100000000000;
count=0;
mat=[];
[rows cols] = size(map);
for i =1:numel(xunit)

    time_A=(map(rrow,rcol));
    if  round(yunit(i)) <= rows &&  round(xunit(i)) <= cols && round(yunit(i)) > 0 &&  round(xunit(i)) > 0 
    time_B=(map(round(yunit(i)),round(xunit(i)))); %y-x switch due to image vs axis difference in matlab
    else time_B =(time_A -1) 
    end
    p1=[rcol,rrow];
    p2=[round(xunit(i)),round(yunit(i))];
    tim=time_B-time_A;
    if tim > 0 && (p1(1)-p2(1)) ~= 0 && (p1(2)-p2(2)) ~=0;
        xcn=p1(1)-p2(1)
        ycn=p1(2)-p2(2)
        %if xcn ~= xco 
        %    if ycn ~= yco %if statment to stop repeating points with slightly diffrenct angles from th(i)
        dc=sqrt((xcn*xcn)+(ycn*ycn))
        count=count+1;
        mat(count,1)=p1(1)-p2(1); 
        mat(count,2)=p1(2)-p2(2);
        mat(count,3)=tim; %time diffrence 
        mat(count,4)=(dc*0.0001*(str2double(get(g1data.pixelsize,'String'))))/tim*1000; %speed in cm/s
        mat(count,5)=th(i)*(180/pi);
        mat(count,6)=dc; %pixel distance
        xco=xcn;
        yco=ycn;
         %   end
        %end
    end
end

mat=unique(mat,'rows');%get rid of repeat points
outs=get(g1data.velout,'Value');
newcount=1;
outs=2 %ADD IN OUTLIER REMOVAL FOR THIS AND OTHER MESAURES
if outs == 2
    singlevels=deleteoutliers(mat(:,4));
    velmax=max(singlevels);
    velmin=min(singlevels);
    for i=1:length(mat(:,4))
        if mat(i,4) <= velmax && mat(i,4) >= velmin
            mat2(newcount,1)=mat(i,1);
            mat2(newcount,2)=mat(i,2);
            mat2(newcount,3)=mat(i,3);
            mat2(newcount,4)=mat(i,4);
            mat2(newcount,5)=mat(i,5);
            mat2(newcount,6)=mat(i,6);
            newcount=newcount+1;
        end
    end
    mat=[];
    mat=mat2;
end

if outs == 3 || outs == 4 || outs == 5 || outs == 6 || outs == 7 
   onedev=std(mat(:,4));
   velmax=mean(mat(:,4))+(outs-2)*onedev;
   velmin=mean(mat(:,4))-(outs-2)*onedev;
    for i=1:length(mat(:,4))
        if mat(i,4) < velmax && mat(i,4) > velmin
            mat2(newcount,1)=mat(i,1);
            mat2(newcount,2)=mat(i,2);
            mat2(newcount,3)=mat(i,3);
            mat2(newcount,4)=mat(i,4);
            mat2(newcount,5)=mat(i,5);
            mat2(newcount,6)=mat(i,6);
            newcount=newcount+1;
        end
    end
    mat=[];
    mat=mat2;
end

[vel_quick, ind_quick] = max(mat(:,4));
[vel_slow, ind_slow] = min(mat(:,4));

Conduction_Velocity1(j)=mat(ind_slow,4);
Conduction_Velocity2(j)=mat(ind_quick,4);
CV_Slow=Conduction_Velocity1;
CV_Fast=Conduction_Velocity2;
APD=apdvals;
    T=table(Section,Cycle_Length,APD,CV_Slow,CV_Fast);
    writetable(T,file,'Delimiter',',','WriteVariableNames',true);
actmap=map


end
   %% Single Dist

if cvopt == 6
veldistance=str2double(get(handles.veldistance,'String'));
veldistancepix=round(veldistance*10000/(str2double(get(g1data.pixelsize,'String'))));

rrow=handles.rrow;
rcol=handles.rcol;
[map,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
    cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),handles.cvimages,g1data.mask,get(g1data.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(g1data.velalgo,'Value'),...
          str2double(get(handles.mint,'String')),str2double(get(handles.maxt,'String')),str2double(get(handles.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),0,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String')));
th = 0:pi/180:2*pi;
xunit = veldistancepix * cos(th) + rcol;
yunit = veldistancepix * sin(th) + rrow;
xco=100000000000;
yco=100000000000;
count=0;
mat=[];
[rows cols] = size(map);
for i =1:numel(xunit)

    time_A=(map(rrow,rcol));
    if  round(yunit(i)) <= rows &&  round(xunit(i)) <= cols && round(yunit(i)) > 0 &&  round(xunit(i)) > 0 
    time_B=(map(round(yunit(i)),round(xunit(i)))); %y-x switch due to image vs axis difference in matlab
    else time_B =(time_A -1) 
    end
    p1=[rcol,rrow];
    p2=[round(xunit(i)),round(yunit(i))];
    tim=time_B-time_A;
    if tim > 0 && (p1(1)-p2(1)) ~= 0 && (p1(2)-p2(2)) ~=0;
        xcn=p1(1)-p2(1)
        ycn=p1(2)-p2(2)
        %if xcn ~= xco 
        %    if ycn ~= yco %if statment to stop repeating points with slightly diffrenct angles from th(i)
        dc=sqrt((xcn*xcn)+(ycn*ycn))
        count=count+1;
        mat(count,1)=p1(1)-p2(1); 
        mat(count,2)=p1(2)-p2(2);
        mat(count,3)=tim; %time diffrence 
        mat(count,4)=(dc*0.0001*(str2double(get(g1data.pixelsize,'String'))))/tim*1000; %speed in cm/s
        mat(count,5)=th(i)*(180/pi);
        mat(count,6)=dc; %pixel distance
        xco=xcn;
        yco=ycn;
         %   end
        %end
    end
end

mat=unique(mat,'rows');%get rid of repeat points
[vel_slow, ind_slow] = min(mat(:,4));

zeroangle=[mat(ind_slow,5)];
vangle(:,2)=mat(:,4);
vangle(:,1)=(mat(:,5)-zeroangle);
mat

vangle = sortrows(vangle,1);
for i=1:length(vangle(:,1))-1
    if(vangle(i+1,1)) == vangle (i,1)
        vangle(i,1)=375 %get rid of repeat angles
    end
end
dcount=1
angind=[];
for d=1:1:721
angind=find(vangle(:,1) == (d-361))
if isempty(angind) == 0;
    wholetable(d,j+1)=vangle(angind,2);
    dcount=dcount+1;
else wholetable(d,j+1)=NaN;
end
angind=[];
end
vangle;
vangle=[];
        T=table(wholetable);
writetable(T,file,'Delimiter',',','WriteVariableNames',false);
actmap=map;
end
    end
end

if saveopt == 1
             %% Modify variable names and save table
        varnames=[];
       if cvopt == 1 || cvopt == 7
       if cvopt == 7
       varnames=cell(numel(apdnames)+2,1)
       varnames{1}='Section';varnames{2}='Cycle_Length';
       end
       if cvopt==1
       varnames=cell(numel(apdnames)+9,1)
       varnames{1}='Section';varnames{2}='Cycle_Length';varnames{numel(apdnames)+3}='CV';varnames{numel(apdnames)+4}='CV_Stdev';varnames{numel(apdnames)+5}='Act_Time';
       end
       for i=3:numel(apdnames)+2;
           varnames{i}=apdnames{i-2};
       end
       varnames{8}=['B2BAPD'];varnames{9}=['B2BAPDSTD'];varnames{10}=['B2BCV'];varnames{11}=['B2BCVSTD'];
       T
       varnames
       %T.Properties.VariableNames = varnames;
       writetable(T,file,'Delimiter',',','WriteVariableNames',true)
       end
end
%% Save APDs from points/rois
   if saveopt == 2 || saveopt == 3
     [a,b,ext] = fileparts(filename);
     if opt1 == 1 
     apdT=table(APDROI);
     apdfilename=[b,'_apd',ext];
     apdfile=[pathname,apdfilename];
     writetable(apdT,apdfile,'Delimiter',',','WriteVariableNames',true)
     end
     
     if opt2 == 1 
     tauT=table(tauROI);
     taufilename=[b,'_tau',ext];
     taufile=[pathname,taufilename];
     writetable(tauT,taufile,'Delimiter',',','WriteVariableNames',true)
     end
     
     if opt3 == 1 
     DIT=table(DIROI);
     DIfilename=[b,'_DI',ext];
     DIfile=[pathname,DIfilename];
     writetable(DIT,DIfile,'Delimiter',',','WriteVariableNames',true)
     end
     
     if opt4 == 1 
     ttpT=table(ttpROI);
     ttpfilename=[b,'_ttp',ext];
     ttpfile=[pathname,ttpfilename];
     writetable(ttpT,ttpfile,'Delimiter',',','WriteVariableNames',true)
     end
     
     if opt5 == 1 
     maxuT=table(maxuROI);
     maxufilename=[b,'_maximum_upstroke_velocity',ext];
     maxufile=[pathname,maxufilename];
     writetable(maxuT,maxufile,'Delimiter',',','WriteVariableNames',true)
     end
     
     if get(handles.savesigs,'Value') == 1
        ROIsignal=handles.av_sig';
        for s=1:length(ROIsignal(1,:))
           if get(g1data.invertopt,'Value') == 1
               ROIsignal(:,s)=imcomplement(ROIsignal(:,s))
           end
            ROIsignal(:,s)= ROIsignal(:,s)-min( ROIsignal(:,s));
            if get(handles.normsig,'Value') == 1
            ROIsignal(:,s)= ROIsignal(:,s)./max( ROIsignal(:,s));
            end
        end
        time=[1:length(ROIsignal(:,1))].*(1/str2double(get(g1data.framerate,'String')));
        time=time';
        sigfilename=[b,'_signal',ext];
        sigfile=[pathname,sigfilename];
        sigT=table(time,ROIsignal);
        writetable(sigT,sigfile,'Delimiter',',','WriteVariableNames',true)
     end
   end

%% Save APD maps
if saveopt == 4
        save(file, 'APDmap')
end
    
if saveopt == 6
    writetable(T,file,'Delimiter',',','WriteVariableNames',true)
end
if saveopt == 7
    writetable(T,file,'Delimiter',',','WriteVariableNames',true)
end
delete(wb)



function getMousePositionOnImage(hObject, eventdata, handles)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
cvopt=get(handles.cvopt,'Value');
saveopt=get(handles.saveas,'Value');

if saveopt == 1
if cvopt == 5 || cvopt == 6
    
%reset image
    axes(handles.axes1)
    cla
    jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [amap] = handles.amap;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(amap));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(amap, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
imshow(g1data.im,[])

hold on
cursorPoint = get(handles.axes1, 'CurrentPoint');
handles.curX = cursorPoint(1,1);
handles.curY = cursorPoint(1,2);
handles.curX=30;
handles.curY=30;
xLimits = get(handles.axes1, 'xlim');
yLimits = get(handles.axes1, 'ylim');

if (handles.curX > min(xLimits) && handles.curX < max(xLimits) && handles.curY > min(yLimits) && handles.curY < max(yLimits))

veldistance=str2double(get(handles.veldistance,'String'));
veldistancepix=round(veldistance*10000/(str2double(get(g1data.pixelsize,'String'))));
handles.rrow=floor(handles.curY);
handles.rcol=floor(handles.curX);
rcol=handles.rcol;
rrow=handles.rrow;
plot(rcol,rrow,'k+','MarkerSize',20,'LineWidth',3);
th = 0:pi/180:2*pi;
xunit = veldistancepix * cos(th) + rcol;
yunit = veldistancepix * sin(th) + rrow;
h = plot(xunit, yunit,'k','LineWidth',3);
guidata(hObject, handles);
end
end

if cvopt == 4
     if isempty(handles.pointA) == 1
    axes(handles.axes1)
    jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [amap] = handles.amap;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(amap));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(amap, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
     end
imshow(g1data.im,[])
hold on
imshow(g1data.im,[])
cursorPoint = get(handles.axes1, 'CurrentPoint');
handles.curX = cursorPoint(1,1);
handles.curY = cursorPoint(1,2);
%overide to centre of image

xLimits = get(handles.axes1, 'xlim');
yLimits = get(handles.axes1, 'ylim');


A=handles.pointA
B=handles.pointB
C=handles.pointC

if (handles.curX > min(xLimits) && handles.curX < max(xLimits) && handles.curY > min(yLimits) && handles.curY < max(yLimits))     
    if isempty(handles.pointA) == 1
               handles.pointA=[round(handles.curX), round(handles.curY)]
               axes(handles.axes1);
               hold on
               plot(handles.pointA(1),handles.pointA(2),'k+','MarkerSize',20,'LineWidth',3);
           elseif isempty(handles.pointA) == 0 && isempty(handles.pointB) == 1
               handles.pointB=[round(handles.curX), round(handles.curY)]
               dp1=handles.pointB-handles.pointA;
q1=quiver(handles.pointA(1),handles.pointA(2),dp1(1),dp1(2),0);
q1.Color='k';
q1.LineWidth=3;
q1.MaxHeadSize=0.5;
elseif isempty(handles.pointA) == 0 && isempty(handles.pointB) == 0 && isempty(handles.pointC) == 1
               handles.pointC=[round(handles.curX), round(handles.curY)];
               dp2=handles.pointC-handles.pointA;
q2=quiver(handles.pointA(1),handles.pointA(2),dp2(1),dp2(2),0);
q2.Color='r';
q2.LineWidth=3;
q2.MaxHeadSize=0.5;
           end
end
guidata(hObject, handles);
end
end

if saveopt == 2
    cursorPoint = get(handles.axes1, 'CurrentPoint');
handles.curX = cursorPoint(1,1);
handles.curY = cursorPoint(1,2);
xLimits = get(handles.axes1, 'xlim');
yLimits = get(handles.axes1, 'ylim');
      set(handles.removepoints,'Visible','on')
    set(handles.removerois,'Visible','off')
    axes(handles.axes1)
    imshow(g1data.im,[])
    for i=1:size(g1data.boundaries,1) 
plot(g1data.boundaries{i}(:,2),g1data.boundaries{i}(:,1),'r','LineWidth',2);
end
     if (handles.curX > min(xLimits) && handles.curX < max(xLimits) && handles.curY > min(yLimits) && handles.curY < max(yLimits))     
        handles.apdpointcount=handles.apdpointcount+1
        handles.apdpoints(handles.apdpointcount,1)=round(handles.curX);
        handles.apdpoints(handles.apdpointcount,2)=round(handles.curY);
        si=[];
        for i = 1:handles.apdpointcount
        a=jet(handles.apdpointcount)
        if i == 1
            a(i,:)=[0 0 1]
        end
        if i == handles.apdpointcount
            a(i,:)=[1 0 0]
        end
        plot(handles.apdpoints(i,1),handles.apdpoints(i,2),'+','MarkerSize',20,'LineWidth',3,'Color',a(i,:));
        if get(handles.roilabel,'Value') == 1
           text(handles.apdpoints(i,1)+1,handles.apdpoints(i,2)-2,['point ',num2str(i)],'Color',(a(i,:)))
        end
        si(i,:) = (squeeze((g1data.images(handles.apdpoints(i,2),handles.apdpoints(i,1),:))));
        end
        handles.av_sig=si;
     end
end

if saveopt == 3
    cursorPoint = get(handles.axes1, 'CurrentPoint');
handles.curX = cursorPoint(1,1);
handles.curY = cursorPoint(1,2);
xLimits = get(handles.axes1, 'xlim');
yLimits = get(handles.axes1, 'ylim');
shapeopt=get(handles.roishape,'Value');
sizeopt=str2double(get(handles.roisize,'String'));
    set(handles.removepoints,'Visible','off')
    set(handles.removerois,'Visible','on')
    axes(handles.axes1)
    imshow(g1data.im,[])
    for i=1:size(g1data.boundaries,1) 
    plot(g1data.boundaries{i}(:,2),g1data.boundaries{i}(:,1),'r','LineWidth',2);
    end
     if (handles.curX > min(xLimits) && handles.curX < max(xLimits) && handles.curY > min(yLimits) && handles.curY < max(yLimits))     
        handles.roipointcount=handles.roipointcount+1
        handles.roipoints(handles.roipointcount,1)=round(handles.curX);
        handles.roipoints(handles.roipointcount,2)=round(handles.curY);
        for i = 1:handles.roipointcount
        a=jet(handles.roipointcount)
        if i == 1
            a(i,:)=[0 0 1]
        end
        if i == handles.roipointcount
            a(i,:)=[1 0 0]
        end
        if shapeopt == 1
          pos_vec=[handles.roipoints(i,1)-((sizeopt-1)/2),handles.roipoints(i,2)-((sizeopt-1)/2),sizeopt,sizeopt]
          rectangle('Position',pos_vec,'EdgeColor',a(i,:),'LineWidth',3);
          if get(handles.roilabel,'Value') == 1
             text(pos_vec(1)+sizeopt+1,pos_vec(2),['ROI',num2str(i)],'Color',(a(i,:)))
          end
          pcount = 0;
          for r =1:sizeopt
              for c=1:sizeopt
                  pcount=pcount+1;
                  handles.roiareas(pcount,1,i)=handles.roipoints(i,1)+(r-((sizeopt-1)/2));
                  handles.roiareas(pcount,2,i)=handles.roipoints(i,2)+(c-((sizeopt-1)/2));
                  si(pcount,:) = (squeeze((g1data.images(handles.roiareas(pcount,2,i),handles.roiareas(pcount,1,i),:))));
              end
          end
        end
        handles.av_sig(i,:)=sum(si);
     end
     end
     
end
guidata(hObject, handles);
         % --- Executes on selection change in saveas.
function saveas_Callback(hObject, eventdata, handles)
% hObject    handle to saveas (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
cvopt=get(handles.cvopt,'Value');
saveopt=get(handles.saveas,'Value');

if saveopt == 2 || saveopt == 3 || saveopt == 4
    set(handles.winsize,'Enable','off')
    set(handles.binsize,'Enable','off')
    set(handles.veldistance,'Enable','off')
    set(handles.minvel,'Enable','off')
    set(handles.maxvel,'Enable','off')
    set(handles.mint,'Enable','off')
    set(handles.maxt,'Enable','off')
end

if saveopt == 3
   
end


% Hints: contents = cellstr(get(hObject,'String')) returns saveas contents as cell array
%        contents{get(hObject,'Value')} returns selected item from saveas


% --- Executes during object creation, after setting all properties.
function saveas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in removerois.
function removerois_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    h = findobj('Tag','ElectroMap');
g1data = guidata(h);
    handles.roipoints=[];
    handles.roipointcount=0;
    handles.roiareas=[];
    handles.av_sig=[];
    pcount=0;
    axes(handles.axes1)
    imshow(g1data.im,[])
        for i=1:size(g1data.boundaries,1) 
plot(g1data.boundaries{i}(:,2),g1data.boundaries{i}(:,1),'r','LineWidth',2);
end
    guidata(hObject, handles);
% hObject    handle to removerois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in removepoints.
function removepoints_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    h = findobj('Tag','ElectroMap');
g1data = guidata(h);
    handles.apdpoints=[]
    handles.apdpointcount=0;
    handles.av_sig=[];
    axes(handles.axes1)
    imshow(g1data.im,[])
        for i=1:size(g1data.boundaries,1) 
plot(g1data.boundaries{i}(:,2),g1data.boundaries{i}(:,1),'r','LineWidth',2);
end
    guidata(hObject, handles);
% hObject    handle to removepoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in roishape.
function roishape_Callback(hObject, eventdata, handles)
% hObject    handle to roishape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roishape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roishape


% --- Executes during object creation, after setting all properties.
function roishape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roishape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roisize_Callback(hObject, eventdata, handles)
% hObject    handle to roisize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roisize as text
%        str2double(get(hObject,'String')) returns contents of roisize as a double


% --- Executes during object creation, after setting all properties.
function roisize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roisize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in roiana.
function roiana_Callback(hObject, eventdata, handles)
% hObject    handle to roiana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiana contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiana


% --- Executes during object creation, after setting all properties.
function roiana_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savesigs.
function savesigs_Callback(hObject, eventdata, handles)
% hObject    handle to savesigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savesigs


% --- Executes on button press in normsig.
function normsig_Callback(hObject, eventdata, handles)
% hObject    handle to normsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normsig


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
axes(handles.axes1)
GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;ax=axes;clf;
new_handle=copyobj(handles.axes1,fig);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in imdisplay.
function imdisplay_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns imdisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imdisplay


% --- Executes during object creation, after setting all properties.
function imdisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roilabel.
function roilabel_Callback(hObject, eventdata, handles)
% hObject    handle to roilabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roilabel


% --- Executes on button press in apdopt.
function apdopt_Callback(hObject, eventdata, handles)
% hObject    handle to apdopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of apdopt


% --- Executes on button press in Tauopt.
function Tauopt_Callback(hObject, eventdata, handles)
% hObject    handle to Tauopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tauopt


% --- Executes on button press in DIopt.
function DIopt_Callback(hObject, eventdata, handles)
% hObject    handle to DIopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DIopt


% --- Executes on button press in ttpopt.
function ttpopt_Callback(hObject, eventdata, handles)
% hObject    handle to ttpopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ttpopt


% --- Executes on button press in maxuopt.
function maxuopt_Callback(hObject, eventdata, handles)
% hObject    handle to maxuopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxuopt
