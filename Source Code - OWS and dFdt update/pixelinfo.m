function varargout = pixelinfo(varargin)
% function for compare pixels GUI (Name issue will be rectified in future release)
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary
% PIXELINFO MATLAB code for pixelinfo.fig
%      PIXELINFO, by itself, creates a new PIXELINFO or raises the existing
%      singleton*.
%
%      H = PIXELINFO returns the handle to a new PIXELINFO or the handle to
%      the existing singleton*.
%
%      PIXELINFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIXELINFO.M with the given input arguments.
%
%      PIXELINFO('Property','Value',...) creates a new PIXELINFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pixelinfo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pixelinfo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pixelinfo

% Last Modified by GUIDE v2.5 29-Jan-2018 14:10:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pixelinfo_OpeningFcn, ...
                   'gui_OutputFcn',  @pixelinfo_OutputFcn, ...
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


% --- Executes just before pixelinfo is made visible.
function pixelinfo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pixelinfo (see VARARGIN)

% Choose default command line output for pixelinfo
handles.output = hObject
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
handles.t=str2num(get(g1data.t,'String'));
axes(handles.image);
im=g1data.I;
imshow(im);
wb=waitbar(0.5,'Transfering Maps');
%maps 
handles.map=g1data.apdmap;
handles.alll=g1data.apalll;

handles.isomap=g1data.actmap;
handles.vout=g1data.vout;
handles.quivers_Xout=g1data.quivers_Xout;
handles.quivers_Yout=g1data.quivers_Yout;
handles.quivers_vxout=g1data.quivers_vxout;
handles.quivers_vyout=g1data.quivers_vyout;

delete(wb)
handles.points=0;
handles.rows=[];
handles.cols=[];
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pixelinfo wait for user response (see UIRESUME)
% uiwait(handles.pixelinfo);


% --- Outputs from this function are returned to the command line.
function varargout = pixelinfo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[filename,pathname] = uiputfile({'*.txt';'*.mat'}, 'Save signal');
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
handles.tictoc
T=table(handles.tictoc'./1000,handles.zeropixelsignal');
writetable(T,file,'Delimiter',' ','WriteVariableNames',false);

function getMousePositionOnImage(hObject, eventdata, handles)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
handles.cm=['b','r','g','m','c','k'];
g1data = guidata(h);
cursorPoint = get(handles.image, 'CurrentPoint');
curX = cursorPoint(1,1);
curY = cursorPoint(1,2);
row=floor(curY);
col=floor(curX);
xLimits = get(handles.image, 'xlim');
yLimits = get(handles.image, 'ylim');

if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits))
   handles.points=handles.points+1;
handles.rows(handles.points)=row;
handles.cols(handles.points)=col;
    [~,~,num] = size(g1data.images(:,:,:));
    num;
for i=1:num  %original signal
     pixelsignal(i)=g1data.images(row,col,i);
     presignal(i)=g1data.preimages(row,col,i);
end
exposure=(1/str2num(get(g1data.framerate,'String')));
time=([1:num].*exposure);
inversion=get(g1data.invertopt, 'Value');
    pixelsignal=imcomplement(pixelsignal);


zeropre=imcomplement(presignal);
mini2=min(zeropre);
zeropre=zeropre-mini2;

%temporal filter
tfilt=get(g1data.tfilt, 'Value')
if tfilt == 2
    pixelsignal=double(pixelsignal);
    pixelsignal=sgolayfilt(pixelsignal, 3,11);
end
if tfilt == 3
   d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
   pixelsignal=filtfilt(d,pixelsignal);
end
handles.tictoc=time;

%image
axes(handles.image)
hold on
plot(col,row,'+','MarkerSize',20,'LineWidth',5,'Color',handles.cm(handles.points));


%signal
mini=min(pixelsignal);
zeropixelsignal=pixelsignal-mini;
zeropre=zeropre*(max(zeropixelsignal)/max(zeropre));
handles.zeropre=zeropre;
handles.zeropixelsignal=zeropixelsignal;
axes(handles.signal);
hold on
plot(handles.tictoc,handles.zeropixelsignal,handles.cm(handles.points),'LineWidth',2)
xlabel('Time (ms)');
ylabel('Fluorescence Intensity (arb units)');
axis tight 
%APD
axes(handles.APD)
mapsbabyonepix((get(g1data.aptime1,'Value')),str2num(get(g1data.framerate,'String')),30,g1data.mask,g1data.images,g1data.averageBeat,row,col,handles.cm(handles.points),str2num(get(g1data.beforeGUI,'String')),str2num(get(g1data.afterGUI,'String')),get(g1data.apdbl,'Value'),str2num(get(g1data.apdblnum,'String')),str2num(get(g1data.taustart,'String')),str2num(get(g1data.taufinish,'String')),get(handles.normalise,'Value'),get(g1data.tfilt,'Value'));
xlabel('Time (ms)')
ylabel('Fluorescence Intensity (arb units)');
axis tight

else
disp('Cursor is outside bounds of image.');
end

guidata(hObject, handles);


% --- Executes on button press in rawon.
function rawon_Callback(hObject, eventdata, handles)
% hObject    handle to rawon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hint: get(hObject,'Value') returns toggle state of rawon


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.points = 0;
handles.rows=[];
handles.cols=[];
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
axes(handles.image); cla; 
guidata(hObject, handles);
imchoice_Callback(hObject, eventdata, handles)
axes(handles.signal); cla;
axes(handles.APD); cla;
guidata(hObject, handles);


% --- Executes on selection change in imchoice.
function imchoice_Callback(hObject, eventdata, handles)
% hObject    handle to imchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
axes(handles.image)
cla
im=g1data.I;
% Image
if get(handles.imchoice,'Value') == 1
imshow(im);
end

if get(handles.imchoice,'Value') == 2
    jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    map=handles.map;
    map(isnan(map)) = 0;
    imshow(map,[], 'InitialMagnification', 800,'colormap',jetcolormap);
title(['APD', num2str(handles.t), ' Distribution']);
apscale=get(g1data.apdscale,'Value');

if apscale == 1
   caxis([min(handles.alll) max(handles.alll)]) 
end
if apscale == 2
caxis([str2num(get(g1data.cmin,'String')) str2num(get(g1data.cmax,'String'))]);
end
end

if get(handles.imchoice,'Value') == 3
      jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [map] = handles.isomap;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(map));
elseif get(g1data.isoopt,'Value') == 2
mini=str2num(get(g1data.isomin,'String'));
maxi=str2num(get(g1data.isomax,'String'));
end
imshow(map, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
end

if get(handles.imchoice,'Value') == 4
          jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [map] = handles.isomap;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(map));
elseif get(g1data.isoopt,'Value') == 2
mini=str2num(get(g1data.isomin,'String'));
maxi=str2num(get(g1data.isomax,'String'));
end
imshow(map, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);

scal = 0.5;
hold on
quiver(handles.quivers_Xout,handles.quivers_Yout,scal*handles.quivers_vxout,scal*handles.quivers_vyout,0,'k');
hold off
end

%Points
isempty(handles.rows)
if isempty(handles.rows) == 0
   hold on
   for i=1:length(handles.rows)
plot(handles.cols(i),handles.rows(i),'+','MarkerSize',20,'LineWidth',5,'Color',handles.cm(i));
   end
hold off
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns imchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imchoice


% --- Executes during object creation, after setting all properties.
function imchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in avsave.
function avsave_Callback(hObject, eventdata, handles)
% hObject    handle to avsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in fullsave.
function fullsave_Callback(hObject, eventdata, handles)
% hObject    handle to fullsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in normalise.
function normalise_Callback(hObject, eventdata, handles)
% hObject    handle to normalise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalise
