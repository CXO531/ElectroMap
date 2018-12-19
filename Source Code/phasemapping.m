function varargout = phasemapping(varargin)
% Main function for running phase mapping GUI
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, Please see 'licsence.txt' at ...
 
% Last Updated -
 
% Update Summary

% PHASEMAPPING MATLAB code for phasemapping.fig
%      PHASEMAPPING, by itself, creates a new PHASEMAPPING or raises the existing
%      singleton*.
%
%      H = PHASEMAPPING returns the handle to a new PHASEMAPPING or the handle to
%      the existing singleton*.
%
%      PHASEMAPPING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASEMAPPING.M with the given input arguments.
%
%      PHASEMAPPING('Property','Value',...) creates a new PHASEMAPPING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phasemapping_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phasemapping_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phasemapping

% Last Modified by GUIDE v2.5 08-Jun-2017 14:56:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phasemapping_OpeningFcn, ...
                   'gui_OutputFcn',  @phasemapping_OutputFcn, ...
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


% --- Executes just before phasemapping is made visible.
function phasemapping_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phasemapping (see VARARGIN)
handles.output = hObject;
h = findobj('Tag','ElectroMap');
wb=waitbar(0.5, 'Performing Phase Calculations');
%% get data from ElectroMap
g1data = guidata(h);
handles.row=[];
handles.col=[];
handles.pointcount=0;
handles.images=g1data.images;
handles.mask=g1data.mask;
handles.cm=['b','r','g','m','y'];
handles.exposure=1/str2num(get(g1data.framerate,'String'));
[handles.rows, handles.cols, handles.num]= size(handles.images(:,:,:));
section_choice=get(g1data.listbox2,'Value');
handles.CL=round(g1data.avgCL(2,section_choice),-1); 
B = 1/handles.CL*ones(handles.CL,1);

%% compute phase values using hilbert transform 
for row = 1:handles.rows
    for col =1:handles.cols
        pixelsignal=handles.images(row,col,1:handles.num); % get pixel signal
        out = filter(B,1,pixelsignal);    
        pixelsignal=double(pixelsignal);
        shiftedsignal=pixelsignal-out; %average filter signal 
        shiftedpixelsignal=shiftedsignal(handles.CL:length(shiftedsignal)); %remove first and last sections < CL used for shift
        hsig=hilbert(shiftedpixelsignal); % perform hilbert transfrom  
        handles.phases(row,col,:)=-1*angle(hsig); %find phase angle for each pixel in each frame
        handles.freq(row,col,:)=gradient((-1*angle(hsig)));
        handles.xs(row,col,:)=real(hsig);
        handles.ys(row,col,:)=imag(hsig); %find components of phase angle
    end
end

axes(handles.Videoaxes)
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
A=handles.phases(:,:,1);
A(handles.mask == 0) =NaN;
ddd=[0 0 0;jet];
imshow(A, [-pi pi],'ColorMap',ddd,'InitialMagnification', 800);
    
% Choose default command line output for phasemapping
handles.output = hObject;
delete(wb)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes phasemapping wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = phasemapping_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.Videoaxes)
handles = guidata(hObject);
i=get(handles.slider1,'Value');
i=i*(handles.num-handles.CL);
i=floor(i);
if i ==0
    i=1;
end
A=handles.phases(:,:,i);
    handles.mask(handles.mask == 0) =NaN;
    ddd=[0 0 0;jet];
imshow(A, [-pi pi],'ColorMap',ddd,'InitialMagnification', 800);

if isempty(handles.row)==0
    axes(handles.phasedelay); cla;
    axes(handles.phase); cla;
    axes(handles.Videoaxes)
    for k=1:length(handles.row)
        hold on
        axes(handles.Videoaxes)
         plot(handles.col(k),handles.row(k),'+','MarkerSize',20,'LineWidth',5,'Color',handles.cm(k));

    x=squeeze(handles.xs(handles.row(k),handles.col(k),:));
    y=squeeze(handles.ys(handles.row(k),handles.col(k),:));
    a=squeeze(handles.phases(handles.row(k),handles.col(k),:));
    freq=squeeze(handles.freq(handles.row(k),handles.col(k),:));
    a=sgolayfilt(a,3,11)
    freq=sgolayfilt(freq,3,11)
    axes(handles.phasedelay);
    hold on
    plot(x(1:i),y(1:i),'+','Color',handles.cm(k))
    axes(handles.phase)
    plot(1:length(a),a,'Color',handles.cm(k))
    end
end
axes(handles.phase)
hold on
ax=gca;
hx = plot([i i], ylim, 'Color', 'k');
yticks([-pi 0 pi ])
yticklabels({'-\pi','0','\pi'})
%pause(0.1)
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
answer = questdlg('Warning! Video can not be pasued currently, would you still like to play?', ...
	'Video Warning', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
   for i=1:handles.num-handles.CL
    handles.frame=i;
    set(handles.slider1,'Value',handles.frame/(handles.num-handles.CL));
    slider1_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
end
    case 'No'
end

guidata(hObject, handles);

% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
handles = guidata(hObject);
    handles.row=[];
    handles.col=[];
    handles.pointcount=0;
    axes(handles.phase); cla;
    axes(handles.phasedelay);cla;
    guidata(hObject, handles);
    slider1_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in phasecalc.
function phasecalc_Callback(hObject, eventdata, handles)
% hObject    handle to phasecalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns phasecalc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from phasecalc


% --- Executes during object creation, after setting all properties.
function phasecalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phasecalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in rotorcalc.
function rotorcalc_Callback(hObject, eventdata, handles)
% hObject    handle to rotorcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rotorcalc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rotorcalc


% --- Executes during object creation, after setting all properties.
function rotorcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotorcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function getMousePositionOnImage(hObject, eventdata, handles)
handles = guidata(hObject);
disp('hi')
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
cursorPoint = get(handles.Videoaxes, 'CurrentPoint');
curX = cursorPoint(1,1);
curY = cursorPoint(1,2);
rrow=floor(curY);
rcol=floor(curX);
xLimits = get(handles.Videoaxes, 'xlim');
yLimits = get(handles.Videoaxes, 'ylim');
[rows, cols, num] = size(handles.images(:,:,:))
if rows == 1
        rrow = 1 
        curY= min(yLimits)+((max(yLimits)-min(yLimits))/2);
    end
    if cols == 1
        rcol = 1 
        curX= min(xLimits)+((max(xLimits)-min(xLimits))/2);
    end 
    
if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits)) && handles.mask(rrow,rcol)~=0
    disp('hiiiii')
    handles.pointcount=handles.pointcount+1;
    if handles.pointcount == 6
        handles.pointcount = 1;
        
    end
    handles.row(handles.pointcount)=rrow;
    handles.col(handles.pointcount)=rcol;
    
    axes(handles.Videoaxes)
    guidata(hObject, handles);
    slider1_Callback(hObject, eventdata, handles)
    
    
i=get(handles.slider1,'Value');
i=i*(handles.num-handles.CL);
i=floor(i);

x=squeeze(handles.xs(handles.row(handles.pointcount),handles.col(handles.pointcount),:));
    y=squeeze(handles.ys(handles.row(handles.pointcount),handles.col(handles.pointcount),:));
    a=squeeze(handles.phases(handles.row(handles.pointcount),handles.col(handles.pointcount),:));
     a=sgolayfilt(a,3,11)
    axes(handles.phasedelay);
     plot(x(1:i),y(1:i),'+','Color',handles.cm(handles.pointcount))
     axes(handles.phase);
     hold on
     plot(1:length(a),a,'Color',handles.cm(handles.pointcount))


yticks([-pi 0 pi ])
yticklabels({'-\pi','0','\pi'})
xlabel('Time(s)')
ylabel('Phase')
ax=gca;
yticks([-pi 0 pi ])
yticklabels({'-\pi','0','\pi'})
end
guidata(hObject, handles);
