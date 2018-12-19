function varargout = alternangui(varargin)

% Function for running Alternan GUI 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary
% ALTERNANGUI MATLAB code for alternangui.fig
%      ALTERNANGUI, by itself, creates a new ALTERNANGUI or raises the existing
%      singleton*.
%
%      H = ALTERNANGUI returns the handle to a new ALTERNANGUI or the handle to
%      the existing singleton*.
%
%      ALTERNANGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALTERNANGUI.M with the given input arguments.
%
%      ALTERNANGUI('Property','Value',...) creates a new ALTERNANGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alternangui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alternangui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alternangui

% Last Modified by GUIDE v2.5 26-Sep-2017 15:52:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @alternangui_OpeningFcn, ...
    'gui_OutputFcn',  @alternangui_OutputFcn, ...
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


% --- Executes just before alternangui is made visible.
function alternangui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alternangui (see VARARGIN)

% Choose default command line output for alternangui
handles.output = hObject;
wb=waitbar(0.5,'Performing Alternans Calcualtions');
handles.output = hObject;
h = findobj('Tag','ElectroMap');
%% get data from ElectroMap
g1data = guidata(h);
handles.row=[];
handles.col=[];
handles.section_choice=get(g1data.listbox2,'Value');
handles.section=g1data.section;
handles.pointcount=0;
handles.images=g1data.images;
handles.mask=g1data.mask;
handles.cm=['b','r','g','m','y']
handles.framerate=str2num(get(g1data.framerate,'String'));
handles.AP=str2num(get(g1data.t,'String'));
handles.exposure=1/str2num(get(g1data.framerate,'String'));
[handles.rows, handles.cols, handles.num]= size(handles.images(:,:,:));
handles.time=[1:1:handles.num];
handles.time=(handles.time)*handles.exposure;
minpeakdist = str2num(get(g1data.minpeak,'String'));
handles.minpeakdist = ceil(minpeakdist/(1/str2num(get(g1data.framerate,'String'))));
handles.tfilt=get(g1data.tfilt,'Value');
[handles.pks,handles.locs] = findpeaks(g1data.averages, 'MINPEAKHEIGHT', g1data.peakheight, 'MINPEAKDISTANCE', g1data.minpeakdist);
handles.before=str2num(get(g1data.beforeGUI,'String'));
handles.after=str2num(get(g1data.afterGUI,'String'));
axes(handles.movieaxes)
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
handles.roichosen=0;
handles.isZoomed=0;
handles.mc=1;

%% setup small and big locs matrixs
%small
handles.ovsmalllocs=NaN(handles.rows,handles.cols,handles.num);handles.ovsmallPks=NaN(handles.rows,handles.cols,handles.num);handles.ovsmallloads_t=NaN(handles.rows,handles.cols,handles.num);handles.ovsmallloads=NaN(handles.rows,handles.cols,handles.num);
%big
handles.ovbiglocs=NaN(handles.rows,handles.cols,handles.num);handles.ovbigPks=NaN(handles.rows,handles.cols,handles.num);handles.ovbigloads_t=NaN(handles.rows,handles.cols,handles.num);handles.ovbigloads=NaN(handles.rows,handles.cols,handles.num);
%% do calcs

% Pixel maps
for row=1:handles.rows
    for col=1:handles.cols
        if handles.mask(row,col) ~= 0
            sig = (squeeze((handles.images(row,col,:))));
            [norms,pixa,t_pixa,pixl,t_pixl,pixap,amp_m,load_m,apd_m,...
                smalllocs,biglocs,smallPks,bigPks,smallloads_t,bigloads_t,smallloads,bigloads]...
                =alfn(sig,handles.framerate,handles.pks,handles.locs,handles.before,handles.after,handles.tfilt,handles.AP);
            handles.normsig(row,col,:)=squeeze(norms);
            handles.pixamp(row,col,:)=squeeze(pixa);
            handles.t_pixamp(row,col,:)=squeeze(t_pixa);
            handles.pixload(row,col,:)=squeeze(pixl);
            handles.t_pixload(row,col,:)=squeeze(t_pixl);
            handles.pixapd(row,col,:)=squeeze(pixap);
            handles.amp_map(row,col,:)=squeeze(amp_m);
            handles.load_map(row,col,:)=squeeze(load_m);
            handles.apd_map(row,col,:)=squeeze(apd_m);
            
            %small stuff
            for i =1:length(smalllocs)
                handles.ovsmalllocs(row,col,i)=squeeze(smalllocs(i));
                handles.ovsmallPks(row,col,i)=squeeze(smallPks(i));
                handles.ovsmallloads_t(row,col,i)=squeeze(smallloads_t(i));
                handles.ovsmallloads(row,col,i)=squeeze(smallloads(i));
            end
            %big stuff
            
            for i =1:length(biglocs)
                handles.ovbiglocs(row,col,i)=squeeze(biglocs(i));
                handles.ovbigPks(row,col,i)=squeeze(bigPks(i));
                handles.ovbigloads_t(row,col,i)=squeeze(bigloads_t(i));
                handles.ovbigloads(row,col,i)=squeeze(bigloads(i));
            end
        end
    end
end

%Averages
[handles.av_normsig,handles.av_pixamp,handles.av_t_pixamp,handles.av_pixload,handles.av_t_pixload,handles.av_pixapd,...
handles.av_amp_map,handles.av_load_map,handles.av_apd_map,...
handles.av_smalllocs,handles.av_biglocs,handles.av_smallPks,handles.av_bigPks,handles.av_smallloads_t,handles.av_bigloads_t,handles.av_smallloads,handles.av_bigloads]...
=alfn(imcomplement(g1data.averages),handles.framerate,handles.pks,handles.locs,handles.before,handles.after,handles.tfilt,handles.AP);

%results box
handles.L=mean(handles.av_bigPks)
handles.D=mean(handles.av_smallloads)
handles.S=mean([handles.av_smallPks-handles.D])

handles.release_alt=1-(handles.S/handles.L);
handles.load_alt=handles.D/handles.L;


rstr=[num2str(mean(handles.av_amp_map)),'+/-',num2str(std(handles.av_amp_map)),' %'];
rstr=[rstr newline num2str(handles.release_alt)];
rstr=[rstr newline num2str(handles.load_alt)];
rstr=[rstr newline];
rstr=[rstr newline [num2str(mean(handles.av_apd_map)),'+/-',num2str(std(handles.av_apd_map)),' ms']];

set(handles.resultstext,'String',rstr);

if size(handles.amp_map,1) ~= size(handles.mask,1) || size(handles.amp_map,2) ~= size(handles.mask,2)
    for row=1:handles.rows
        for col=1:handles.cols
            if handles.mask(row,col) == 0
                handles.normsig(row,col,:)=zeros(1,length(squeeze(norms)));
                handles.pixamp(row,col,:)=zeros(1,length(squeeze(pixa)));
                handles.t_pixamp(row,col,:)=zeros(1,length(squeeze(t_pixa)));
                handles.pixload(row,col,:)=zeros(1,length(squeeze(pixl)));
                handles.t_pixload(row,col,:)=zeros(1,length(squeeze(t_pixl)));
                handles.pixapd(row,col,:)=zeros(1,length(squeeze(pixap)));
                handles.amp_map(row,col,:)=zeros(1,length(squeeze(amp_m)));
                handles.load_map(row,col,:)=zeros(1,length(squeeze(load_m)));
                handles.apd_map(row,col,:)=zeros(1,length(squeeze(apd_m)));
            end
        end
    end
end

%medfiit maps
for i =1:length(handles.amp_map(1,1,:))
    handles.amp_map(:,:,i)=medfilt2(handles.amp_map(:,:,i));
    handles.load_map(:,:,i)=medfilt2(handles.load_map(:,:,i));
    handles.apd_map(:,:,i)=medfilt2(handles.apd_map(:,:,i));
end

set(handles.slider1,'Max',length(handles.amp_map(1,1,:)));
set(handles.slider1,'SliderStep',[1/(length(handles.amp_map(1,1,:))) 0.1]);
set(handles.apdmax,'String',num2str(max(max(max(handles.apd_map)))));
set(handles.apdmin,'String',num2str(min(min(min(handles.apd_map)))));
guidata(hObject, handles);

slider1_Callback(hObject, eventdata, handles)
Tracesource_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
delete(wb)
guidata(hObject, handles);

% UIWAIT makes alternangui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alternangui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in save_movie.
function save_movie_Callback(hObject, eventdata, handles)
handles.output = hObject;
[filename,pathname] = uiputfile({'*.gif'}, 'Save Movie');
fname=[pathname,filename];
imax=get(handles.slider1,'Max');
for i =1:imax
    set(handles.slider1,'Value',i)
    slider1_Callback(hObject, eventdata, handles)
    drawnow
    axes(handles.axes4);
    h=gca;
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,fname,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,fname,'gif','WriteMode','append');
    end
    handles = guidata(hObject);
end
% hObject    handle to save_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_values.
function save_values_Callback(hObject, eventdata, handles)
% hObject    handle to save_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
[filename,pathname] = uiputfile({'*.csv';'*.txt';'*.mat'}, 'Save Values');
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
if get(handles.Tracesource,'Value') == 1
    Beat=[1:length(handles.av_pixamp)]';
    Load=[handles.av_pixload]';
    Peak=[handles.av_pixamp];
    APD=[handles.av_pixapd];
    T=table(Beat,Load,Peak,APD);
    if strcmp('.csv',ext) == 1 || strcmp('.txt',ext) == 1
        writetable(T,file,'Delimiter',',','WriteVariableNames',1);
    end
end

if get(handles.Tracesource,'Value') == 3
    Beat=[1:length(handles.av_pixamp)]';
    Load=[handles.roi_pixload]';
    Peak=[handles.roi_pixamp];
    APD=[handles.roi_pixapd];
    T=table(Beat,Load,Peak,APD);
    if strcmp('.csv',ext) == 1 || strcmp('.txt',ext) == 1
        writetable(T,file,'Delimiter',',','WriteVariableNames',1);
    end
end
handles = guidata(hObject);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
axes(handles.movieaxes)
i=get(handles.slider1,'Value');
if i == 0
    i= 0.1
end
i=ceil(i)
handles.mc=i;
%% Which map
if get(handles.popupmenu1,'Value') == 1
    A=handles.amp_map(:,:,i);
elseif get(handles.popupmenu1,'Value') == 2
    A=handles.load_map(:,:,i);
elseif get(handles.popupmenu1,'Value') == 3 && i < length(handles.apd_map(1,1,:))
    A=handles.apd_map(:,:,i);
end
A(handles.mask == 0) = NaN;
%ddd=[[1 1 1];jet];
ddd=jet;
A
if get(handles.popupmenu1,'Value') == 1 || get(handles.popupmenu1,'Value') == 2
    if get(handles.amp_scale,'Value') == 1
        imshow(A, [],'ColorMap',ddd,'InitialMagnification', 800);
    end
    if get(handles.amp_scale,'Value') == 2 && get(handles.popupmenu1,'Value') == 1
        imshow(A, [min(min(min(handles.amp_map(:,:,:)))) max(max(max(handles.amp_map(:,:,:))))],'ColorMap',ddd,'InitialMagnification', 800);
    end
    if get(handles.amp_scale,'Value') == 2 && get(handles.popupmenu1,'Value') == 2
        imshow(A, [min(min(min(handles.load_map(:,:,:)))) max(max(max(handles.load_map(:,:,:))))],'ColorMap',ddd,'InitialMagnification', 800);
    end
    if get(handles.amp_scale,'Value') == 3
        imshow(A, [str2num(get(handles.ampmin,'String')) str2num(get(handles.ampmax,'String'))],'ColorMap',ddd,'InitialMagnification', 800);
    end
end

if get(handles.popupmenu1,'Value') == 3
    if get(handles.apd_scale,'Value') == 1
        imshow(A, [],'ColorMap',ddd,'InitialMagnification', 800);
    end
    if get(handles.apd_scale,'Value') == 2
        imshow(A, [min(min(min(handles.apd_map(:,:,:)))) max(max(max(handles.apd_map(:,:,:))))],'ColorMap',ddd,'InitialMagnification', 800);
    end
    if get(handles.apd_scale,'Value') == 3
        imshow(A, [str2num(get(handles.apdmin,'String')) str2num(get(handles.apdmax,'String'))],'ColorMap',ddd,'InitialMagnification', 800);
    end
end

Tracesource_Callback(hObject, eventdata, handles)

%% colorbar

freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
colormap('jet');
hcb.Location = 'southoutside'
ax = gca;
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;


hcb.TicksMode='manual';hcb.TickLabelsMode='manual';

%amp and load maps
if get(handles.popupmenu1,'Value') == 1 || get(handles.popupmenu1,'Value') == 2
    if get(handles.amp_scale,'Value') == 2
        if get(handles.popupmenu1,'Value') == 1
            cmin=min(min(min(handles.amp_map(:,:,:)))); cmax = max(max(max(handles.amp_map(:,:,:))));
        elseif  get(handles.popupmenu1,'Value') == 2
            cmin=min(min(min(handles.load_map(:,:,:)))); cmax = max(max(max(handles.load_map(:,:,:))));
        end
    elseif get(handles.amp_scale,'Value') == 1
        cmin=min(min(min(A(:,:)))); cmax = max(max(max(A(:,:))));
    elseif get(handles.amp_scale,'Value') == 3
        cmin=str2num(get(handles.ampmin,'String')); cmax=str2num(get(handles.ampmax,'String'));
    end
    stepp=(cmax-cmin)/5;
    hcb.TickLabels=[cmin:stepp:cmax];
    hcb.Ticks=[0.01,0.2:0.2:1];
    hcb.Label.String='Amplitude Alternan Ratio';
    axis off
end

if get(handles.popupmenu1,'Value') == 3
    if get(handles.apd_scale,'Value') == 1
        cmin=min(min(min(A(:,:)))); cmax = max(max(max(A(:,:))));
    elseif get(handles.apd_scale,'Value') == 2
        cmin=min(min(min(handles.apd_map(:,:,:)))); cmax = max(max(max(handles.apd_map(:,:,:))));
    elseif get(handles.apd_scale,'Value') == 3
        cmin=str2num(get(handles.apdmin,'String')); cmax=str2num(get(handles.apdmax,'String'));
    end
    stepp=(cmax-cmin)/5
    hcb.TickLabels=[cmin:stepp:cmax]
    hcb.Ticks=[0.01,0.2:0.2:1];
    hcb.Label.String='APD/CaD Alternan Amplitude (ms)';
    axis off
end


%%overaly POI
% add poi
axes(handles.movieaxes);
if handles.pointcount ~=0 && get(handles.Tracesource,'Value') == 2
    hold on
    plot(handles.col,handles.row,'+','MarkerSize',20,'LineWidth',5,'Color',handles.cm(2));
end

%%overlay ROI
if get(handles.Tracesource,'Value') == 3 && handles.roichosen == 1
    axes(handles.movieaxes);
    hold on
    for i=1:size(handles.boundaries,1)
        plot(handles.boundaries{i}(:,2),handles.boundaries{i}(:,1),'r','LineWidth',2);
    end
end
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


% --- Executes on button press in playbutton.
function playbutton_Callback(hObject, eventdata, handles)
% hObject    handle to playbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
handles = guidata(hObject);
imax=get(handles.slider1,'Max');
for i =1:imax
    set(handles.slider1,'Value',i)
    slider1_Callback(hObject, eventdata, handles)
    pause(0.05)
    handles = guidata(hObject);
end
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_trace.
function save_trace_Callback(hObject, eventdata, handles)
% hObject    handle to save_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function getMousePositionOnImage(hObject, eventdata, handles)
handles = guidata(hObject);
disp('hi')
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
cursorPoint = get(handles.movieaxes, 'CurrentPoint');
curX = cursorPoint(1,1);
curY = cursorPoint(1,2);
rrow=floor(curY);
rcol=floor(curX);
xLimits = get(handles.movieaxes, 'xlim');
yLimits = get(handles.movieaxes, 'ylim');
[rows, cols, num] = size(handles.images(:,:,:));

if get(handles.Tracesource,'Value') == 2
    if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits)) && handles.mask(rrow,rcol)~=0
        handles.pointcount=1;
        if handles.pointcount == 6
            handles.pointcount = 1;
            
        end
        handles.row=rrow
        handles.col=rcol
        rrow
        rcol
        guidata(hObject, handles);
        
        
        axes(handles.movieaxes)
        cla
        handles.output = hObject;
        slider1_Callback(hObject, eventdata, handles)
        handles = guidata(hObject);
        %     axes(handles.movieaxes)
        %     plot(handles.col(handles.pointcount),handles.row(handles.pointcount),'+','MarkerSize',20,'LineWidth',5,'Color',handles.cm(2));
        
    end
end


Tracesource_Callback(hObject, eventdata, handles)

guidata(hObject, handles);


% --- Executes on selection change in Tracetype.
function Tracetype_Callback(hObject, eventdata, handles)
% hObject    handle to Tracetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
Tracesource_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: contents = cellstr(get(hObject,'String')) returns Tracetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Tracetype


% --- Executes during object creation, after setting all properties.
function Tracetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tracetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in amp_scale.
function amp_scale_Callback(hObject, eventdata, handles)
% hObject    handle to amp_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: contents = cellstr(get(hObject,'String')) returns amp_scale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from amp_scale


% --- Executes during object creation, after setting all properties.
function amp_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amp_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ampmin_Callback(hObject, eventdata, handles)
% hObject    handle to ampmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: get(hObject,'String') returns contents of ampmin as text
%        str2double(get(hObject,'String')) returns contents of ampmin as a double


% --- Executes during object creation, after setting all properties.
function ampmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ampmax_Callback(hObject, eventdata, handles)
% hObject    handle to ampmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: get(hObject,'String') returns contents of ampmax as text
%        str2double(get(hObject,'String')) returns contents of ampmax as a double


% --- Executes during object creation, after setting all properties.
function ampmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apd_scale.
function apd_scale_Callback(hObject, eventdata, handles)
% hObject    handle to apd_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: contents = cellstr(get(hObject,'String')) returns apd_scale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apd_scale


% --- Executes during object creation, after setting all properties.
function apd_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apd_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apdmin_Callback(hObject, eventdata, handles)
% hObject    handle to apdmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: get(hObject,'String') returns contents of apdmin as text
%        str2double(get(hObject,'String')) returns contents of apdmin as a double


% --- Executes during object creation, after setting all properties.
function apdmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apdmax_Callback(hObject, eventdata, handles)
% hObject    handle to apdmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Hints: get(hObject,'String') returns contents of apdmax as text
%        str2double(get(hObject,'String')) returns contents of apdmax as a double


% --- Executes during object creation, after setting all properties.
function apdmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Tracesource.
function Tracesource_Callback(hObject, eventdata, handles)
% hObject    handle to Tracesource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = hObject;
%slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
pp=handles.mc+1;

%% AVERAGES (both zoom and no zoom)
if get(handles.Tracesource,'Value') == 1
    handles.L=mean(handles.av_bigPks);
    handles.D=mean(handles.av_smallloads);
    handles.S=mean([handles.av_smallPks-handles.D]);
    handles.release_alt=1-(handles.S/handles.L)
    handles.load_alt=handles.D/handles.L
    traceopt=get(handles.Tracetype,'Value');
    
    if traceopt == 1
        axes(handles.trace)
        cla reset
        plot(handles.av_normsig,'k');
        hold on
        ylim([0 1])
        plot(handles.av_t_pixamp,handles.av_pixamp,'ob');
        plot(handles.av_t_pixload,handles.av_pixload,'xb');
        hx = plot([(handles.locs(pp+1)-handles.before) (handles.locs(pp+1)-handles.before)], ylim, 'Color', 'k');
        hx2 = plot([(handles.locs(pp+1)+handles.after) (handles.locs(pp+1)+handles.after)], ylim, 'Color', 'k');
    end
    
    if traceopt == 2
        axes(handles.trace)
        cla reset
        yyaxis right
        plot(handles.av_normsig,'k');
        hold on
        plot(handles.av_t_pixamp,handles.av_pixamp,'ob');
        plot(handles.av_t_pixload,handles.av_pixload,'xb');
        yyaxis left
        plot(handles.av_t_pixamp(2:end),handles.av_amp_map,'rx');
        plot(handles.av_t_pixamp(2:end),handles.av_load_map,'ro');
        hx = plot([(handles.locs(pp+1)-handles.before) (handles.locs(pp+1)-handles.before)], ylim, 'Color', 'k');
        hx2 = plot([(handles.locs(pp+1)+handles.after) (handles.locs(pp+1)+handles.after)], ylim, 'Color', 'k');
    end
    
    if traceopt == 3
        axes(handles.trace)
        cla reset
        yyaxis right
        plot(handles.av_normsig,'k');
        ax=gca; ax.YColor='k'
        xlabel('Time (ms)')
        ylabel('Normalised Fluorescence')
        yyaxis left
        plot(handles.av_t_pixamp,handles.av_pixapd,'xb-');
        hold on
        plot(handles.av_t_pixamp(2:end),handles.av_apd_map,'rx-');
        ylabel('APD/APD Alternan (ms)')
        ax=gca; ax.YColor= 'r';
        hx = plot([(handles.locs(pp+1)-handles.before) (handles.locs(pp+1)-handles.before)], ylim, 'Color', 'k');
        hx2 = plot([(handles.locs(pp+1)+handles.after) (handles.locs(pp+1)+handles.after)], ylim, 'Color', 'k');
    end
    
    
    if traceopt == 4
        axes(handles.trace)
        cla reset
        plot(handles.av_normsig,'k');
        hold on
        for i=1:length(handles.av_smalllocs)
            plot([handles.av_smalllocs(i)-handles.before:handles.av_smalllocs(i)+handles.after],handles.av_normsig(handles.av_smalllocs(i)-handles.before:handles.av_smalllocs(i)+handles.after),'b');
        end
        for i=1:length(handles.av_biglocs)
            ix = find(handles.av_smallloads_t>handles.av_biglocs(i));
            if isempty(ix) == 0
                ix = handles.av_smallloads_t(ix(1))-handles.av_biglocs(i);
                if ix <= handles.after;
                    next_load=ix
                else next_load=handles.after
                end
            else next_load=handles.after
            end
            plot([handles.av_biglocs(i)-handles.before:handles.av_biglocs(i)+next_load],handles.av_normsig(handles.av_biglocs(i)-handles.before:handles.av_biglocs(i)+next_load),'r'); %to next load time rather than +after to aviod unwanted
        end
        
        plot(handles.av_smallloads_t,handles.av_smallloads,'ob')
        plot(handles.av_bigloads_t,handles.av_bigloads,'or')
        hx = plot([(handles.locs(pp+1)-handles.before) (handles.locs(pp+1)-handles.before)], ylim, 'Color', 'k');
        hx2 = plot([(handles.locs(pp+1)+handles.after) (handles.locs(pp+1)+handles.after)], ylim, 'Color', 'k');
        
        ax=gca;
        xx=get(ax,'Xlim');
        yy=get(ax,'Ylim');
        
        hL = plot([xx(1) xx(2)], [handles.L handles.L],'LineStyle', '--', 'Color', 'r');
        hD = plot([xx(1) xx(2)], [handles.D handles.D],'LineStyle', '--', 'Color', 'k');
        hS = plot([xx(1) xx(2)], [handles.S handles.S],'LineStyle', '--', 'Color', 'b');
        
    end
    %update results box
    rstr=[num2str(mean(handles.av_amp_map)),'+/-',num2str(std(handles.av_amp_map)),' %'];
    rstr=[rstr newline num2str(handles.release_alt)];
    rstr=[rstr newline num2str(handles.load_alt)];
    rstr=[rstr newline [' ']];
    rstr=[rstr newline [num2str(mean(handles.av_apd_map)),'+/-',num2str(std(handles.av_apd_map)),' ms']];
    set(handles.resultstext,'String',rstr);
end


%% POI no zoom
if get(handles.Tracesource,'Value') == 2 && handles.isZoomed == 0;
    traceopt=get(handles.Tracetype,'Value')
    sig = (squeeze((handles.images(handles.row,handles.col,:))));
    [norms,pixa,t_pixa,pixl,t_pixl,pixap,amp_m,load_m,apd_m,...
        smalllocs,biglocs,smallPks,bigPks,smallloads_t,bigloads_t,smallloads,bigloads]...
        =alfn(sig,handles.framerate,handles.pks,handles.locs,handles.before,handles.after,handles.tfilt,handles.AP);
    
    % stuff for big/small peak loads
    handles.pixsmalllocs=squeeze(smalllocs);
    handles.pixbiglocs=squeeze(biglocs);
    handles.pixsmallPks=squeeze(smallPks);
    handles.pixbigPks=squeeze(bigPks);
    handles.pixsmallloads_t=squeeze(smallloads_t);
    handles.pixbigloads_t=squeeze(bigloads_t);
    handles.pixsmallloads=squeeze(smallloads);
    handles.pixbigloads=squeeze(bigloads);
    
    handles.poiL=mean(handles.pixbigPks);
    handles.poiD=mean(handles.pixsmallloads);
    handles.poiS=mean([handles.pixsmallPks-handles.D]);
    handles.poirelease_alt=1-(handles.poiS/handles.poiL);
    handles.poiload_alt=handles.poiD/handles.poiL;
    
    if traceopt == 1
        axes(handles.trace)
        cla reset
        plot(squeeze(handles.normsig(handles.row,handles.col,:)),'k');
        hold on
        plot(squeeze(handles.t_pixamp(handles.row,handles.col,:)),squeeze(handles.pixamp(handles.row,handles.col,:)),'ob');
        plot(squeeze(handles.t_pixload(handles.row,handles.col,:)),squeeze(handles.pixload(handles.row,handles.col,:)),'xb');
    end
    
    if traceopt == 2
        axes(handles.trace)
        cla reset
        yyaxis right
        plot(squeeze(handles.normsig(handles.row,handles.col,:)),'k');
        hold on
        plot(squeeze(handles.t_pixamp(handles.row,handles.col,:)),squeeze(handles.pixamp(handles.row,handles.col,:)),'ob');
        plot(squeeze(handles.t_pixload(handles.row,handles.col,:)),squeeze(handles.pixload(handles.row,handles.col,:)),'xb');
        yyaxis left
        plot(squeeze(handles.t_pixamp(handles.row,handles.col,2:end)),squeeze(handles.amp_map(handles.row,handles.col,:)),'rx');
        plot(squeeze(handles.t_pixamp(handles.row,handles.col,2:end)),squeeze(handles.load_map(handles.row,handles.col,:)),'ro');
        
    end
    
    if traceopt == 3
        axes(handles.trace)
        cla reset
        yyaxis right
        plot(squeeze(handles.normsig(handles.row,handles.col,:)),'k');
        ax=gca; ax.YColor='k'
        xlabel('Time (ms)')
        ylabel('Normalised Fluorescence')
        yyaxis left
        plot(squeeze(handles.t_pixamp(handles.row,handles.col,:)),squeeze(handles.pixapd(handles.row,handles.col,:)),'xb-');
        hold on
        plot(squeeze(handles.t_pixamp(handles.row,handles.col,2:end)),squeeze(handles.apd_map(handles.row,handles.col,:)),'rx-');
        ylabel('APD/APD Alternan (ms)')
        ax=gca; ax.YColor= 'r';
    end
    
    
    if traceopt == 4
        axes(handles.trace)
        cla reset
        plot(norms,'k');
        hold on
        for i=1:length(handles.pixsmalllocs)
            plot([handles.pixsmalllocs(i)-handles.before:handles.pixsmalllocs(i)+handles.after],norms(handles.pixsmalllocs(i)-handles.before:handles.pixsmalllocs(i)+handles.after),'b');
        end
        for i=1:length(handles.pixbiglocs)
            ix = find(handles.pixsmallloads_t>handles.pixbiglocs(i));
            if isempty(ix) == 0
                length(handles.pixsmallloads_t)
                length(handles.pixbiglocs)
                ix1 = handles.pixsmallloads_t(ix(1))-handles.pixbiglocs(i);
                if ix1 <= handles.after;
                    next_load=ix1
                else next_load=handles.after
                end
            else next_load=handles.after
            end
            plot([handles.pixbiglocs(i)-handles.before:handles.pixbiglocs(i)+next_load],norms(handles.pixbiglocs(i)-handles.before:handles.pixbiglocs(i)+next_load),'r'); %to next load time rather than +after to aviod unwanted
        end
        
        
        plot(handles.pixsmallloads_t,handles.pixsmallloads,'ob')
        plot(handles.pixbigloads_t,handles.pixbigloads,'or')
        plot(handles.pixsmalllocs,handles.pixsmallPks,'xb')
        plot(handles.pixbiglocs,handles.pixbigPks,'xr')
        ax=gca;
        xx=get(ax,'Xlim');
        yy=get(ax,'Ylim');
        
        hL = plot([xx(1) xx(2)], [handles.poiL handles.poiL],'LineStyle', '--', 'Color', 'r');
        hD = plot([xx(1) xx(2)], [handles.poiD handles.poiD],'LineStyle', '--', 'Color', 'k');
        hS = plot([xx(1) xx(2)], [handles.poiS handles.poiS],'LineStyle', '--', 'Color', 'b');
        
    end
    %set(handles.poiamp,'String',[num2str(mean(squeeze(handles.amp_map(handles.row,handles.col,:)))),' +/- ', num2str(std(squeeze(handles.amp_map(handles.row,handles.col,:)))),'%']);
    %set(handles.poiload,'String',[num2str(mean(squeeze(handles.load_map(handles.row,handles.col,:)))),' +/- ', num2str(std(squeeze(handles.load_map(handles.row,handles.col,:)))),'%']);
    %set(handles.poiapd,'String',[num2str(mean(squeeze(handles.apd_map(handles.row,handles.col,:)))),' +/- ', num2str(std(squeeze(handles.apd_map(handles.row,handles.col,:)))),'ms']);
    
    
    rstr=[num2str(mean(handles.amp_map(handles.row,handles.col,:))),'+/-',num2str(std(handles.amp_map(handles.row,handles.col,:))),' %'];
    rstr=[rstr newline num2str(handles.poirelease_alt)];
    rstr=[rstr newline num2str(handles.poiload_alt)];
    rstr=[rstr newline];
    rstr=[rstr newline [num2str(mean(handles.apd_map(handles.row,handles.col,:))),'+/-',num2str(std(handles.apd_map(handles.row,handles.col,:))),' ms']];
    
    set(handles.resultstext,'String',rstr);
    
end

%% ROI
if get(handles.Tracesource,'Value') == 3
    check=handles.roichosen
    if handles.roichosen == 0;
        figure,
        imshow(handles.images(:,:,1), [],'InitialMagnification', 800)
        title('Select ROI and press enter');
        handles.BWpoly = [];
        handles.BWpoly=uint32(roipoly);
        polyfig=gcf;
        close(polyfig);
        
        % create average signal from this area
        for j=1:handles.num
            A1=handles.images(:,:,j);
            B = A1.*handles.BWpoly;
            avs(:,j) = mean2(B);
        end
        avs=avs-min(avs);
        avs=avs/max(avs);
        handles.roisig=avs;
        handles.roichosen=1;
        
        [handles.boundaries] = bwboundaries(handles.BWpoly);
        slider1_Callback(hObject, eventdata, handles);
    end
    
    [handles.roi_normsig,handles.roi_pixamp,handles.roi_t_pixamp,handles.roi_pixload,handles.roi_t_pixload,handles.roi_pixapd,handles.roi_amp_map,handles.roi_load_map,handles.roi_apd_map,...
        handles.roi_smalllocs,handles.roi_biglocs,handles.roi_smallPks,handles.roi_bigPks,handles.roi_smallloads_t,handles.roi_bigloads_t,handles.roi_smallloads,handles.roi_bigloads]...
        =alfn(handles.roisig,handles.framerate,handles.pks,handles.locs,handles.before,handles.after,handles.tfilt,handles.AP)
    handles.roiL=mean(handles.roi_bigPks);
    
    handles.roiD=mean(handles.roi_smallloads);
    handles.roiS=mean([handles.roi_smallPks-handles.roiD]);
    handles.roi_release_alt=1-(handles.roiS/handles.roiL);
    handles.roi_load_alt=handles.roiD/handles.roiL;
    
    traceopt=get(handles.Tracetype,'Value');
    
    
    %roi results
    rstr=[num2str(mean(handles.roi_amp_map)),'+/-',num2str(std(handles.roi_amp_map)),' %'];
    rstr=[rstr newline num2str(handles.roi_release_alt)];
    rstr=[rstr newline num2str(handles.roi_load_alt)];
    rstr=[rstr newline];
    rstr=[rstr newline [num2str(mean(handles.roi_apd_map)),'+/-',num2str(std(handles.roi_apd_map)),' ms']];
    
    set(handles.resultstext,'String',rstr);
    
    if traceopt == 1
        axes(handles.trace)
        cla reset
        plot(handles.roi_normsig,'k');
        hold on
        ylim([0 1])
        plot(handles.roi_t_pixamp,handles.roi_pixamp,'ob');
        plot(handles.roi_t_pixload,handles.roi_pixload,'xb');
    end
    
    if traceopt == 2
        axes(handles.trace)
        cla reset
        yyaxis right
        plot(handles.roi_normsig,'k');
        hold on
        plot(handles.roi_t_pixamp,handles.roi_pixamp,'ob');
        plot(handles.roi_t_pixload,handles.roi_pixload,'xb');
        hold on
        yyaxis left
        plot(handles.roi_t_pixamp(2:end),handles.roi_amp_map,'rx');
        plot(handles.roi_t_pixamp(2:end),handles.roi_load_map,'ro');
    end
    
    if traceopt == 3
        axes(handles.trace)
        cla reset
        yyaxis right
        plot(handles.roi_normsig,'k');
        ax=gca; ax.YColor='k'
        xlabel('Time (ms)')
        ylabel('Normalised Fluorescence')
        yyaxis left
        plot(handles.roi_t_pixamp,handles.roi_pixapd,'xb-');
        hold on
        plot(handles.roi_t_pixamp(2:end),handles.roi_apd_map,'rx-');
        ylabel('APD/APD Alternan (ms)')
        ax=gca; ax.YColor= 'r';
    end
    
    if traceopt == 4
        axes(handles.trace)
        cla reset
        plot(handles.roi_normsig,'k');
        hold on
        for i=1:length(handles.roi_smalllocs)
            plot([handles.roi_smalllocs(i)-handles.before:handles.roi_smalllocs(i)+handles.after],handles.roi_normsig(handles.roi_smalllocs(i)-handles.before:handles.roi_smalllocs(i)+handles.after),'b');
        end
        for i=1:length(handles.roi_biglocs)
            ix = find(handles.roi_smallloads_t>handles.roi_biglocs(i));
            if isempty(ix) == 0
                ix = handles.roi_smallloads_t(ix(1))-handles.roi_biglocs(i);
                if ix <= handles.after;
                    next_load=ix
                else next_load=handles.after
                end
            else next_load=handles.after
            end
            plot([handles.roi_biglocs(i)-handles.before:handles.roi_biglocs(i)+next_load],handles.roi_normsig(handles.roi_biglocs(i)-handles.before:handles.roi_biglocs(i)+next_load),'r'); %to next load time rather than +after to aviod unwanted
        end
        
        plot(handles.roi_smallloads_t,handles.roi_smallloads,'xb')
        plot(handles.roi_bigloads_t,handles.roi_bigloads,'xr')
        plot(handles.roi_smalllocs,handles.roi_smallPks,'ob')
        plot(handles.roi_biglocs,handles.roi_bigPks,'or')
    end
end
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns Tracesource contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Tracesource


% --- Executes during object creation, after setting all properties.
function Tracesource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tracesource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in removepoint.
function removepoint_Callback(hObject, eventdata, handles)
% hObject    handle to removepoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
handles.pointcount=0;
slider1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

% --- Executes on button press in removeregion.
function removeregion_Callback(hObject, eventdata, handles)
handles.output = hObject;
handles.roichosen=0;
slider1_Callback(hObject, eventdata, handles)
guidata(hObject);
% hObject    handle to removeregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in zoomanal.
function zoomanal_Callback(hObject, eventdata, handles)
% hObject    handle to zoomanal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('__________________________________')
handles
axes(handles.trace);
handles.newlim=get(gca,'XLim');
handles.newlim(1)=floor(handles.newlim(1));
handles.newlim(2)=ceil(handles.newlim(2));
n1=handles.newlim(1)-handles.before;
n2=handles.newlim(2)+handles.after;
tracesource=get(handles.Tracesource,'Value');
if tracesource == 1
    
    %% find indicies between n1 and n2
    %AMP AND APD
    amp_ind1=find(handles.av_t_pixamp>n1);amp_ind2=find(handles.av_t_pixamp<n2);
    amp_ind1=amp_ind1(1)
    amp_ind2=amp_ind2(end)
    if amp_ind2 > length(handles.av_amp_map)
        amp_ind2=length(handles.av_amp_map)
    end
    %Large Peaks
    big_ind1=find(handles.av_biglocs>n1);big_ind2=find(handles.av_biglocs<n2);
    big_ind1=big_ind1(1);
    big_ind2=big_ind2(end);
    %Small Peaks
    small_ind1=find(handles.av_smalllocs>n1);small_ind2=find(handles.av_smalllocs<n2);
    small_ind1=small_ind1(1) ;
    small_ind2=small_ind2(end);
    %Small Loads
    smallloads_ind1=find(handles.av_smallloads_t>n1);smallloads_ind2=find(handles.av_smallloads_t<n2);
    smallloads_ind1=smallloads_ind1(1);
    smallloads_ind2=smallloads_ind2(end);
    
    %% cut to new lims
    handles.zav_amp_map=handles.av_amp_map(amp_ind1:amp_ind2);
    handles.zav_apd_map=handles.av_apd_map(amp_ind1:amp_ind2);
    handles.zav_bigPks=handles.av_bigPks(big_ind1:big_ind2);
    handles.zav_smallPks=handles.av_smallPks(small_ind1:small_ind2);
    handles.zav_smallloads=handles.av_smallloads(smallloads_ind1:smallloads_ind2);
    
    %% LDS and release alts
    handles.zL=mean(handles.zav_bigPks);
    handles.zD=mean(handles.zav_smallloads);
    handles.zS=mean([handles.zav_smallPks-handles.zD]);
    handles.zrelease_alt=1-(handles.zS/handles.zL)
    handles.zload_alt=handles.zD/handles.zL
    
    %% Update results string
    rstr=[num2str(mean(handles.zav_amp_map)),'+/-',num2str(std(handles.zav_amp_map)),' %'];
    rstr=[rstr newline num2str(handles.zrelease_alt)];
    rstr=[rstr newline num2str(handles.zload_alt)];
    rstr=[rstr newline [' ']];
    rstr=[rstr newline [num2str(mean(handles.zav_apd_map)),'+/-',num2str(std(handles.zav_apd_map)),' ms']];
    set(handles.resultstext,'String',rstr);
    
end

if tracesource == 2
    %% find indicies between n1 and n2
    %AMP AND APD
    amp_ind1=find(squeeze(handles.t_pixamp(handles.row,handles.col,:))>n1);amp_ind2=find(squeeze(handles.t_pixamp(handles.row,handles.col,:))<n2)
    amp_ind1=amp_ind1(1);
    amp_ind2=amp_ind2(end);
    if amp_ind2 > length(handles.av_amp_map)
        amp_ind2=length(handles.av_amp_map)
    end
    %Large Peaks
    big_ind1=find(squeeze(handles.ovbiglocs(handles.row,handles.col,:))>n1);big_ind2=find(squeeze(handles.ovbiglocs(handles.row,handles.col,:))<n2);
    big_ind1=big_ind1(1);
    big_ind2=big_ind2(end);
    %Small Peaks
    small_ind1=find(squeeze(handles.ovsmalllocs(handles.row,handles.col,:))>n1);small_ind2=find(squeeze(handles.ovsmalllocs(handles.row,handles.col,:))<n2);
    small_ind1=small_ind1(1) ;
    small_ind2=small_ind2(end);
    %Small Loads
    smallloads_ind1=find(squeeze(handles.ovsmallloads_t(handles.row,handles.col,:))>n1);smallloads_ind2=find(squeeze(handles.ovsmallloads_t(handles.row,handles.col,:))<n2);
    smallloads_ind1=smallloads_ind1(1);
    smallloads_ind2=smallloads_ind2(end);
    
    %% cut to new lims
    handles.zamp_map=squeeze(handles.amp_map(handles.row,handles.col,amp_ind1:amp_ind2));
    handles.zapd_map=squeeze(handles.apd_map(handles.row,handles.col,amp_ind1:amp_ind2));
    handles.zbigPks=squeeze(handles.ovbigPks(handles.row,handles.col,big_ind1:big_ind2));
    handles.zsmallPks=squeeze(handles.ovsmallPks(handles.row,handles.col,small_ind1:small_ind2));
    handles.zsmallloads=squeeze(handles.ovsmallloads(handles.row,handles.col,smallloads_ind1:smallloads_ind2));
    
    %% LDS and release alts
    handles.zL=mean(handles.zbigPks);
    handles.zD=mean(handles.zsmallloads);
    handles.zS=mean([handles.zsmallPks-handles.zD]);
    handles.zrelease_alt=1-(handles.zS/handles.zL)
    handles.zload_alt=handles.zD/handles.zL
    
    %% Update results string
    rstr=[num2str(mean(handles.zamp_map)),'+/-',num2str(std(handles.zamp_map)),' %'];
    rstr=[rstr newline num2str(handles.zrelease_alt)];
    rstr=[rstr newline num2str(handles.zload_alt)];
    rstr=[rstr newline [' ']];
    rstr=[rstr newline [num2str(mean(handles.zapd_map)),'+/-',num2str(std(handles.zapd_map)),' ms']];
    set(handles.resultstext,'String',rstr);
    
end

if tracesource == 3
    
    %% find indicies between n1 and n2
    %AMP AND APD
    amp_ind1=find(handles.roi_t_pixamp>n1);amp_ind2=find(handles.roi_t_pixamp<n2);
    amp_ind1=amp_ind1(1)
    amp_ind2=amp_ind2(end)
    if amp_ind2 > length(handles.roi_amp_map)
        amp_ind2=length(handles.roi_amp_map)
    end
    %Large Peaks
    big_ind1=find(handles.roi_biglocs>n1);big_ind2=find(handles.roi_biglocs<n2);
    big_ind1=big_ind1(1);
    big_ind2=big_ind2(end);
    %Small Peaks
    small_ind1=find(handles.roi_smalllocs>n1);small_ind2=find(handles.roi_smalllocs<n2);
    small_ind1=small_ind1(1) ;
    small_ind2=small_ind2(end);
    %Small Loads
    smallloads_ind1=find(handles.roi_smallloads_t>n1);smallloads_ind2=find(handles.roi_smallloads_t<n2);
    smallloads_ind1=smallloads_ind1(1);
    smallloads_ind2=smallloads_ind2(end);
    
    %% cut to new lims
    handles.zroi_amp_map=handles.roi_amp_map(amp_ind1:amp_ind2);
    handles.zroi_apd_map=handles.roi_apd_map(amp_ind1:amp_ind2);
    handles.zroi_bigPks=handles.roi_bigPks(big_ind1:big_ind2);
    handles.zroi_smallPks=handles.roi_smallPks(small_ind1:small_ind2);
    handles.zroi_smallloads=handles.roi_smallloads(smallloads_ind1:smallloads_ind2);
    
    %% LDS and release alts
    handles.zL=mean(handles.zroi_bigPks);
    handles.zD=mean(handles.zroi_smallloads);
    handles.zS=mean([handles.zroi_smallPks-handles.zD]);
    handles.zrelease_alt=1-(handles.zS/handles.zL)
    handles.zload_alt=handles.zD/handles.zL
    
    %% Update results string
    rstr=[num2str(mean(handles.zroi_amp_map)),'+/-',num2str(std(handles.zroi_amp_map)),' %'];
    rstr=[rstr newline num2str(handles.zrelease_alt)];
    rstr=[rstr newline num2str(handles.zload_alt)];
    rstr=[rstr newline [' ']];
    rstr=[rstr newline [num2str(mean(handles.zroi_apd_map)),'+/-',num2str(std(handles.zroi_apd_map)),' ms']];
    set(handles.resultstext,'String',rstr);
    
end
handles = guidata(hObject);
% Hint: get(hObject,'Value') returns toggle state of checkbox1

function checkbox1_Callback(hObject, eventdata, handles)
handles.output = hObject;
cho1=get(handles.checkbox1,'Value')
if cho1 == 0
    zoom xon
else zoom off
end


% --- Executes on button press in meanmap.
function meanmap_Callback(hObject, eventdata, handles)
% hObject    handle to meanmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
cho = questdlg('What kind of map?', ...
    'Map choice', ...
    'Amplitude','Release','Load','Load')
ddd=[[1,1,1];jet];

prompt = {'Enter start time (ms):','Enter end time (ms):'};
dlg_title = 'Frame Select';
num_lines = 1;
defaultans = {'0',num2str(length(handles.av_normsig))};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
f1=str2num(answer{1});f2=str2num(answer{2});
a1=find(handles.locs>f1);
a1=a1(1);
a2=find(handles.locs<f2);
a2=a2(end);
if a2>length(handles.amp_map(1,1,:))
    a2=length(handles.amp_map(1,1,:));
end
switch cho
    case 'Amplitude'
        for i =1:handles.rows
            for j=1:handles.cols
                length(handles.amp_map(1,1,:))
                a=handles.amp_map(i,j,a1:a2);
                A(i,j)=nanmean(a(a~=0));
            end
        end
        figure,
        if get(handles.popupmenu1,'Value') == 1 || get(handles.popupmenu1,'Value') == 2
            if get(handles.amp_scale,'Value') == 1
                imshow(A, [],'ColorMap',ddd,'InitialMagnification', 800);
            end
            if get(handles.amp_scale,'Value') == 2 && get(handles.popupmenu1,'Value') == 1
                imshow(A, [min(min(min(handles.amp_map(:,:,:)))) max(max(max(handles.amp_map(:,:,:))))],'ColorMap',ddd,'InitialMagnification', 1600);
            end
            if get(handles.amp_scale,'Value') == 2 && get(handles.popupmenu1,'Value') == 2
                imshow(A, [min(min(min(handles.load_map(:,:,:)))) max(max(max(handles.load_map(:,:,:))))],'ColorMap',ddd,'InitialMagnification', 1600);
            end
            if get(handles.amp_scale,'Value') == 3
                imshow(A, [str2num(get(handles.ampmin,'String')) str2num(get(handles.ampmax,'String'))],'ColorMap',ddd,'InitialMagnification', 1600);
            end
        end
        if get(handles.amp_scale,'Value') == 2
            if get(handles.popupmenu1,'Value') == 1
                cmin=min(min(min(handles.amp_map(:,:,:)))); cmax = max(max(max(handles.amp_map(:,:,:))));
            elseif  get(handles.popupmenu1,'Value') == 2
                cmin=min(min(min(handles.load_map(:,:,:)))); cmax = max(max(max(handles.load_map(:,:,:))));
            end
        elseif get(handles.amp_scale,'Value') == 1
            cmin=min(min(min(A(:,:)))); cmax = max(max(max(A(:,:))));
        elseif get(handles.amp_scale,'Value') == 3
            cmin=str2num(get(handles.ampmin,'String')); cmax=str2num(get(handles.ampmax,'String'));
        end
        %hcb=colorbar;
        %colormap('jet');
        %colormap('jet');
        %hcb.Label.String=' Average Amplitude Alternan (%)';
        
        for i =1:handles.rows
            for j=1:handles.cols
                length(handles.apd_map(1,1,:))
                b=handles.apd_map(i,j,a1:a2);
                B(i,j)=nanmean(b(b~=0));
            end
        end
        figure,
        imshow(B, [0 10],'ColorMap',ddd,'InitialMagnification', 1600);
        
        %amp and load maps
        
        
    case 'Release'
        for i=1:handles.rows
            for j=1:handles.cols
                Ds=[];
                Ls=[];
                Ss=[];
                if handles.mask(i,j) ~=0
                    %%find big and smallPks inside frames
                    s1=[];s2=[];b1=[];b2=[];
                    s1=find(handles.ovsmallloads_t(i,j,:)>f1);
                    s1=s1(1);
                    s2=find(handles.ovsmallloads_t(i,j,:)<f2);
                    s2=s2(end);
                    
                    b1=find(handles.ovbigloads_t(i,j,:)>f1);
                    b1=b1(1);
                    b2=find(handles.ovbigloads_t(i,j,:)<f2);
                    b2=b2(end);
                    
                    Ls=squeeze(handles.ovbigPks(i,j,b1:b2));
                    Ds=squeeze(handles.ovsmallloads(i,j,s1:s2));
                    Ss=squeeze(handles.ovsmallPks(i,j,s1:s2));
                    
                    
                    L=nanmean(Ls);
                    D=nanmean(Ds);
                    S1=nanmean(Ss);
                    S=S1-D;
                    A(i,j)=(1-S/L);
                else A(i,j)=NaN;
                end
            end
        end
        A=medfilt2(A);
        figure,
        imshow(A, [0 1],'ColorMap',ddd,'InitialMagnification', 1600);
        
    case 'Load'
        for i=1:handles.rows
            for j=1:handles.cols
                Ds=[];
                Ls=[];
                Ss=[];
                if handles.mask(i,j) ~=0
                    %%find big and smallPks inside frames
                    s1=[];s2=[];b1=[];b2=[];
                    s1=find(handles.ovsmallloads_t(i,j,:)>f1);
                    s1=s1(1);
                    s2=find(handles.ovsmallloads_t(i,j,:)<f2)
                    s2=s2(end);
                    
                    b1=find(handles.ovbigloads_t(i,j,:)>f1);
                    b1=b1(1);
                    b2=find(handles.ovbigloads_t(i,j,:)<f2);
                    b2=b2(end);
                    
                    Ls=squeeze(handles.ovbigPks(i,j,b1:b2));
                    Ds=squeeze(handles.ovsmallloads(i,j,s1:s2))
                    Ss=squeeze(handles.ovsmallPks(i,j,s1:s2));
                    
                    
                    L=nanmean(Ls)
                    D=nanmean(Ds)
                    S1=nanmean(Ss);
                    S=S1-D;
                    A(i,j)=(D/L);
                else A(i,j)=NaN;
                end
            end
        end
        A=medfilt2(A);
        figure,
        imshow(A, [-0.05 .3],'ColorMap',ddd,'InitialMagnification', 1600);
end



function apdcad_Callback(hObject, eventdata, handles)
% hObject    handle to apdcad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apdcad as text
%        str2double(get(hObject,'String')) returns contents of apdcad as a double


% --- Executes during object creation, after setting all properties.
function apdcad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdcad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
