function varargout = conduction(varargin)

% Function for running conduction GUI. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries

% Release Date - 
% For licence information, please see 'licsence.txt' at ...

% Last Updated -

% Update Summary

% CONDUCTION MATLAB code for conduction.fig
%      CONDUCTION, by itself, creates a new CONDUCTION or raises the existing
%      singleton*.
%
%      H = CONDUCTION returns the handle to a new CONDUCTION or the handle to
%      the existing singleton*.
%
%      CONDUCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONDUCTION.M with the given input arguments.
%
%      CONDUCTION('Property','Value',...) creates a new CONDUCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before conduction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to conduction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help conduction

% Last Modified by GUIDE v2.5 02-Dec-2017 17:09:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @conduction_OpeningFcn, ...
                   'gui_OutputFcn',  @conduction_OutputFcn, ...
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


% --- Executes just before conduction is made visible.
function conduction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to conduction (see VARARGIN)

% Choose default command line output for conduction
handles.output = hObject;
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
axes(handles.actmap);
wb=waitbar(0.5,'Transfering Activation Map');


handles.map=g1data.actmap;
handles.act_t=g1data.act_t;
handles.vout=g1data.vout;
handles.quivers_Xout=g1data.quivers_Xout;
handles.quivers_Yout=g1data.quivers_Yout;
handles.quivers_vxout=g1data.quivers_vxout;
handles.quivers_vyout=g1data.quivers_vyout;

delete(wb)



set(handles.actopt,'Value',get(g1data.velalgo,'Value'))
axes(handles.actmap);
jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [map] = handles.map;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(map));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(map, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
%set activation curve settings from electromap
set(handles.actmin,'String',get(g1data.actmin,'String'));
set(handles.actmax,'String',get(g1data.actmax,'String'));
set(handles.actquote,'String',get(g1data.actquote,'String'));

%
handles.point_A=[];
handles.point_B=[]; 

%set local vector settings from electromap

EMnums=g1data.rdata
set(handles.localmean,'String',num2str(EMnums(2,1)));
set(handles.localout,'Value',get(g1data.velout,'Value'))
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes conduction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = conduction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in actopt.
function actopt_Callback(hObject, eventdata, handles)
% hObject    handle to actopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
axes(handles.actmap);
[handles.map,~,~,handles.act_t,~,~,~,~,~,~,handles.vout,handles.quivers_Xout,handles.quivers_Yout,handles.quivers_vxout,handles.quivers_vyout,~,~,~]...
=cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),g1data.cvimages,g1data.mask,get(handles.localout,'Value'),str2double(get(g1data.minvel,'String')),str2double(get(g1data.maxvel,'String')),get(handles.actopt,'Value'),...
0,1000,str2double(get(g1data.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),1,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String'))); 



jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [map] = handles.map;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(map));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(map, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
guidata(hObject, handles);
actcurve_Callback(hObject, eventdata, handles);
local_Callback(hObject, eventdata, handles);
getMousePositionOnImage(hObject, eventdata, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns actopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from actopt


% --- Executes during object creation, after setting all properties.
function actopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in single.
function single_Callback(hObject, eventdata, handles)
% hObject    handle to single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
if get(handles.single,'value') == 1 %%need this in here to rememeber curosr points
       set(handles.actcurve,'Value',0) ;
set(handles.local,'Value',0); 
       axes(handles.actmap)
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of single


% --- Executes on button press in local.
function local_Callback(hObject, eventdata, handles)
% hObject    handle to local (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
if get(handles.local, 'Value') == 1
set(handles.actcurve,'Value',0);
set(handles.single,'Value',0);

axes(handles.actmap);

map=handles.map;
set(handles.localmean,'String',[num2str(mean(handles.vout)),'cm/s'])
set(handles.localmin,'String',[num2str(min(handles.vout)),'cm/s'])
set(handles.localmax,'String',[num2str(max(handles.vout)),'cm/s'])



jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
    [map] = handles.map;
if get(g1data.isoopt,'Value') == 1
mini=0;
maxi=max(max(map));
elseif get(g1data.isoopt,'Value') == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end
imshow(map, [0 maxi], 'InitialMagnification', 800,'colormap',jetcolormap),
caxis([mini maxi]);
scal = 0.5;
hold on
quiver(handles.quivers_Xout,handles.quivers_Yout,scal*handles.quivers_vxout,scal*handles.quivers_vyout,0,'k')
hold off
%quiversc=[vx,vy,v,Xpos,Ypos,angle]
size(handles.quivers_vxout)
size(handles.quivers_vyout)
size(handles.vout)
quiverv=[handles.quivers_vxout, handles.quivers_vyout,handles.vout,handles.quivers_Xout,handles.quivers_Yout];
for i =1:length(handles.quivers_vxout)
    vangle(i)=(atan(quiverv(i,1)/quiverv(i,2))*(180/pi));
    if quiverv(i,1) < 0
        vangle(i)=vangle(i)+180;
    end
    vangle(i)=vangle(i)+90;
end
length(handles.quivers_vxout)
vangle=vangle'
quiverv=[quiverv,vangle];
binsize=10;
hva=hist(vangle,[binsize/2:binsize:360-binsize/2])
anglevs=zeros(1,length(hva));
for i =1:length(handles.quivers_vxout)
       j=ceil(quiverv(i,6)/binsize);
       anglevs(j)=anglevs(j)+quiverv(i,3);
end
for i=1:length(hva)
    anglevs360(i)=anglevs(i)/hva(i);
end

goto180 = get(handles.opt180,'Value')

if goto180 == 1
   bin180=length(hva)/2
   for i=1:length(hva)/2 
       hva180(i)=hva(i)+hva(i+bin180)
       anglevs180(i)=anglevs(i)+anglevs(i+bin180)
       anglevs180(i)=anglevs180(i)/hva180(i);
   end
end
       

axes(handles.axes2)
cla
if get(handles.localvelopt,'Value') == 1 
histogram(handles.vout,str2double(get(g1data.binnumber,'String')));
hold on
axis tight
title('Velocity Distribution')
xlabel('Conduction Velocity (cm/s)')
ylabel('Number of Pixels')
hold off
elseif get(handles.localvelopt,'Value') == 2
if goto180 == 0
[asasasasa]=hist(vangle,[binsize/2:binsize:360-binsize/2])
plot([binsize/2:binsize:360-binsize/2],asasasasa)
    % for .csv output
    handles.angless=binsize/2:binsize:360-binsize/2;
    handles.hvals=asasasasa;
title('Angular Distribution')
xlabel('Angle of Velocity Vector')
ylabel('Number of Pixels')
xlim([0 360])
xticks([0 45 90 135 180 225 270 315 360])
end
if goto180 == 1
    plot([binsize/2:binsize:180-binsize/2],hva180)
    % for .csv output
    handles.angless=binsize/2:binsize:180-binsize/2;
    handles.hvals=hva180;
title('Angular Distribution')
xlabel('Angle of Velocity Vector')
ylabel('Number of Pixels')
xlim([0 180])
xticks([0 30 60 90 120 150 180])
end
elseif get(handles.localvelopt,'Value') == 3
    if goto180 == 0
    plot([binsize/2:binsize:360-binsize/2],anglevs360)
    handles.angless=binsize/2:binsize:360-binsize/2;
    handles.hvals=anglevs360;
xlabel('Angle')
ylabel('Conduction Velocity (cm/s)')
xlim([0 360])
xticks([0 45 90 135 180 225 270 315 360])
    end
    if goto180 == 1
        plot([binsize/2:binsize:180-binsize/2],anglevs180)
        handles.angless=binsize/2:binsize:180-binsize/2;
        handles.hvals=anglevs180;
        xlabel('Angle')
ylabel('Conduction Velocity (cm/s)')
xlim([0 180])
xticks([0 30 60 90 120 150 180])
    end
  
end
set(handles.localmean,'String',[num2str(mean(handles.vout)),'cm/s'])
set(handles.localmin,'String',[num2str(min(handles.vout)),'cm/s'])
set(handles.localmax,'String',[num2str(max(handles.vout)),'cm/s'])
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of local


% --- Executes on button press in actcurve.
function actcurve_Callback(hObject, eventdata, handles)
% hObject    handle to actcurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);

if get(handles.actcurve,'Value') ==1
set(handles.single,'Value',0);
set(handles.local,'Value',0) ;
tim=handles.act_t;
repolmin=min(tim); %save for repol calcs
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;
axes(handles.axes2);
cla
hold on
plot(xbins,tissueact,'k')
% dtis=diff(tissueact);
% dtis=dtis-min(dtis);
% dtis=dtis/(max(dtis)-min(dtis));
% plot(xbins(2:end),dtis*100,'rx')
handles.xbins=xbins;
handles.tissueact=tissueact;
xlabel('Time (ms)');
ylabel('Tissue Activated (%)');
title('Activation Curve');
ax=gca;
axis tight
ylim([0 100]);
xx=get(ax,'Xlim');
yy=get(ax,'Ylim');
hold on

actmax=str2double(get(handles.actmax,'String'));
actmin=str2double(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
Imax=Imax(1);
else Imax=allpts;
end
Imin=Imin(1);
Imin=Imin(1);

if actmax < 100
timmax=Imax*0.01;
else timmax=max(xbins);
end
tim100=max(xbins);
timmin=Imin*0.01;


closestmax=tissueact(Imax);

closestmin=tissueact(Imin);

hamin = plot([timmin timmin], [yy(1) closestmin] ,'LineStyle', '--', 'Color', 'b'); 
hbmin = plot([xx(1) timmin], [closestmin closestmin],'LineStyle', '--', 'Color', 'b');

hamax = plot([timmax timmax], [yy(1) closestmax] ,'LineStyle', '--', 'Color', 'r'); 
hbmax = plot([xx(1) timmax], [closestmax closestmax],'LineStyle', '--', 'Color', 'r');

timdiff=timmax-timmin;
if actmin == 0
    timdiff=timmax;
end
if get(handles.normalact,'Value')==1
    normfac=225/allpts;
    timdiff=timdiff*normfac;
    tim100=tim100*normfac;
end
set(handles.act100, 'String', [num2str(tim100),' ms']);
set(handles.actquote, 'String', [num2str(timdiff),' ms']);
hold off
end

 
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of actcurve


% --- Executes on button press in export_values.
function export_values_Callback(hObject, eventdata, handles)
% hObject    handle to export_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[filename,pathname] = uiputfile({'*.csv';})
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
if get(handles.actcurve,'Value') == 1
    T=table(handles.xbins',handles.tissueact');
end
if get(handles.local,'Value') == 1 && get(handles.localvelopt,'Value') == 1
    T=table(handles.vout);
end
if get(handles.local,'Value') == 1 && get(handles.localvelopt,'Value') > 1
    size(handles.angless)
    size(handles.hvals)
    T=table(handles.angless',handles.hvals');
end
if get(handles.single,'Value')==1
    T=table(handles.vangle);
end
    writetable(T,file,'Delimiter',',','WriteVariableNames',false);
% --- Executes on button press in export_graph.
function export_graph_Callback(hObject, eventdata, handles)
% hObject    handle to export_graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;ax=axes;clf;
new_handle=copyobj(handles.axes2,fig);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])

% --- Executes on button press in exportmap.
function exportmap_Callback(hObject, eventdata, handles)
% hObject    handle to exportmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;ax=axes;clf;
new_handle=copyobj(handles.actmap,fig);
jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])

function normarea_Callback(hObject, eventdata, handles)
% hObject    handle to normarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of normarea as text
%        str2double(get(hObject,'String')) returns contents of normarea as a double


% --- Executes during object creation, after setting all properties.
function normarea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in distance.
function distance_Callback(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of distance



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

function getMousePositionOnImage(hObject, eventdata, handles)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
if get(handles.singleopt,'Value') == 1
    set(handles.text12,'String','Mean CV:')
    set(handles.text13,'String','Minimum:')
    set(handles.text14,'String','Maximum:')
    set(handles.clearAB,'Visible','Off');
cursorPoint = get(handles.actmap, 'CurrentPoint');
handles.curX = cursorPoint(1,1);
handles.curY = cursorPoint(1,2);
xLimits = get(handles.actmap, 'xlim');
yLimits = get(handles.actmap, 'ylim');

if (handles.curX > min(xLimits) && handles.curX < max(xLimits) && handles.curY > min(yLimits) && handles.curY < max(yLimits))
    if get(handles.single,'value') == 1
        axes(handles.actmap)

veldistance=str2double(get(handles.veldistance,'String'));
veldistancepix=round(veldistance*10000/(str2double(get(g1data.pixelsize,'String'))));

map=handles.map;
m=min(map(map>0));
[rows cols] = size(map);
count=0;
for r=1:rows
    for c=1:cols
        if map(r,c) == m
            count=count+1;
            minind(count,1)=r;
            minind(count,2)=c;
        end
    end
end
actrow=round(mean(minind(:,1)));
actcol=round(mean(minind(:,2)));



isochoice=get(g1data.isoopt,'Value');
if isochoice == 1
mini=0;
maxi=max(max(map));
elseif isochoice == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end

imshow(map, [0 maxi], 'InitialMagnification', 800),
hold on
jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
caxis([mini maxi]);
freezeColors 
hold on
plot(actcol,actrow,'r+','MarkerSize',20,'LineWidth',3);
hold off
rrow=floor(handles.curY);
rcol=floor(handles.curX); %This way around due to matlab images being (1,1) in top left
if rrow(1) == 0
    rrow(1) = 1
end

if rcol(1) == 0
    rcol(1) = 1
end

hold on
plot(rcol,rrow,'k+','MarkerSize',20,'LineWidth',3);
th = 0:pi/180:2*pi;
xunit = veldistancepix * cos(th) + rcol;
yunit = veldistancepix * sin(th) + rrow;
h = plot(xunit, yunit,'k','LineWidth',3);
rcol
rrow
mat=[];
count=0;
xco=[100000000000];
yco=[100000000000];
for i =1:numel(xunit)

    time_A=(map(rrow,rcol));
    if  round(yunit(i)) <= rows &&  round(xunit(i)) <= cols && round(yunit(i)) > 0 &&  round(xunit(i)) > 0 
    time_B=(map(round(yunit(i)),round(xunit(i)))); %y-x switch due to image vs axis difference in matlab
    else time_B =(time_A -1) 
    end
    p1=[rcol,rrow];
    p2=[round(xunit(i)),round(yunit(i))];
    tim=time_B-time_A;
    if tim > 0
        xcn=p1(1)-p2(1)
        ycn=p1(2)-p2(2)
        if get(handles.graphhold,'Value') == 1
        if xcn ~= xco 
            if ycn ~= yco %if statment to stop repeating points with slightly diffrenct angles from th(i)
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
            end
        end
        end
        if get(handles.graphhold,'Value') == 0
                 dc=sqrt((xcn*xcn)+(ycn*ycn))
        count=count+1;
        mat(count,1)=p1(1)-p2(1); 
        mat(count,2)=p1(2)-p2(2);
        mat(count,3)=tim; %time diffrence 
        mat(count,4)=(dc*0.0001*(str2double(get(g1data.pixelsize,'String'))))/tim*1000; %speed in cm/s
        mat(count,5)=th(i)*(180/pi);
        mat(count,6)=dc; %pixel distance
        end
        
    end
end
mat=unique(mat,'rows');%get rid of repeat points
outs=get(handles.singleout,'Value');
newcount=1;
if isempty(mat) == 1
    errordlg('No postive conduction vectorrs found at this distance. Consider changing distance or starting point')
end
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


p2quick=[mat(ind_quick,1),mat(ind_quick,2)];
p2slow=[mat(ind_slow,1),mat(ind_slow,2)];
hold on

q=quiver(p1(1),p1(2),-p2quick(1),-p2quick(2),0);
q.Color='k';
q.LineWidth=3;



q=quiver(p1(1),p1(2),-p2slow(1),-p2slow(2),0); %-y due to matlab imshow and 
q.Color='r';
q.LineWidth=3;
hold off

tim_avg=mean(mat(:,3));
speed_quick=mat(ind_quick,4);
speed_slow=mat(ind_slow,4);
speed_avg=mean(mat(:,4));

zeroangle=[mat(ind_slow,5)];
vangle(:,2)=mat(:,4);
vangle(:,1)=(mat(:,5)-zeroangle);
vangle = sortrows(vangle,1)
axes(handles.axes2);
       cla reset
[~,indmin]=min(vangle(:,2));

axes(handles.axes2)       
hold on
plot(vangle(:,1),vangle(:,2),'kx','LineWidth',2,'MarkerSize',15);
plot(vangle(indmin,1),vangle(indmin,2),'rx','LineWidth',2,'MarkerSize',15);
handles.vangle=vangle;
xlabel('Angle from slowest conduction (degrees)');
ylabel('Single vector conduction velocity (cm/s)');
title('Single Vector Conduction Velocity') ;
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
set(handles.singlemean,'String',[num2str(speed_avg),'cm/s'])
set(handles.singleminimum,'String',[num2str(speed_slow),'cm/s'])
set(handles.singlemaximum,'String',[num2str(speed_quick),'cm/s'])

    end
end
end

if get(handles.singleopt,'Value') == 2
    set(handles.text12,'String','Speed:')
    set(handles.text13,'String','Distance:')
    set(handles.text14,'String','Time:')
    set(handles.clearAB,'Visible','On');
cursorPoint = get(handles.actmap, 'CurrentPoint');
handles.curX = cursorPoint(1,1);
handles.curY = cursorPoint(1,2);
xLimits = get(handles.actmap, 'xlim');
yLimits = get(handles.actmap, 'ylim');
A=handles.point_A
B=handles.point_B
if (handles.curX > min(xLimits) && handles.curX < max(xLimits) && handles.curY > min(yLimits) && handles.curY < max(yLimits))
       if get(handles.single,'value') == 1
           if isempty(handles.point_A) == 1
               handles.point_A=[round(handles.curX), round(handles.curY)]
               axes(handles.actmap);
               hold on
               plot(handles.point_A(1),handles.point_A(2),'k+','MarkerSize',20,'LineWidth',3);
           elseif isempty(handles.point_A) == 0 && isempty(handles.point_B) == 1
               handles.point_B=[round(handles.curX), round(handles.curY)]
                              axes(handles.actmap);
map=handles.map;
m=min(map(map>0));
[rows cols] = size(map);
count=0;
for r=1:rows
    for c=1:cols
        if map(r,c) == m
            count=count+1;
            minind(count,1)=r;
            minind(count,2)=c;
        end
    end
end
actrow=round(mean(minind(:,1)));
actcol=round(mean(minind(:,2)));
               
isochoice=get(g1data.isoopt,'Value');
if isochoice == 1
mini=0;
maxi=max(max(map));
elseif isochoice == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end

imshow(map, [0 maxi], 'InitialMagnification', 800),
hold on
jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
caxis([mini maxi]);
freezeColors 

               hold on
               plot(handles.point_A(1),handles.point_A(2),'k+','MarkerSize',20,'LineWidth',3);
               hold off

axes(handles.actmap)
hold on

time_A=map(handles.point_A(2),handles.point_A(1));
time_B=map(handles.point_B(2),handles.point_B(1)); %times in ms
dist=sqrt((handles.point_A(2)-handles.point_B(2))^2+(handles.point_A(1)-handles.point_B(1))^2)*str2double(get(g1data.pixelsize,'String')) %distance in um
speed=dist*0.0001/(time_B-time_A)*1000; %speed in cm/s
speed=sqrt(speed*speed)
time=sqrt((time_B-time_A)^2)
hold on
dp=handles.point_B-handles.point_A;
q=quiver(handles.point_A(1),handles.point_A(2),dp(1),dp(2),0);
q.Color='k';
q.LineWidth=3;
q.MaxHeadSize=0.5;
set(handles.singlemean,'String',[num2str(speed),' cm/s']);
set(handles.singleminimum,'String',[num2str(dist/10000),' cm']);
set(handles.singlemaximum,'String',[num2str(time),' ms']);
handles.point_B=[];
guidata(hObject, handles);
           end
       end
end
end
guidata(hObject, handles);


% --- Executes on selection change in singleout.
function singleout_Callback(hObject, eventdata, handles)
% hObject    handle to singleout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
getMousePositionOnImage(hObject, eventdata, handles)
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns singleout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from singleout


% --- Executes during object creation, after setting all properties.
function singleout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function actmin_Callback(hObject, eventdata, handles)
% hObject    handle to actmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
actcurve_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of actmin as text
%        str2double(get(hObject,'String')) returns contents of actmin as a double


% --- Executes during object creation, after setting all properties.
function actmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function actmax_Callback(hObject, eventdata, handles)
% hObject    handle to actmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
actcurve_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of actmax as text
%        str2double(get(hObject,'String')) returns contents of actmax as a double


% --- Executes during object creation, after setting all properties.
function actmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normalact.
function normalact_Callback(hObject, eventdata, handles)
% hObject    handle to normalact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
actcurve_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of normalact


% --- Executes on selection change in localout.
function localout_Callback(hObject, eventdata, handles)
% hObject    handle to localout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
[handles.map,~,~,handles.act_t,~,~,~,~,~,~,handles.vout,handles.quivers_Xout,handles.quivers_Yout,handles.quivers_vxout,handles.quivers_vyout,~,~,~]...
=cvmap(str2double(get(g1data.pixelsize,'String')),str2double(get(g1data.framerate,'String')),g1data.cvimages,g1data.mask,get(handles.localout,'Value'),str2double(get(g1data.minvel,'String')),str2double(get(g1data.maxvel,'String')),get(handles.actopt,'Value'),...
0,1000,str2double(get(g1data.winsize,'String')),str2double(get(g1data.beforeGUI,'String')),str2double(get(g1data.wint,'String')),1,str2double(get(g1data.t,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2double(get(g1data.splineN,'String'))); 
guidata(hObject, handles);
local_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns localout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from localout


% --- Executes during object creation, after setting all properties.
function localout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to localout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in graphhold.
function graphhold_Callback(hObject, eventdata, handles)
% hObject    handle to graphhold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of graphhold


% --- Executes on selection change in singleopt.
function singleopt_Callback(hObject, eventdata, handles)
% hObject    handle to singleopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns singleopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from singleopt


% --- Executes during object creation, after setting all properties.
function singleopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearAB.
function clearAB_Callback(hObject, eventdata, handles)
% hObject    handle to clearAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
map=handles.map;
axes(handles.actmap)
               
isochoice=get(g1data.isoopt,'Value');
if isochoice == 1
mini=0;
maxi=max(max(map));
elseif isochoice == 2
mini=str2double(get(g1data.isomin,'String'));
maxi=str2double(get(g1data.isomax,'String'));
end

imshow(map, [0 maxi], 'InitialMagnification', 800),
hold on
jetcolormap = (colormap('jet'));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
caxis([mini maxi]);
freezeColors 

handles.point_A=[];
set(handles.singlemean,'String',' ');
set(handles.singleminimum,'String',' ' );
set(handles.singlemaximum,'String',' ');

guidata(hObject, handles);


% --- Executes on selection change in localvelopt.
function localvelopt_Callback(hObject, eventdata, handles)
% hObject    handle to localvelopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns localvelopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from localvelopt


% --- Executes during object creation, after setting all properties.
function localvelopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to localvelopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in opt180.
function opt180_Callback(hObject, eventdata, handles)
% hObject    handle to opt180 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of opt180


% --- Executes on selection change in reopt.
function reopt_Callback(hObject, eventdata, handles)
% hObject    handle to reopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns reopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from reopt


% --- Executes during object creation, after setting all properties.
function reopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
