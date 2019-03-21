function varargout = ElectroMap(varargin)

% Main function for running ElectroMap. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Maintained by Chris O'Shea - Email CXO531@bham.ac.uk for any queries
% Version 1.0
% Release Date - 
% For license information, please see 'license.txt' at ...

% Last Updated -

% Update Summary





% ElectroMap MATLAB code for ElectroMap.fig
%      ElectroMap, by itself, creates a new ElectroMap or raises the existing
%      singleton*.
%
%      H = ElectroMap returns the handle to a new ElectroMap or the handle to
%      the existing singleton*.
%
%      ElectroMap('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ElectroMap.M with the given input arguments.
%
%      ElectroMap('Property','Value',...) creates a new ALLM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ElectroMap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      NUstop.  All inputs are passed to ElectroMap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)"
%
% See also: GUIDE, GUIDATA, GUIHANDLESF

% Edit the above text to modify the response to help ElectroMap

% Last Modified by GUIDE v2.5 19-Dec-2018 13:47:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
'gui_Singleton',  gui_Singleton, ...
'gui_OpeningFcn', @ElectroMap_OpeningFcn, ...
'gui_OutputFcn',  @ElectroMap_OutputFcn, ...
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


% --- Executes just before ElectroMap is made visible.
function ElectroMap_OpeningFcn(hObject, ~, handles, varargin)
% Fucnction for setting values on interface intialisation 
handles.output = hObject;
set(handles.invertopt,'Value',1); %inversion of signal
set(handles.sfilt,'Value',2); %spatial filtering (gaussian)
set(handles.velout, 'Value', 4); %Velociy outlier removal 

handles.bgon=1; %background on/off switch
handles.bgcol='w'; %backgroung colour

handles.folder_name=[]; %filled when directory folder chosen
handles.fname='opening of GUI hold'; 
handles.lastprocessedfname='opening of GUI hold'; %holds for loaded and processed files

handles.drawcon=0;
handles.conbon=[];
handles.medifilt=1;

masterdir=cd;
if isdeployed == 0
addpath(masterdir,'-frozen');
end

handles.ttpstart=10;
handles.ttpend=90; %time to peak defualt settings 

handles.roinum=1;
handles.roisum=0; %roi defualt settings

handles.snrt1=10;
handles.snrt2=30;

handles.pbefore=5;
handles.pafter=5;

handles.herefromroiload=0; %switch for load roi from .txt file
set(handles.manthresh,'Enable','off') %slider off unless manual trheshold level set
handles.rect=[];
handles.loadedmask=[];

set(handles.resegment,'Enable','off') 
set(handles.B2B,'Enable','off')
set(handles.pushprocess,'Enable','off') 
set(handles.producemaps,'Enable','off')
handles.herefromsegmentpush=0;
handles.filming = 0; %switch for saving maps to video files

%axis not visible until something is in them
axes(handles.mapaxes); axis off
axes(handles.bgimage);axis off
axes(handles.axes2);axis off
axes(handles.cb); axis off
axes(handles.imageaxes); axis off
zoom xon
handles.isZoomed=0;



handles.fmin=0.5;handles.fmax=10;handles.fbin=0.05;handles.dfwin=0; %Frequency mapping defaults 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ElectroMap wait for user response (see UIRESUME)
% uiwait(handles.ElectroMap);


% --- Outputs from this function are returned to the command line.
function varargout = ElectroMap_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushselect.
function pushselect_Callback(hObject, ~, ~)
% hObject    handle to pushselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
%% Populate listbox with .mat and .tif files
set(handles.listbox1,'Value',1); %set value to 1 on each new file chosen to stop error
[handles.folder_name]=uigetdir;
workingdir=handles.folder_name;
if isdeployed == 0
cd(workingdir);
end

%get all files in directory
allfiles = dir(handles.folder_name);
if isdeployed == 0
addpath(handles.folder_name);
end
file_list= {};
count=0;

%Find tif and mat files
for i=1:length(allfiles)
k = strfind(allfiles(i).name, '.TIF');
d = strfind(allfiles(i).name, '.tif');
m = strfind(allfiles(i).name, '.mat');
%l = isfolder(allfiles(i).name);
if isempty(k) ~= 1 %&& l ~= 1
count=count+1;
file=allfiles(i).name;
file_list{count}=file;
end
if isempty(d) ~= 1 %&& l ~= 1
count=count+1;
file=allfiles(i).name;
file_list{count}=file;
end
if isempty(m) ~= 1 %&& l ~= 1
count=count+1;
file=allfiles(i).name;
file_list{count}=file;
end
end
set(handles.listbox1,'String',file_list);
guidata(hObject, handles);

% --- Executes on selection change in listbox1.

function listbox1_Callback(~, ~, ~)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, ~, ~)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushload.
function pushload_Callback(hObject, ~, ~)

handles = guidata(hObject);
set(handles.listbox2,'Value',1);
set(handles.pushprocess,'Enable','on') 
set(handles.producemaps,'Enable','off')
%% Get image info from GUI
%image
handles.threshop=get(handles.threshopt,'Value'); %thershold choice
handles.threshman=(get(handles.manthresh,'Value')); %manual thershold setting 
imchoice=get(handles.imagedisp,'Value'); %image to display
cropchoice=get(handles.cropbox,'Value'); %crop image setting 
handles.cropchoice=cropchoice; %0 means no crop, 1 mean new crop, 2 means crop from before

if cropchoice == 1
handles.rect = [];
end

quinnieopt=get(handles.squareROI,'Value'); %custom ROI setting 

%% file info
chosenfilecontents=cellstr(get(handles.listbox1,'String'));
choice=get(handles.listbox1,'Value');
fname=chosenfilecontents{choice};
if ispc == 1 %change of file setting it mac or pc
handles.fnamenew=[handles.folder_name,'\',fname];
else
handles.fnamenew=[handles.folder_name,'/',fname];
end
tf = strcmp(handles.fname,handles.fnamenew);
if tf == 0
handles.rect = [];
handles.images=[];
handles.lastprocessedfname='Pressing load or changing threshold hold';
end
handles.fname=handles.fnamenew;



%% If custom roi chosen, get rid off manual thresholding
if quinnieopt == 1 || handles.herefromroiload == 1
handles.threshop = 2;
handles.threshman = -50000;
if handles.herefromroiload == 1
quinnieopt = 0;
end
end
%% load, crop and reset opt, threshold image
inversion=get(handles.invertopt,'Value');
camopt=0;
axes(handles.imageaxes)
set(handles.resegment,'Enable','off')
set(handles.B2B,'Enable','off')

%% Load new image using OMimload function
if tf == 0
[num_images,handles.newrect,mask,im,handles.I,boundaries,handles.camopt,handles.frame1,handles.fluoim,handles.rois,handles.rhsn] = OMimload(handles.fname,cropchoice,quinnieopt,handles.threshop,handles.threshman,handles.rect,inversion,camopt,get(handles.imagedisp,'Value'),handles.roinum,handles.roisum);
handles.mask=[];
handles.mask=mask;
handles.im=im; handles.num_images=num_images; handles.rect=handles.newrect;
%handles.loadedmask=[];
%handles.herefromroiload = 0;
if isempty(handles.loadedmask) == 0 && handles.herefromroiload == 1
handles.mask=[];
handles.I=[];
mask=handles.loadedmask;
mask=uint16(mask);
boundaries = bwboundaries(mask);
handles.mask=mask;
if size(im,1) ~= size(mask,1) || size(im,2) ~= size(mask,2)
if abs(size(im,1)-size(mask,1)) <= 2 && abs(size(im,2)-size(mask,2)) <= 2
choice = questdlg('ROI dimensons do not match Image but only slighty off. Would you like to reshape ROI?', ...
'ROI mismatch', ...
'Yes','No','Yes');
switch choice
case 'Yes'
[rows,cols]=size(im);
[rows2,cols2]=size(mask);


if rows<=rows2
    rows3=rows;
else
    rows3=rows2;
end

if cols<=cols2
    cols3=cols;
else
    cols3=cols2;
end

newmask=zeros(rows,cols);
for r=1:rows3
for c=1:cols3
newmask(r,c)=mask(r,c);
end
end

mask=uint16(newmask);
end
else
handles.herefromroiload = 0;
handles.loadedmask=[];
guidata(hObject,handles)
h=errordlg('Loaded ROI dimensions do not match Image');
waitfor(h)
end
end
handles.I=im.*mask;
handles.mask=mask;
end
set(handles.cropbox,'Value',0);
handles.boundaries=boundaries;
end

%% Rethreshold loaded image set 
if tf == 1
[num_images,handles.newrect,mask,im,handles.I,boundaries,handles.camopt,handles.frame1,handles.fluoim,handles.rois,handles.rhsn] = OMimload(handles.fname,cropchoice,quinnieopt,handles.threshop,handles.threshman,handles.rect,inversion,camopt,get(handles.imagedisp,'Value'),handles.roinum,handles.roisum);;
handles.mask=[];
handles.mask=mask;
handles.im=im; handles.num_images=num_images; handles.rect=handles.newrect;
%handles.loadedmask=[];
%handles.herefromroiload = 0;
if isempty(handles.loadedmask) == 0 && handles.herefromroiload == 1
handles.mask=[];
handles.I=[];
mask=handles.loadedmask;
mask=uint16(mask);

boundaries = bwboundaries(mask);
handles.mask=mask;
if size(im,1) ~= size(mask,1) || size(im,2) ~= size(mask,2)
if abs(size(im,1)-size(mask,1)) <= 2 && abs(size(im,2)-size(mask,2)) <= 2
choice = questdlg('ROI dimensons do not match Image but only slighty off. Would you like to reshape ROI?', ...
'ROI mismatch', ...
'Yes','No','Yes');
switch choice
case 'Yes'
[rows,cols]=size(im);
[rows2,cols2]=size(mask);
if rows2>rows && cols2>cols
newmask=zeros(rows,cols);
for r=1:rows
for c=1:cols
newmask(r,c)=mask(r,c);
end
end
end
if rows>rows2 && cols>cols2
newmask=zeros(rows,cols);
for r=1:rows2
for c=1:cols2
newmask(r,c)=mask(r,c);
end
end
end

mask=uint16(newmask);
end
else
handles.herefromroiload = 0;
handles.loadedmask=[];
guidata(hObject,handles);
h=errordlg('Loaded ROI dimensions do not match Image');
waitfor(h)
end
end
handles.I=im.*mask;
handles.mask=mask;
end
set(handles.cropbox,'Value',0);
handles.boundaries=boundaries;
end
%% change pic in GUI
axes(handles.imageaxes);
cla;
if imchoice == 1
imshow(handles.frame1,[],'InitialMagnification', 400)
colormap('gray')
freezeColors
hold on
for i=1:size(boundaries,1)
plot(boundaries{i}(:,2),boundaries{i}(:,1),'r','LineWidth',2);
end
hold off
end
if imchoice == 2
imshow(handles.fluoim,[],'InitialMagnification', 400)
colormap('jet')
freezeColors
hold on
for i=1:size(boundaries,1)
plot(boundaries{i}(:,2),boundaries{i}(:,1),'k','LineWidth',2);
end
hold off
end

set(handles.pushprocess,'Enable','on') 
set(handles.producemaps,'Enable','off')

guidata(hObject, handles);

% --- Executes on button press in pushprocess.
function pushprocess_Callback(hObject, ~, handles)
% hObject    handle to pushprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = guidata(hObject);
%% get processing options from GUI
%segmentation
minpeakdist = str2double(get(handles.minpeak,'String'));
minpeakdist = ceil(minpeakdist/(1/str2double(get(handles.framerate,'String'))));
handles.minpeakdist=minpeakdist;
segchoice = get(handles.segchoice,'Value');
div=str2double(get(handles.segsize,'String'));
minboundary=str2double(get(handles.minbound,'String'));
minmumofpeaks=str2double(get(handles.minnum,'String'));
handles.avgCL=[];
%Baseline
BLopt=(get(handles.BLopt,'Value'));

%filtering
tfilt=get(handles.tfilt,'Value');
sfilt=get(handles.sfilt,'Value');
sfiltsize=str2double(get(handles.sfiltsize,'String'));

%outlieropts
handles.outlier=get(handles.apdout,'Value');
handles.outliervel=get(handles.velout,'Value');

%inversion
inversion=get(handles.invertopt,'Value');

% frame removal
handles.frameremove=get(handles.removef,'Value');

%% is same file check and process images
chosenfilecontents=cellstr(get(handles.listbox1,'String'));
choice=get(handles.listbox1,'Value');
newfname=chosenfilecontents{choice};
tf = strcmp(newfname,handles.lastprocessedfname);
if tf == 0
loadnewims=1;
handles.lastprocessedfname=newfname;
elseif tf == 1 && handles.herefromsegmentpush == 0
newsettingschoice = questdlg('Re-Process?', ...
'Re-Load same file', ...
'Yes','No - Just segment','No - Just segment');
switch newsettingschoice
case 'Yes'
loadnewims = 1;
case 'No - Just segment'
loadnewims = 0;
num_images=handles.num_images;
end
elseif tf == 1 && handles.herefromsegmentpush == 1
loadnewims = 0;
end


if loadnewims == 1
if isempty(handles.rect) == 0
handles.cropchoice = 1;
end
handles.averages=[];
[handles.preimages,images,averages,mask] = OMimprocess(handles.fname,handles.im,handles.rect,handles.num_images,handles.cropchoice,handles.mask,sfilt,sfiltsize,inversion,tfilt,handles.frameremove,handles.camopt,str2double(get(handles.sfiltsigma,'String')),handles.pbefore,handles.pafter,handles.rhsn);
handles.waverages=averages;
handles.images=images;
handles.averages=averages;
num_images=handles.num_images;
handles.mask=mask;
end
set(handles.resegment,'Enable','on')
set(handles.B2B,'Enable','on')
if loadnewims == 0
averages=handles.averages;
images=handles.images;
num_images=handles.num_images;
end


%% Baseline Drift Correction

%Top hat filter 
if BLopt == 1 || BLopt == 4
th_len=str2double(get(handles.thlen,'String'));
th_len=(th_len)/str2double(get(handles.framerate,'String'));
th_len=round(th_len);

se = strel('line', th_len, 0.5);
BLAV = imopen(averages, se);
%figure, plot(BL)
end

%Poly 4th degree
if BLopt == 2 || BLopt == 5
[p,~,mu]=polyfit(1:length(averages),averages,4);
BLAV=polyval(p,1:length(averages),[],mu);
end
%Poly 11th degree
if BLopt == 3 || BLopt == 6
[p,~,mu]=polyfit(1:length(averages),averages,11);
BLAV=polyval(p,1:length(averages),[],mu);
end

% No BL correction 
if BLopt == 7
BLAV=min(averages);
end

handles.averages = (averages-BLAV); %Baseline subtraction
%% Remove baseline from each pixel

BLAV=BLAV-min(BLAV);
if BLopt == 4 || BLopt == 5 || BLopt == 6
for t = 1:size(images,3)
images(:,:,t)=images(:,:,t)+BLAV(t);
end

end

wb=waitbar(0.5,'Removing Baseline');

signal=zeros(1,size(images,3));
if BLopt == 1 || BLopt == 2 || BLopt == 3
for row=1:size(images,1) %%MY BL REMOVAL
for col=1:size(images,2)
for frame = 1:size(images,3)
signal(frame)=images(row,col,frame);
end
if inversion == 1
signal=imcomplement(signal);
end
if BLopt == 1
se = strel('line', th_len, 0.5);
BL = imopen(signal, se);
end

if BLopt == 2
[p,~,mu]=polyfit(1:length(signal),signal,4);
BL=polyval(p,1:length(images(row,col,:)),[],mu);
end

if BLopt == 3
[p,~,mu]=polyfit(1:length(signal),signal,11);
BL=polyval(p,1:length(signal),[],mu);
end


for frame = 1:size(images,3)
images(row,col,frame)=images(row,col,frame)+BL(frame);
end
images(row,col,:)=images(row,col,:) - min(images(row,col,:)); %make all mins zero
end

end

end




handles.images=images;
waitbar(0.95,wb,'Segmenting Signal');
wholav=handles.averages;

%% Display Signal
schoice = 1
if schoice == 1
handles.averages=handles.waverages;
end
if schoice == 2 || schoice == 3 || schoice == 4
if schoice == 2 || schoice == 4
figure,
imshow(handles.frame1, [],'InitialMagnification', 800) 
title('Make your selection and press enter');
[~,rec]=imcrop;
cropfig=gcf;
close(cropfig)
rec=floor(rec);
r1=rec(2);
c1=rec(1);
if r1 == 0 
r1=1;
end
if c1 == 0 
c1=1;
end
r2=floor(rec(2)+rec(4));
c2=floor(rec(1)+rec(3));
[ar,ac,numim]=size(handles.images);
rmask=zeros(ar,ac);

for r=r1:r2
for c=c1:c2
rmask(r,c)=1;   
end
end
rmask=uint16(rmask);
newav=zeros(1,numim);
for j=1:numim
class(handles.images)
class(rmask)
roiim=handles.images(:,:,j).*rmask;
newav(j)=sum(sum(roiim));
end

% figure,
newav=imcomplement(newav);
newav=newav-min(newav);
end
if schoice == 3 || schoice == 4
if schoice == 3
newav=handles.waverages;
end
dnewav=smooth(newav);
dnewav=smooth(diff(dnewav));
dnewav(1:10)=0;
dnewav=dnewav-min(dnewav);
dnewav(1:10)=0;
newav=[0,dnewav'];
end
%save overall average for later
wholav=handles.averages;
handles.averages=newav;
end
set(handles.listbox2,'Value',1)

axes(handles.axes2)
plot(handles.averages)
drawnow()

%% BLremoval
%% Baseline Drift Correction
if schoice == 1 || schoice == 2
if schoice == 1
handles.averages=wholav;
end
if BLopt == 1 || BLopt == 4
th_len=str2double(get(handles.thlen,'String'));
th_len=(th_len)/str2double(get(handles.framerate,'String'));
th_len=round(th_len);

se = strel('line', th_len, 0.5);
BLAV = imopen(handles.averages, se);
%figure, plot(BL)
end

if BLopt == 2 || BLopt == 5
[p,~,mu]=polyfit(1:length(handles.averages),handles.averages,4);
BLAV=polyval(p,1:length(handles.averages),[],mu);
end

if BLopt == 3 || BLopt == 6
[p,~,mu]=polyfit(1:length(handles.averages),handles.averages,11);
BLAV=polyval(p,1:length(handles.averages),[],mu);
end

if BLopt == 7
BLAV=min(handles.averages);
end

handles.averages = (handles.averages-BLAV);
end
axes(handles.axes2)
plot(handles.averages)
drawnow()

%% DETECT PEAKS
before=str2double(get(handles.beforeGUI,'String'));
before=before*str2double(get(handles.framerate,'String'));
before=round(before);
after=str2double(get(handles.afterGUI,'String'));
after=after*str2double(get(handles.framerate,'String'));
after=round(after);
handles.locs=[];handles.q2locs=[];handles.avgCL=[];
[handles.locs,~,handles.q2locs,handles.avgCL,handles.numofpeaksoverall,handles.peakheight]=Omseg2...
(handles.averages,str2double(get(handles.peakhigh,'String')),minpeakdist,str2double(get(handles.peakhigh,'String')),minpeakdist,minboundary,segchoice,minmumofpeaks,num_images,div,before,after);

%% Zoomed Section
axes(handles.axes2);
origInfo = getappdata(gca, 'matlab_graphics_resetplotview');
handles.isZoomed = 0;
if isempty(origInfo)
handles.isZoomed = 0;
else
handles.isZoomed = 1;
end
exposure=1/str2double(get(handles.framerate,'String'));
handles.newlim=get(gca,'XLim')/exposure;
handles.newlim(1)=floor(handles.newlim(1));
handles.newlim(2)=ceil(handles.newlim(2));
%% axes2

handles.avgCL=handles.avgCL.*(1/str2double(get(handles.framerate,'String')));
axes(handles.axes2)
cla
CM =['b','r','g','y','c','m','k'];
exposure=1/str2double(get(handles.framerate,'String'));
handles.averagestime=(0:1:(length(handles.averages)-1))*exposure;
plot(handles.averagestime, handles.averages,'k'),
xlabel('time (ms) \rightarrow');
ylabel('Fluorescence Intensity');
xlim([0 length(handles.averages)*exposure]);
hold on
plot(handles.averagestime(handles.locs),handles.averages(handles.locs), 'or');
before=round(str2double(get(handles.beforeGUI,'String'))/exposure);
after=round(str2double(get(handles.afterGUI,'String'))/exposure);
if length(handles.locs) < 2
 handles.q2locs=handles.locs;
 handles.avgCL(2,1)=0;
end
for i = 1:length(handles.q2locs(:,1))
c=mod(i,6);
if c == 0
c=6;
end
handles.q2locs;
A=(handles.q2locs(i,:));
if isempty(A) == 1
errordlg('No constant cycle length regions found. Please adjust pre-process settings')
end
if min(A(A>0)) < before %if first peak v.close to beginning it is ignored to stop dim error
k = find(A);
A(k(1))=0;
end
tstart=min(A(A>0))-before;
if tstart == 0
tstart = 1;
end
tend=max(A)+after;
if tend > length(handles.averagestime)
if length(A)>1
tend = A(end-1)+after;
elseif length(A) == 1
tend=length(handles.averagestime);
c=7;
end
end

if tend>length(handles.averagestime)
tend=length(handles.averagestime);
end

plot(handles.averagestime(tstart:tend),handles.averages(tstart:tend),'color',CM(c));
end
%end
colorbar off
hold off

set(gca,'FontSize',8);
ax=gca;
line(get(ax,'XLim'),[handles.peakheight handles.peakheight],'Color','b')
%populate section listbox

section = {};

if handles.numofpeaksoverall == 1
section{1}=('N/A');
else
for i = 1:length(handles.q2locs(:,1))
length(handles.q2locs(:,1));
handles.q2locs;
handles.avgCL;
section{i}=[num2str(i),' (',num2str((handles.avgCL(2,i))),'ms)'];
end
end

%% add zoomed section
if handles.isZoomed == 1
handles.q2locs
if length(handles.q2locs(1,:)) < 2 %for single peak seg
handles.q2locs(:,end+1)=0;
end
newline=zeros(1,length(handles.q2locs(1,:)));
newline(1)=handles.newlim(1);
newline(2)=handles.newlim(2);
handles.q2locs=[handles.q2locs;newline];
newsection='Zoomed Section';
section{length(section)+1}=newsection;
end

handles.section=section;
%% Couple sections to auto windows
%peak count set to 1 because of first ignored peak, should be changes
handles.winopt=1;


if length(handles.locs) == 1 || length(handles.locs) == 2
if length(handles.locs) == 1
handles.q2locs(1,1)=handles.locs;
end
if length(handles.locs) == 2
handles.q2locs=[];
handles.q2locs=handles.locs;
CL2=handles.locs(2)-handles.locs(1)*exposure;
handles.section;
handles.section{1}=[num2str(CL2),'ms'];
section=handles.section;
end
end
set(handles.listbox2,'String',section);
axes(handles.mapaxes);
pretty=get(handles.colmap,'String');
jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
delete(wb);

if schoice == 2 || schoice == 3
handles.averages=wholav;
end
set(handles.pushprocess,'Enable','on') 
set(handles.producemaps,'Enable','on')
guidata(hObject, handles);

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, ~)
handles = guidata(hObject);
axes(handles.axes2)
%handles.filming = 0;
exposure=1/str2double(get(handles.framerate,'String'));
handles.averagestime=(0:1:(length(handles.averages(1,:))-1))*exposure;
plot(handles.averagestime, handles.averages,'k'),
xlabel('time (ms) \rightarrow');
ylabel('Fluorescence Intensity');
xlim([0 length(handles.averages)*exposure]);
hold on
plot(handles.averagestime(handles.locs),handles.averages(handles.locs), 'or');
before=round(str2double(get(handles.beforeGUI,'String'))/exposure);
after=round(str2double(get(handles.afterGUI,'String'))/exposure);
section_choice=get(handles.listbox2,'Value');
A=(handles.q2locs(section_choice,:));

if isempty(A) == 1
errordlg('No constant cycle length regions found. Please adjust pre-process settings')
end
if min(A(A>0)) < before %if first peak v.close to beginning it is ignored to stop dim error
k = find(A);
A(k(1))=0;
end
tstart=min(A(A>0))-before;
if tstart == 0
tstart = 1;
end

tend=max(A)+after;
if tend > length(handles.averagestime)
tend = length(handles.averagestime);
end
plot(handles.averagestime(tstart:tend),handles.averages(tstart:tend),'color','r');
ax=gca;
line(get(ax,'XLim'),[handles.peakheight handles.peakheight],'Color','b')
set(gca,'FontSize',8);


producemaps_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, ~, ~)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tfilt.
function segchoice_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function segchoice_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function segsize_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function segsize_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BLopt.
function BLopt_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function BLopt_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tfilt.
function tfilt_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function tfilt_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function minpeak_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function minpeak_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function minnum_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function minnum_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function minbound_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function minbound_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in invertopt.
function invertopt_Callback(~, ~, ~)
% hObject    handle to invertopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertopt


% --- Executes on selection change in sfilt.
function sfilt_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function sfilt_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function sfiltsize_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function sfiltsize_CreateFcn(hObject, ~, ~)
% hObject    handle to sfiltsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Mapchoice.
function Mapchoice_Callback(hObject, ~, ~)
% hObject    handle to Mapchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
choice=get(handles.Mapchoice,'Value');
axes(handles.bgimage);
if handles.bgon == 0
imshow(handles.frame1,[],'InitialMagnification', 400)
colormap('gray')
freezeColors
else
cla reset
axis off
end
[rows,cols]=size(handles.frame1)
axes(handles.mapaxes);
cla reset
axis off
drawnow()
imshow(zeros(rows,cols))
isochoice=get(handles.isoopt,'Value');
CVmap=[];
if choice == 10
cla reset
set(handles.meanchange,'String','');
set(handles.textchange,'String','');
[map,~,alll,~]=snrs(handles.averageBeat,handles.mask,handles.snrt1,handles.snrt2,(get(handles.tfilt,'Value')));
if isempty(map) == 1
map=0;
alll=0;
end
dmap=map; %displau map. values altered so displyaed so bg=white etc. Display and expoterd stats NOT based on this map. 
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
colormap(jetcolormap);
if get(handles.apdscale,'Value') == 1
dmap(dmap==0)= NaN;
him=imshow(dmap,[], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([floor(min(min(alll))) ceil(max(max(alll)))])
end
if get(handles.apdscale,'Value') == 2
dmap(dmap==0)=NaN;
him=imshow(dmap,[str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))])
end
title('SNR Map')
freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
pretty=get(handles.colmap,'String');
colormap(colormap(pretty{get(handles.colmap,'Value')}));
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
if get(handles.apdscale,'Value') == 1
stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
hcb.TickLabels=(floor(min(min(alll))):stepp:ceil(max(max(alll))));
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Signal/Noise';
axis off
end
if get(handles.apdscale,'Value') == 2
stepp=(str2double(get(handles.cmax,'String'))-str2double(get(handles.cmin,'String')))/5;
hcb.TickLabels=[str2double(get(handles.cmin,'String')):stepp:str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Signal/Noise';
axis off
end
pos = get(hcb,'Position');
hcb.Label.Position=[pos(1) pos(2)-1.2];
stdall=std(alll);
palll=prctile(alll,[5,50,95]);
handles.rdata(4,1)=mean(alll);handles.rdata(4,2)=stdall;handles.rdata(4,3)=stdall/sqrt(numel(alll));handles.rdata(4,4)=stdall*stdall;handles.rdata(4,5)=((palll(3)-palll(1))/palll(2));

rownames{1}='APD';rownames{2}='CV';rownames{3}='Amp';rownames{4}='SNR';
axes(handles.mapaxes);
set(handles.rtable,'RowName',rownames);
set(handles.rtable,'data',handles.rdata);
end
if choice == 1 
cla reset
set(handles.meanchange,'String','');
set(handles.textchange,'String','');
t=str2double(get(handles.t,'String'));

map=handles.apdmap;
alll=handles.apalll;
if isempty(map) == 1
map=0;
alll=0;
end
dmap=map; %displau map. values altered so displyaed so bg=white etc. Display and expoterd stats NOT based on this map. 
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));

colormap(jetcolormap);
if get(handles.apdscale,'Value') == 1
dmap(dmap==0)= NaN;
him=imshow(dmap,[], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([floor(min(min(alll))) ceil(max(max(alll)))])
end
if get(handles.apdscale,'Value') == 2
dmap(dmap==0)=NaN;
him=imshow(dmap,[str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))])
end
freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
pretty=get(handles.colmap,'String');
colormap(colormap(pretty{get(handles.colmap,'Value')}));
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
if get(handles.apdscale,'Value') == 1
stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Duration (ms)';
axis off
end
if get(handles.apdscale,'Value') == 2
stepp=(str2double(get(handles.cmax,'String'))-str2double(get(handles.cmin,'String')))/5;
hcb.TickLabels=[str2double(get(handles.cmin,'String')):stepp:str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Duration (ms)';
axis off
end
pos = get(hcb,'Position');
hcb.Label.Position=[pos(1) pos(2)-1.2];
stdall=std(alll);
palll=prctile(alll,[5,50,95]);

handles.rdata(1,1)=mean(alll);handles.rdata(1,2)=stdall;handles.rdata(1,3)=stdall/sqrt(numel(alll));handles.rdata(1,4)=stdall*stdall;handles.rdata(1,5)=((palll(3)-palll(1))/palll(2));

[~,~,allSNRr,allSNRdb]=snrs(handles.averageBeat,handles.mask,handles.snrt1,handles.snrt2,(get(handles.tfilt,'Value')));

handles.rdata(4,1)=mean(allSNRr);handles.rdata(4,2)=mean(allSNRdb);
rownames=get(handles.rtable,'RowName');
rownames{4}='SNR';
axes(handles.mapaxes);
title(['APD', num2str(t)]);
set(handles.rtable,'RowName',rownames);
set(handles.rtable,'data',handles.rdata);
end

if choice == 2
cla
set(handles.meanchange,'String','');
set(handles.textchange,'String','');
map=handles.actmap;
if isochoice == 1
mini=0;
maxi=max(max(map));
elseif isochoice == 2
mini=str2double(get(handles.isomin,'String'));
maxi=str2double(get(handles.isomax,'String'));
end
dmap=map;
dmap(dmap==0)= NaN;
him=imshow(dmap, [0 maxi], 'InitialMagnification', 800);
hold on
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
colormap(jetcolormap);
caxis([mini maxi]);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
title('Activation Map');


freezeColors

%colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
stepp=(maxi-mini)/5;
hcb.TickLabels=[mini:stepp:maxi];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Time of activation (ms)';
axis off
handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
rownames=get(handles.rtable,'RowName');
rownames{4}='';
set(handles.rtable,'RowName',rownames);
set(handles.meanchange,'String',rownames{4});
set(handles.rtable,'data',handles.rdata);
end
if choice == 3
cla
set(handles.meanchange,'String','');
set(handles.textchange,'String','');

map=handles.actmap;
quivers_X=handles.quivers_X;
quivers_Y=handles.quivers_Y;
quivers_vx=handles.quivers_vx;
quivers_vy=handles.quivers_vy;

if isochoice == 1
mini=0;
maxi=max(max(map));
elseif isochoice == 2
mini=str2double(get(handles.isomin,'String'));
maxi=str2double(get(handles.isomax,'String'));
end
dmap=map;
dmap(dmap==0)= NaN;
him=imshow(dmap, [0 maxi], 'InitialMagnification', 800);
hold on
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
colormap(jetcolormap);
caxis([mini maxi]);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
scal = str2double(get(handles.scal,'String'));
hold on
quiver(quivers_X,quivers_Y,scal*quivers_vx,scal*quivers_vy,0,'k')

% Compile CV map
% factor=str2double(get(handles.pixelsize,'String'))/10;
% [rows,cols] = size(map);
% CVXmap=zeros(rows,cols);
% CVXmap(sub2ind(size(CVXmap),quivers_Y,quivers_X)) = quivers_vx;
% sCVXmap=CVXmap.*CVXmap;
% CVYmap=zeros(rows,cols);
% CVYmap(sub2ind(size(CVYmap),quivers_Y,quivers_X)) = quivers_vy;
% sCVYmap=CVYmap.*CVYmap;
% %construct cv map
% CVmap=sqrt(sCVXmap+sCVYmap);
% CVmap=CVmap*factor;

hold off
title('Activation Map');
freezeColors
%colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual'; hcb.TickLabelsMode='manual';
stepp=(maxi-mini)/5;
hcb.TickLabels=[mini:stepp:maxi];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Time of activation (ms)';
axis off
handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
rownames=get(handles.rtable,'RowName');
rownames{4}='';
set(handles.rtable,'RowName',rownames);
set(handles.meanchange,'String',rownames{4});
set(handles.rtable,'data',handles.rdata);
map=CVmap;
end
if choice == 4
cla
set(handles.meanchange,'String','');
set(handles.textchange,'String','');
map=handles.actmap;
quivers_Xout=handles.quivers_Xout;
quivers_Yout=handles.quivers_Yout;
quivers_vxout=handles.quivers_vxout;
quivers_vyout=handles.quivers_vyout;

if isochoice == 1
mini=0;
maxi=max(max(map));
elseif isochoice == 2
mini=str2double(get(handles.isomin,'String'));
maxi=str2double(get(handles.isomax,'String'));
end
axes(handles.mapaxes)
dmap=map;
dmap(dmap==0)= NaN;
him=imshow(dmap, [0 maxi], 'InitialMagnification', 800);
hold on
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
colormap(jetcolormap);
caxis([mini maxi]);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
scal = str2double(get(handles.scal,'String'));
hold on
quiver(quivers_Xout,quivers_Yout,scal*quivers_vxout,scal*quivers_vyout,0,'k')
hold off
% Compile CV map
% factor=str2double(get(handles.pixelsize,'String'))/10;
% [rows cols] = size(map);
% CVmap=zeros(rows,cols);
% CVXmap=zeros(rows,cols);
% CVXmap(sub2ind(size(CVXmap),quivers_Yout,quivers_Xout)) = quivers_vxout;
% sCVXmap=CVXmap.*CVXmap;
% CVYmap=zeros(rows,cols);
% CVYmap(sub2ind(size(CVYmap),quivers_Yout,quivers_Xout)) = quivers_vyout;
% sCVYmap=CVYmap.*CVYmap;
% %construct cv map
% CVmap=sqrt(sCVXmap+sCVYmap);
% CVmap=CVmap*factor;
title('Activation Map');

freezeColors

axes(handles.cb);
cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
stepp=(maxi-mini)/5;
hcb.TickLabels=[mini:stepp:maxi];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Time of activation (ms)';
axis off
handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
rownames=get(handles.rtable,'RowName');
rownames{4}='';
set(handles.rtable,'RowName',rownames);
set(handles.meanchange,'String',rownames{4});
set(handles.rtable,'data',handles.rdata);
map=CVmap;
end

if choice == 5
cla
axes(handles.mapaxes);
axis off
%wb=waitbar(0.9,'Calculating Frequencies');
tic
[map]=domfreq(handles.mask,handles.imagerange,str2double(get(handles.framerate,'String')),handles.fmin,handles.fmax,handles.fbin,handles.dfwin,get(handles.tfilt,'Value'));
toc
dfs=map(map>0);
dmap=map;
if get(handles.apdscale,'Value') == 1
dmap(dmap==0)= NaN;
him=imshow(dmap,[min(min(dfs)) max(max(dfs))], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([floor(min(min(dfs))) ceil(max(max(dfs)))])
title('Dominant Frequency Map')
axes(handles.cb);

cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
hcb.TickLabels=[min(min(dfs)) max(max(dfs))];
hcb.Ticks=[0 1];
hcb.Label.String='Dominant Frequency (Hz)';
end


pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
colormap(jetcolormap);

if get(handles.apdscale,'Value') == 2
dmap(dmap==0)= NaN;
him=imshow(dmap,[str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))])
axes(handles.cb);

cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
hcb.TickLabels=[str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0 1];
hcb.Label.String='Dominant Frequency (Hz)';
end
freezeColors



axis off
dfs=map(map>0);
dfsp=prctile(dfs,[5,50,95]);
handles.rdata(4,1)=mean(dfs);handles.rdata(4,2)=std(dfs);handles.rdata(4,3)=std(dfs)/numel(dfs);handles.rdata(4,4)=std(dfs)*std(dfs);handles.rdata(4,4)=std(dfs)*std(dfs);handles.rdata(4,5)=((dfsp(3)-dfsp(1))/dfsp(2));
rownames=get(handles.rtable,'RowName');
rownames{4}='DF';
set(handles.rtable,'RowName',rownames);
set(handles.rtable,'data',handles.rdata);
end
if choice == 6
wb=waitbar(0.5,'Producing Diastolic Map');
section_choice=get(handles.listbox2,'Value');
A=handles.q2locs(section_choice,:);
A=A(A~=0);
frame_1=A(1);
frame_last=A(end);
exposure = 1/str2double(get(handles.framerate,'String'));
after=str2double(get(handles.afterGUI,'String'));
after=round(after/exposure);

if frame_last+after>handles.num_images
frame_last=A(end-1);
end
[map,~,alll]=DInt(str2double(get(handles.framerate,'String')),handles.I,handles.images,handles.outlier,str2double(get(handles.cmin,'String')),str2double(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.afterGUI,'String')),str2double(get(handles.minpeak,'String')),frame_1,frame_last,str2double(get(handles.t,'String')));
if isempty(alll) == 1
map=0;
alll=0;
end
numel(alll)
map(isnan(map)) = 0;
axes(handles.mapaxes);
imshow(map,[], 'InitialMagnification', 800);
title(['Diastolic Interval Distribution']);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
if get(handles.apdscale,'Value') == 1
caxis([floor(min(min(alll))) ceil(max(max(alll)))])
end
if get(handles.apdscale,'Value') == 2
caxis([str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))])
end
freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;

hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
if get(handles.apdscale,'Value') == 1

stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
hcb.TickLabels=(floor(min(min(alll))):stepp:ceil(max(max(alll))));
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Diastolic Interval (ms)';
axis off
end
if get(handles.apdscale,'Value') == 2
stepp=(str2double(get(handles.cmax,'String'))-str2double(get(handles.cmin,'String')))/5;
hcb.TickLabels=[str2double(get(handles.cmin,'String')):stepp:str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Diastolic Interval (ms)';
axis off

end
dis=alll;
disper=prctile(dis,[5,50,95]);
handles.rdata(4,1)=mean(dis);handles.rdata(4,2)=std(dis);handles.rdata(4,3)=std(dis)/numel(dis);handles.rdata(4,4)=std(dis)*std(dis);handles.rdata(4,4)=std(dis)*std(dis);handles.rdata(4,5)=((disper(3)-disper(1))/disper(2));
rownames=get(handles.rtable,'RowName');
rownames{4}='DI';
set(handles.rtable,'RowName',rownames);
set(handles.rtable,'data',handles.rdata);
delete(wb)
end
if choice == 7
cla
t=str2double(get(handles.t,'String'));

[map,~,alll]=ttpnew(handles.ttpstart,handles.ttpend,str2double(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2double(get(handles.cmin,'String')),str2double(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2double(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2double(get(handles.apdblnum,'String')));
dmap=map;
title('TTP');
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
colormap(jetcolormap);
dmap(dmap==0)= NaN;
if get(handles.apdscale,'Value') == 1
dmap(dmap==0)= NaN;
him=imshow(dmap,[], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([floor(min(min(alll))) ceil(max(max(alll)))])
end
if get(handles.apdscale,'Value') == 2
dmap(dmap==0)=NaN;
him=imshow(dmap,[str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))], 'InitialMagnification', 800);
set(him, 'AlphaData', ~isnan(dmap))
if handles.bgon == 1
axis on
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
end
caxis([str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))])
end
title('Time to peak Map');
freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
if get(handles.apdscale,'Value') == 1
max(max(alll))
min(min(alll))
stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Duration (ms)';
axis off
end
if get(handles.apdscale,'Value') == 2
stepp=(str2double(get(handles.cmax,'String'))-str2double(get(handles.cmin,'String')))/5;
hcb.TickLabels=[str2double(get(handles.cmin,'String')):stepp:str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Time to Peak (ms)';
axis off
end
ttps=alll;
ttpsp=prctile(ttps,[5,50,95]);
handles.rdata(4,1)=mean(ttps);handles.rdata(4,2)=std(ttps);handles.rdata(4,3)=std(ttps)/numel(ttps);handles.rdata(4,4)=std(ttps)*std(ttps);handles.rdata(4,4)=std(ttps)*std(ttps);handles.rdata(4,5)=((ttpsp(3)-ttpsp(1))/ttpsp(2));
rownames=get(handles.rtable,'RowName');
rownames{4}='TTP';
set(handles.rtable,'RowName',rownames);
set(handles.rtable,'data',handles.rdata);


end

if choice == 8
cla
wb=waitbar(0.5,'Producing tau map');
[map,~,alll,~,~,~,~]=tautest3(str2double(get(handles.taustart,'String')),str2double(get(handles.taufinish,'String')),str2double(get(handles.framerate,'String')),handles.I,handles.images,handles.averageBeat,handles.outlier,str2double(get(handles.cmin,'String')),str2double(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2double(get(handles.r2cut,'String')));
map(isnan(map)) = 0;
cla
axes(handles.mapaxes);
imshow(map,[], 'InitialMagnification', 800);
title('Relaxation Constant Map');
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
if get(handles.apdscale,'Value') == 1
caxis([floor(min(min(alll))) ceil(max(max(alll)))])
end
if get(handles.apdscale,'Value') == 2
caxis([str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))])
end
freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
if get(handles.apdscale,'Value') == 1
max(max(alll));
min(min(alll));
stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Relaxation Constant (ms)';
axis off
end
if get(handles.apdscale,'Value') == 2
stepp=(str2double(get(handles.cmax,'String'))-str2double(get(handles.cmin,'String')))/5;
hcb.TickLabels=[str2double(get(handles.cmin,'String')):stepp:str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Relaxation Constant (ms)';
axis off
end
delete(wb)
taus=alll;
tausp=prctile(taus,[5,50,95]);
handles.rdata(4,1)=mean(taus);handles.rdata(4,2)=std(taus);handles.rdata(4,3)=std(taus)/numel(taus);handles.rdata(4,4)=std(taus)*std(taus);handles.rdata(4,4)=std(taus)*std(taus);handles.rdata(4,5)=((tausp(3)-tausp(1))/tausp(2));
rownames=get(handles.rtable,'RowName');
rownames{4}='Tau';
set(handles.rtable,'RowName',rownames);
set(handles.rtable,'data',handles.rdata);

end
if choice == 9
[map]=fluo_map(str2double(get(handles.framerate,'String')),handles.I,handles.images,get(handles.tfilt,'Value'),handles.averageBeat);
map(isnan(map)) = 0;
alll=map(map>0);
cla
axes(handles.mapaxes)
if get(handles.apdscale,'Value') == 1
imshow(map,[floor(min(min(alll))) ceil(max(max(alll)))], 'InitialMagnification', 800);
end
if get(handles.apdscale,'Value') == 2
imshow(map,[str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))], 'InitialMagnification', 800);
end
title('Amplitude Map');
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);

freezeColors
% colorbar
axes(handles.cb);
cla reset
hcb=colorbar;
hcb.Location = 'southoutside';
cpos = hcb.Position;
cpos(4) = 4*cpos(4);
hcb.Position = cpos;
if get(handles.apdscale,'Value') == 1
stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Relaxation Constant (ms)';
axis off
end
if get(handles.apdscale,'Value') == 2
stepp=(str2double(get(handles.cmax,'String'))-str2double(get(handles.cmin,'String')))/5;
hcb.TickLabels=[str2double(get(handles.cmin,'String')):stepp:str2double(get(handles.cmax,'String'))];
hcb.Ticks=[0.01,0.2:0.2:1];
hcb.Label.String='Relaxation Constant (ms)';
axis off
end
hcb.Label.String='Signal Level';
axis off
handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
rownames=get(handles.rtable,'RowName');
rownames{4}='';
set(handles.meanchange,'String',rownames{4});
set(handles.rtable,'data',handles.rdata);
end



% Draw Contours
drawcon = handles.drawcon;
if drawcon == 1
axes(handles.mapaxes);
hold on
for j = 0:str2double(handles.conbon):max(max(map))+str2double(handles.conbon)
mapmask=(map<=j);
A=map.*mapmask;
[cons]=bwboundaries(A);
for i=1:size(cons,1)
plot(cons{i}(:,2),cons{i}(:,1),'k','LineWidth',1);
end
end
end

%    handles.mask=handles.hold_mask; %keep and reinsate overall mask at end
handles.holdmap=map;
handles.holdcvmap=CVmap;
drawnow()
guidata(hObject, handles);
drawnow()

% Hints: contents = cellstr(get(hObject,'String')) returns apdvaluechoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdvaluechoice


% --- Executes during object creation, after setting all properties.
function Mapchoice_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mapoptions.
function mapoptions_Callback(~, ~, ~)
% hObject    handle to mapoptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mapoptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mapoptions


% --- Executes during object creation, after setting all properties.
function mapoptions_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportMap.
function ExportMap_Callback(~, ~, handles)

GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;ax=axes;clf;
new_handle=copyobj(handles.mapaxes,fig);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])


% --- Executes on button press in actpoints.
function actpoints_Callback(hObject, ~, ~)
% hObject    handle to actpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if get(handles.actfittimes,'Value') == 1
[~,act_x,act_y,act_t,~,~,~,~,~,~,~,~,~, ~,~,~,~,~] =...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
0,100,str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
[~,act_x,act_y,act_t,~,~,~,~,~,~,~,~, ~,~,~,~,~] =...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
str2double(get(handles.MINt,'String')),str2double(get(handles.MAXt,'String')),str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
figure, hold on
plot3(act_x,act_y,act_t,'.k')
title('activation points');
zlabel('time (ms)', 'FontSize', 20);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
zlim([0 15])
hold off

% --- Executes on button press in velhist.
function velhist_Callback(hObject, ~, ~)
% hObject    handle to velhist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if get(handles.actfittimes,'Value') == 1
[~,~,~,~,~,~,~,~,~,~,vout,~,~,~,~,~,~,~] =...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
0,100,str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
[~,~,~,~,~,~,~,~,~,~,vout,~,~,~,~,~,~,~] =...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
str2double(get(handles.MINt,'String')),str2double(get(handles.MAXt,'String')),str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
figure, histogram(vout,str2double(get(handles.binnumber,'String')));
hold on
%line([mean(vout) mean(vout)],[0 300],'Color','r','Linewidth',3)
xlabel('Conduction Velocity (cm/s)')
ylabel('Number of Pixels')
hold off
% --- Executes on button press in APDdist.
function APDdist_Callback(hObject, ~, ~)
% hObject    handle to APDdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
alll=handles.apalll;
figure, histogram(alll,str2double(get(handles.binnumber,'String')));
hold on
xlabel('Action Potential Duration (ms)')
ylabel('Number of Pixels')
hold off

% --- Executes on selection change in apdout.
function apdout_Callback(hObject, eventdata, ~)
% hObject    handle to apdout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function apdout_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function cmin_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function cmin_CreateFcn(hObject, eventdata, handles)


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function cmax_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function cmax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in velout.
function velout_Callback(~, ~, handles)
% hObject    handle to velout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.velout, 'Value') == 9
set(handles.minvel,'String','0.2')
set(handles.maxvel,'String','0.8')
set(handles.text64,'String','(0-1)');
elseif get(handles.velout, 'Value') == 2 
set(handles.minvel,'String','0')
set(handles.maxvel,'String','100')
set(handles.text64,'String','cm/s');
end
% Hints: contents = cellstr(get(hObject,'String')) returns velout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velout


% --- Executes during object creation, after setting all properties.
function velout_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function minvel_Callback(~, ~, ~)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minvel as text
%        str2double(get(hObject,'String')) returns contents of minvel as a double


% --- Executes during object creation, after setting all properties.
function minvel_CreateFcn(hObject, ~, ~)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function maxvel_Callback(~, ~, ~)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxvel as text
%        str2double(get(hObject,'String')) returns contents of maxvel as a double


% --- Executes during object creation, after setting all properties.
function maxvel_CreateFcn(hObject, ~, ~)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushrecal.
function pushrecal_Callback(hObject, eventdata, ~)
% hObject    handle to pushrecal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
axes(handles.mapaxes);
handles.outlier=get(handles.apdout,'Value');
handles.outliervel=get(handles.velout,'Value');
wb=waitbar(0.4,'Calculating APD');
t=str2double(get(handles.t,'String'));
                            
[apmap,meanapd,alll,onedev]=mapsbaby(get(handles.aptime1,'Value'),str2double(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2double(get(handles.cmin,'String')),str2double(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2double(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2double(get(handles.apdblnum,'String')),handles.medifilt);
stdall=std(alll);
palll=prctile(alll,[5,50,95]);
handles.rdata(4,1)=mean(alll);handles.rdata(4,2)=stdall;handles.rdata(4,3)=stdall/sqrt(numel(alll));handles.rdata(4,4)=stdall*stdall;handles.rdata(4,5)=((palll(3)-palll(1))/palll(2));
waitbar(0.6,wb,'Calculating conduction velocity');
if get(handles.actfittimes,'Value') == 1
[actmap,~,~,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,~,~,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevcv,varicv,SEcv] =...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
0,100,str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
[actmap,~,~,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,~,~,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevcv,varicv,SEcv] =...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
str2double(get(handles.MINt,'String')),str2double(get(handles.MAXt,'String')),str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
%save AP,CV maps and vectors
handles.apdmap=apmap;
handles.apalll=alll;

handles.actmap=actmap;
handles.act_t=act_t;
handles.vout=vout;
handles.quivers_X=quivers_X;
handles.quivers_Y=quivers_Y;
handles.quivers_vx=quivers_vx;
handles.quivers_vy=quivers_vy;
handles.quivers_Xout=quivers_Xout;
handles.quivers_Yout=quivers_Yout;
handles.quivers_vxout=quivers_vxout;
handles.quivers_vyout=quivers_vyout;

%% update rtable
rdata=handles.rdata;
APp=prctile(alll,[5,50,95]);
CVp=prctile(vout,[5,50,95]);
handles.errCV=[onedevcv,SEcv,varicv,((CVp(3)-CVp(1))/CVp(2))];

guidata(hObject, handles);

rdata(1,1)=mean(alll);rdata(1,2)=handles.err(1);rdata(1,3)=handles.err(2);rdata(1,4)=handles.err(3);rdata(1,5)=handles.err(4);
rdata(2,1)=mean(vout);rdata(2,2)=handles.errCV(1);rdata(2,3)=handles.errCV(2);rdata(2,4)=handles.errCV(3);rdata(2,5)=handles.errCV(4);
handles.rdata=rdata;
set(handles.rtable,'data',handles.rdata);

%activation time
tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2double(get(handles.actmax,'String'));
actmin=str2double(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);


if actmax < 100
Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);

if actmax < 100
timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
pixar=str2double(get(handles.pixelsize,'String'));
pixar=pixar*pixar;
pixar=pixar/1000000; %convert to mm2
normfac=pixar*allpts;
timdiff=timdiff/normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);

delete(wb)
%maps
Mapchoice_Callback(hObject, eventdata, handles)
% hold off
%scales the values to work with out camera frame rate and resolution
%factor = pix/exposure*100; % converts to cm/sec
%CV = mean(v)*factor;
% disp(['mean_cv: ', num2str(CV), 'cm/sec']);
%
% text(0,1,['CV: ',num2str(CV), 'cm/sec'], 'Units', 'Normalized')


guidata(hObject, handles);
% --- Executes on button press in pushmapapply.
function pushmapapply_Callback(hObject, eventdata, handles)
% hObject    handle to pushmapapply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mapchoice_Callback(hObject, eventdata, handles)



% --- Executes on button press in producemaps.
function producemaps_Callback(hObject, eventdata, ~)
% hObject    handle to producemaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
%handles.filming=0;
%% Store each peak into array
before=str2double(get(handles.beforeGUI,'String'));
after=str2double(get(handles.afterGUI,'String'));

exposure=1/str2double(get(handles.framerate,'String'));
before = round(before/exposure); %1000 because we are dealing with ms
after = round(after/exposure);
if handles.filming == 0
wb=waitbar(0.1,'Preparing Images');
section_choice=get(handles.listbox2,'Value');
end

handles.filming
if handles.filming == 1
handles.filmcount=handles.filmcount+1;
section_choice=handles.filmcount;
section_choice
end
if strcmp(handles.section{section_choice},'Zoomed Section') == 1 && handles.filming ~= 1
sig=handles.averages(handles.q2locs(section_choice,1):handles.q2locs(section_choice,2));
%Peak Parameters
maxfluo = max(handles.averages);
mo = mode(handles.averages);
handles.peakheight = maxfluo*str2double(get(handles.peakhigh,'String'));
minpeakdist = str2double(get(handles.minpeak,'String'));
minpeakdist = ceil(minpeakdist/(1/str2double(get(handles.framerate,'String'))));
% find peaks
[~,m] = findpeaks(sig, 'MINPEAKHEIGHT', handles.peakheight, 'MINPEAKDISTANCE', minpeakdist)
m;
m=m+handles.q2locs(section_choice,1)-1;
else
m=handles.q2locs(section_choice,:);
end
f=m(m~=0);
peaks = (length(f)); % ignores last peak as the signal may cut out

%% OVERLAYING ALL BEATS TO MAKE AN AVERAGE BEAT
% total action potential duration
APtime = before+after;
[~,~,num]=size(handles.images(:,:,:))
% create empty matrix to fill later on
if APtime <= num
overlay = zeros(size(handles.im,1), size(handles.im,2), APtime);
else
overlay = zeros(size(handles.im,1), size(handles.im,2), num);
end
% skip the first and last AP to forgo any possible errors exceeding matrix
% dimensions
if f(1) <= before
startloc =2;
else startloc =1;
end
if f(end)+after > num
endloc=numel(f)-1;
else endloc = numel(f);
end
locRange = startloc:endloc;
if length(handles.q2locs) >1


end
% fill matrix
if isempty(locRange)== 1
errordlg('Peak too close to start/end of file to be analysed with current window settings, next peak analysed')
end

if f(locRange(end))+after > length(handles.images(1,1,:))
locRange=locRange(1:end-1);
end
wsmat=[];
tic

f1=f(locRange);
f1start=(f1(1)-before);
f1end=(f1(end)-after);
handles.imagerange=handles.images(:,:,f1start:f1end);
if length(locRange) > 1 %only need to overlay if more than 1 peak
wsmat=zeros(size(handles.images,1),size(handles.images,2),numel(-before:after),numel(locRange));
for x = -before:after
if f(locRange)+after < length(handles.images(1,1,:))
overlay(:,:,x+before+1) = sum(handles.images(:,:,f(locRange)+x),3)./numel(f);
overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(handles.mask);
wsmat(:,:,x+before+1,locRange)=handles.images(:,:,f(locRange)+x);
end

end
end
handles.wsmat=wsmat;
toc

if length(locRange) == 1 %1 beat
for x = -before:after
overlay(:,:,x+before+1) = (handles.images(:,:,f(locRange)+x));
overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(handles.mask);
end
end
handles.cvimages=overlay;
inversion=get(handles.invertopt,'Value');


%% WRITE TO TIFF STACK
% % normalise
%cos = 26/11/16 with new BL removal overlay all negative, so min and max
if handles.numofpeaksoverall > 1
minI = min(overlay(:));
maxI = max(overlay(:));
%
averageBeat = overlay - minI;
averageBeat = (2^16-1)*averageBeat./(maxI);
%
%make 16 bit
handles.averageBeat = uint16(averageBeat);
end
if handles.numofpeaksoverall == 1
disp('hi')
handles.averageBeat=handles.images;
%handles.cvimages=handles.images;
end
%% Get numbers
handles.outlier=get(handles.apdout,'Value');
if handles.filming == 0
waitbar(0.4,wb,'Producing APD map');
t=str2double(get(handles.t,'String'));
                               
[apmap,meanapd,alll,onedev,var,SE]=mapsbaby(get(handles.aptime1,'Value'),str2double(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2double(get(handles.cmin,'String')),str2double(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2double(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2double(get(handles.apdblnum,'String')),handles.medifilt);

waitbar(0.6,wb,'Producing Isochronal map');

if get(handles.actfittimes,'Value') == 1
[actmap,~,~,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,~,~,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,varCV,SECV]=...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
0,400,str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
[actmap,~,~,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,~,~,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,varCV,SECV]=...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
str2double(get(handles.MINt,'String')),str2double(get(handles.MAXt,'String')),str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
APp=prctile(alll,[5,50,95]);
CVp=prctile(vout,[5,50,95]);
handles.err=[onedev,SE,var,((APp(3)-APp(1))/APp(2))];
handles.errCV=[onedevCV,SECV,varCV,((CVp(3)-CVp(1))/CVp(2))];


rdata=zeros(4,4);
rdata(1,1)=mean(alll);rdata(1,2)=handles.err(1);rdata(1,3)=handles.err(2);rdata(1,4)=handles.err(3);rdata(1,5)=handles.err(4);
rdata(2,1)=mean(vout);rdata(2,2)=handles.errCV(1);rdata(2,3)=handles.errCV(2);rdata(2,4)=handles.errCV(3);rdata(2,5)=handles.errCV(4);

[fmap]=fluo_map(str2double(get(handles.framerate,'String')),handles.I,handles.images,get(handles.tfilt,'Value'),handles.averageBeat);
fmap(isnan(fmap))=0;
salll=fmap(fmap>0);
sp=prctile(salll,[5,50,95]);
rdata(3,1)=mean(salll);
rdata(3,2)=std(salll);
rdata(3,3)=std(salll)/sqrt(numel(salll));
rdata(3,4)=std(salll)*std(salll);
rdata(3,5)=((sp(3)-sp(1))/sp(2));
rdata(rdata==0)=NaN;
handles.rdata=rdata;
set(handles.rtable,'Data',rdata);
disp('2218')
%save AP,CV maps and vectors
handles.apdmap=apmap;
handles.apalll=alll;

handles.actmap=actmap;
handles.act_t=act_t;
handles.vout=vout;
handles.quivers_X=quivers_X;
handles.quivers_Y=quivers_Y;
handles.quivers_vx=quivers_vx;
handles.quivers_vy=quivers_vy;
handles.quivers_Xout=quivers_Xout;
handles.quivers_Yout=quivers_Yout;
handles.quivers_vxout=quivers_vxout;
handles.quivers_vyout=quivers_vyout;

%activation time
tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2double(get(handles.actmax,'String'));
actmin=str2double(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);


if actmax < 100
Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);

if actmax < 100
timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
pixar=str2double(get(handles.pixelsize,'String'));
pixar=pixar*pixar;
pixar=pixar/1000000; %convert to mm2
normfac=pixar*allpts;
timdiff=timdiff/normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);


if strcmp(handles.section{1},'N/A') == 1
set(handles.CLdisp,'String',('N/A - only one peak'));
end
if strcmp(handles.section{1},'N/A') == 0
if strcmp(handles.section{section_choice},'Zoomed Section') == 1
set(handles.CLdisp,'String','N/A - Custom Section');
else
set(handles.CLdisp,'String',[num2str((handles.avgCL(2,section_choice))),' ms (Frequency = ',num2str(1000/round(handles.avgCL(2,section_choice),-1)),' Hz)']);
end
end
delete(wb)
end
guidata(hObject, handles);
axes(handles.mapaxes);
%% Upadte roi selector


% Make MAPS!!!!!
Mapchoice_Callback(hObject, eventdata, handles)
drawnow()
guidata(hObject, handles);
% --- Executes on selection change in threshopt.
function threshopt_Callback(hObject, ~, handles)
% hObject    handle to threshopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.threshopt,'Value') == 1
set(handles.manthresh,'Enable','off')
end
if get(handles.threshopt,'Value') == 2
set(handles.manthresh,'Enable','on')
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns threshopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from threshopt


% --- Executes during object creation, after setting all properties.
function threshopt_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function manthresh_Callback(hObject, eventdata, ~)
handles = guidata(hObject);
pushload_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function manthresh_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cropimage.
function cropimage_Callback(~, ~, ~)


% --- Executes on selection change in imagedisp.
function imagedisp_Callback(hObject, ~, ~)

handles = guidata(hObject);
imchoice=get(handles.imagedisp,'Value');
axes(handles.imageaxes);
cla;
[boundaries] = bwboundaries(handles.mask);

if imchoice == 1
imshow(handles.frame1,[],'InitialMagnification', 400)
colormap('gray')
freezeColors
hold on
for i=1:size(boundaries,1)
plot(boundaries{i}(:,2),boundaries{i}(:,1),'r','LineWidth',2);
end
hold off
end

if imchoice == 2
imshow(handles.fluoim,[],'InitialMagnification', 400)
colormap('jet')
freezeColors
hold on
for i=1:size(boundaries,1)
plot(boundaries{i}(:,2),boundaries{i}(:,1),'k','LineWidth',2);
end
hold off
end


axes(handles.mapaxes);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns imagedisp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagedisp


% --- Executes during object creation, after setting all properties.
function imagedisp_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cropbox.
function cropbox_Callback(~, ~, ~)



% --- Executes on button press in exportvalues.
function exportvalues_Callback(hObject, ~, ~)

handles = guidata(hObject);
[filename,pathname] = uiputfile({'*.csv';'*.txt';'*.mat'}, 'Save Map values (Isochronal maps with vector overlay can only be saved as .mat files)');
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
choice=get(handles.Mapchoice,'Value');

map=handles.holdmap;
% mat file save
if strcmp('.mat',ext) == 1
if choice == 1
APD_Dist=map;
save(file,'APD_Dist');
end
if choice == 2;
activation_time=map;
save(file,'activation_time');
end
if choice == 3;
activation_time=map;
xpositions=X_pos;
ypositions=Y_pos;
xvelocities=X_vel;
yvelocities=Y_vel;
velocities=total_vel;
save(file,'activation_time','xpositions','ypositions','xvelocities','yvelocities','velocities','fractional_up');
end
if choice == 4;
activation_time=map;
xpositions=Xout_pos;
ypositions=Yout_pos;
xvelocities=Xout_vel;
yvelocities=Yout_vel;
velocities=total_velout;
save(file,'activation_time','xpositions','ypositions','xvelocities','yvelocities','velocities','fractional_up');
end
if choice == 6
DI=handles.holdmap;
save(file,'DI')
end
if choice == 8
tau=handles.holdmap;
save(file,'tau')
end
if choice == 9
[~,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, onedevCV,varCV,SECV]=...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
0,200,str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
cvmean = mean(vout);
CVdev=onedevCV;

save(file,'RT50','t1080','signal_level','gofmap1080','cvmean','CVdev')
disp('done')
end
end

% csv and txt files
if strcmp('.csv',ext) == 1 || strcmp('.txt',ext) == 1
cho = questdlg('How would you like to export the values?', ...
'Export', ...
'Map','List','List');
switch cho
case 'Map'
T=table(map);
writetable(T,file,'Delimiter',',','WriteVariableNames',false);
if choice == 2 || choice == 3 || choice == 4
cvfile=[pathname,'CV_',filename];
T2=table(handles.holdcvmap);
writetable(T2,cvfile,'Delimiter',',','WriteVariableNames',false);
end
case 'List'
listmap=reshape(map,numel(map),1);
listmap=listmap(listmap>0);
listmap=listmap(isnan(listmap)==0);
T=table(listmap);
writetable(T,file,'Delimiter',',','WriteVariableNames',false);
if choice == 2 || choice == 3 || choice == 4
cvfile=[pathname,'CV_',filename];
CVmap=handles.holdcvmap;
listcvmap=reshape(CVmap,numel(CVmap),1);
listcvmap=listcvmap(listcvmap>0);
listcvmap=listcvmap(isnan(listcvmap)==0);
T=table(listcvmap);
writetable(T,cvfile,'Delimiter',',','WriteVariableNames',false);
end
end

end


% --- Executes on button press in act_movie.
function act_movie_Callback(hObject, ~, ~)

handles = guidata(hObject);
meme=0;
gifsave=questdlg('Save File?','Gif Save','Yes','No','Yes');
switch gifsave
case 'Yes'
meme=1;
case 'No'
meme=0;
end
map=handles.actmap;
if meme == 1
[a,b]=uiputfile('*.gif');
filename=[b,a];
end
h=figure;
hold on
imshow(map,[],'InitialMagnification', 800);
map=map;
if get(handles.isoopt,'Value') == 1
mini=0;
maxi=max(max(map)); %ms
end
if get(handles.isoopt,'Value') == 2
maxi=str2double(get(handles.isomax,'String'));
mini=str2double(get(handles.isomin,'String'));
end
maxiall=max(max(map));
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [0, 0, 0];
colormap(jetcolormap);

caxis([0 maxi]);
delay=0.01;
for j = 1:0.1:ceil(maxiall);
mapmask=(map<j);
A=map.*mapmask;
imshow(A, [mini ceil(maxi)], 'Colormap',jetcolormap, 'InitialMagnification', 400)
if meme == 1
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if j == 1
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay);
else
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
end
end
end





function isomin_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function isomin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function isomax_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function isomax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in isoopt.
function isoopt_Callback(hObject, eventdata, handles)

handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function isoopt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in applyiso.
function applyiso_Callback(hObject, eventdata, handles)

handles = guidata(hObject);
listbox2_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)

handles = guidata(hObject);
axes(handles.mapaxes);
choice = questdlg('Video or Image?', ...
'Thing', ...
'Video','Image','Video');
switch choice
case 'Video'
[filename,pathname] = uiputfile({'*.avi'}, 'Save .avi image video file of currently displayed maps across all sections' );
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
%         if isdepolyed == 0
%         cd(pathname)
%         end
handles.filming = 1;
numsec=length(handles.section);
vidobj = VideoWriter(file);
vidobj.FrameRate=1;
open(vidobj);
set(gca,'nextplot','replacechildren');

guidata(hObject, handles);
%wb=waitbar(0.1,'Producing video file','WindowStyle', 'modal');
for k=1:numsec
handles.filmcount = k;
guidata(hObject,handles);
producemaps_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
axes(handles.mapaxes);
currFrame = getframe;
writeVideo(vidobj,currFrame);
0.1+0.9*(k/numsec)
%waitbar((0.1+0.9*(k/numsec)),wb);
end
close(vidobj);
%delete(wb)
set(handles.listbox2,'Value',numsec)
handles.filming = 0;
case 'Image'
handles.filming = 1;
numsec=length(handles.section);
for i=1:numsec
handles.filmcount = i;
handles.b2bimage=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles);
GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;
ax=axes;clf;
new_handle=copyobj(handles.mapaxes,fig);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[i/numsec-0.2 i/numsec-0.2 i/numsec i/numsec])
%set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
guidata(hObject, handles);
end
handles.filming=0
guidata(hObject, handles);
end


% --- Executes on button press in segEP.
function segEP_Callback(~, ~, ~)
segEP

% --- Executes on button press in fold.
function fold_Callback(~, ~, ~)



function framerate_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function pixelsize_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function pixelsize_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function binnumber_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function binnumber_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getpixelinfo.
function getpixelinfo_Callback(hObject, ~, handles)
% hObject    handle to getpixelinfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
compare
guidata(hObject, handles);

% --- Executes on button press in phasemap.
function phasemap_Callback(~, ~, ~)
% hObject    handle to phasemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

phasemapping


% --- Executes on button press in compare.
function compare_Callback(~, ~, ~)
pixelinfo

% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(~, ~, ~)
conduction


function beforeGUI_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function beforeGUI_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function afterGUI_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function afterGUI_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in velalgo.
function velalgo_Callback(hObject, eventdata, handles)

wb=waitbar(0.5,'Producing Isochronal map');
guidata(hObject, handles);
if get(handles.actfittimes,'Value') == 1
[actmap,~,~,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,~,~,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,varCV,SECV]=...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
0,400,str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
[actmap,~,~,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,~,~,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,varCV,SECV]=...
cvmap(str2double(get(handles.pixelsize,'String')),str2double(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2double(get(handles.minvel,'String')),str2double(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
str2double(get(handles.MINt,'String')),str2double(get(handles.MAXt,'String')),str2double(get(handles.winsize,'String')),str2double(get(handles.beforeGUI,'String')),str2double(get(handles.wint,'String')),0,str2double(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2double(get(handles.splineN,'String')));
end

CVp=prctile(vout,[5,50,95]);

handles.errCV=[onedevCV,SECV,varCV,((CVp(3)-CVp(1))/CVp(2))];

rdata=handles.rdata;

rdata(2,1)=mean(vout);rdata(2,2)=handles.errCV(1);rdata(2,3)=handles.errCV(2);rdata(2,4)=handles.errCV(3);rdata(2,5)=handles.errCV(4);

handles.rdata=rdata;
set(handles.rtable,'Data',rdata);
disp('2768')
%save CV maps and vectors

handles.actmap=actmap;
handles.act_t=act_t;
handles.quivers_X=quivers_X;
handles.quivers_Y=quivers_Y;
handles.quivers_vx=quivers_vx;
handles.quivers_vy=quivers_vy;
handles.quivers_Xout=quivers_Xout;
handles.quivers_Yout=quivers_Yout;
handles.quivers_vxout=quivers_vxout;
handles.quivers_vyout=quivers_vyout;
delete(wb)
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns velalgo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velalgo


% --- Executes during object creation, after setting all properties.
function velalgo_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in actopt.
function actopt_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function actopt_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end




function actmin_Callback(hObject, ~, ~)
handles = guidata(hObject);

tim=handles.act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2double(get(handles.actmax,'String'));
actmin=str2double(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
Imax=Imax(1);
else
    Imax=max(tim);
end
Imin=Imin(1);


if actmax < 100
timmax=Imax*0.01;
else
    timmax=Imax;
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
normfac=225/allpts;  %why 225?
timdiff=timdiff*normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of actmin as text
%        str2double(get(hObject,'String')) returns contents of actmin as a double


% --- Executes during object creation, after setting all properties.
function actmin_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function actmax_Callback(hObject, ~, ~)

handles = guidata(hObject);

tim=handles.act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2double(get(handles.actmax,'String'));
actmin=str2double(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
Imax=Imax(1);
else Imax=max(tim);
end

Imin=Imin(1);


if actmax < 100
timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
timdiff=timmax;
end

if get(handles.checkbox8,'Value')==1
normfac=225/allpts;
timdiff=timdiff*normfac;
end
set(handles.actquote, 'String', [num2str(timdiff),' ms']);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function actmax_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, ~, ~)

handles = guidata(hObject);


tim=handles.act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2double(get(handles.actmax,'String'));
actmin=str2double(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
Imax=Imax(1);
else
    Imax=max(tim);
end
Imin=Imin(1);


if actmax < 100
timmax=Imax*0.01;
else
    timmax=Imax;
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
pixar=str2double(get(handles.pixelsize,'String'));
pixar=pixar*pixar;
pixar=pixar/1000000; %convert to mm2
normfac=pixar*allpts;
timdiff=timdiff/normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);




function t_Callback(hObject, eventdata, ~)
handles = guidata(hObject);
t=str2double(get(handles.t,'String'));
[apmap,meanapd,alll,onedev,var,SE]=mapsbaby(get(handles.aptime1,'Value'),str2double(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2double(get(handles.cmin,'String')),str2double(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2double(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2double(get(handles.apdblnum,'String')),handles.medifilt);
%alll=apmap(apmap>0);
APp=prctile(alll,[5,50,95]);
handles.err=[onedev,SE,var,((APp(3)-APp(1))/APp(2))];
%save AP,CV maps and vectors
handles.apdmap=apmap;
handles.apalll=alll;
rdata=handles.rdata;
rdata(1,1)=mean(alll);rdata(1,2)=handles.err(1);rdata(1,3)=handles.err(2);rdata(1,4)=handles.err(3);rdata(1,5)=handles.err(4);
handles.rdata=rdata;
guidata(hObject,handles)
set(handles.rtable,'Data',rdata);
disp('2977')
mapcho=get(handles.Mapchoice,'Value');
if mapcho == 1
Mapchoice_Callback(hObject, eventdata, handles);
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function t_CreateFcn(hObject, ~, ~)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdscale.
function apdscale_Callback(hObject, eventdata, ~)
% hObject    handle to apdscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function apdscale_CreateFcn(hObject, ~, ~)
% hObject    handle to apdscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function MINt_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function MINt_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function MAXt_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function MAXt_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in actfittimes.
function actfittimes_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function actfittimes_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function winsize_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function winsize_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function taustart_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
if choice == 8
Mapchoice_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function taustart_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function r2cut_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
if choice == 8
Mapchoice_Callback(hObject, eventdata, handles)
end



% --- Executes during object creation, after setting all properties.
function r2cut_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function apdblnum_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function apdblnum_CreateFcn(hObject, ~, ~)


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdbl.
function apdbl_Callback(hObject, eventdata, handles)
t_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function apdbl_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function taufinish_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function taufinish_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function peakhigh_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function peakhigh_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function scal_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function scal_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function wint_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function wint_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in altanal.
function altanal_Callback(hObject, ~, ~)

alternangui


% --- Executes on button press in removef.
function removef_Callback(~, ~, ~)


% --- Executes on selection change in aptime1.
function aptime1_Callback(hObject, eventdata, handles)

t_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function aptime1_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in configure.
function configure_Callback(~, ~, handles)

% Construct a questdlg
choice = questdlg('Would you like to load a configuration file or save current settings?', ...
'Config File', ...
'New File', 'Load Settings','Load Settings');
% Handle response
switch choice
case 'New File'
[fname,PathName] = uiputfile('*.txt')
filename=[PathName,fname];

% create file for writing too
fileID = fopen(filename,'w');
dt=datetime('now'); dt=datestr(dt);
fprintf(fileID,'Configuration File for use in ElectroMap\r\n');
fprintf(fileID,['Created: ',dt,'\r\n'])
fprintf(fileID,['Notes:\r\n\r\n\r\n'])
fprintf(fileID,['--- DO NOT EDIT VARIABLE NAMES OR REMOVE ! BELOW THIS POINT, {} = units, () = settings in ElectroMap or Notes ---\r\n\r\n'])

%cell arrays for popdown menu strings
threshopt_string=get(handles.threshopt,'String');
thershopt_string=threshopt_string{get(handles.threshopt,'Value')};
sfilt_string=get(handles.sfilt,'String');
sfilt_string=sfilt_string{get(handles.sfilt,'Value')};
segchoice_string=get(handles.segchoice,'String');
segchoice_string=segchoice_string{get(handles.segchoice,'Value')}
BLopt_string=get(handles.BLopt,'String');
BLopt_string=BLopt_string{get(handles.BLopt,'Value')}
tfilt_string=get(handles.tfilt,'String');
tfilt_string=tfilt_string{get(handles.tfilt,'Value')};
apdbl_string=get(handles.apdbl,'String');
apdbl_string=apdbl_string{get(handles.apdbl,'Value')};
aptime1_string=get(handles.aptime1,'String');
aptime1_string=aptime1_string{get(handles.aptime1,'Value')};
velalgo_string=get(handles.velalgo,'String');
velalgo_string=velalgo_string{get(handles.velalgo,'Value')};
actfittimes_string=get(handles.actfittimes,'String');
actfittimes_string=actfittimes_string{get(handles.actfittimes,'Value')};
velout_string=get(handles.velout,'String');
velout_string=velout_string{get(handles.velout,'Value')};
apdout_string=get(handles.apdout,'String');
apdout_string=apdout_string{get(handles.apdout,'Value')};
apdscale_string=get(handles.apdscale,'String');
apdscale_string=apdscale_string{get(handles.apdscale,'Value')};
% get settings and put into file
fprintf(fileID,['framerate=',get(handles.framerate,'String'),'! {kHz} !\r\n']);
fprintf(fileID,['pixelsize=',get(handles.pixelsize,'String'),'! {ms} !!\r\n']);
fprintf(fileID,['threshopt=',num2str(get(handles.threshopt,'Value')),'! (',thershopt_string,')',' !\r\n']);
fprintf(fileID,['manthresh=',get(handles.manthresh,'String'),'! {Percent} (Change from automatically generated threshold) !\r\n']);
fprintf(fileID,['sfilt=',num2str(get(handles.sfilt,'Value')),'! (',sfilt_string,') !\r\n']);
fprintf(fileID,['sfiltsize=',get(handles.sfiltsize,'String'),'! {Pixels} !\r\n']);
fprintf(fileID,['segchoice=',num2str(get(handles.segchoice,'Value')),'! (',segchoice_string,') !\r\n']);
fprintf(fileID,['segsize=',get(handles.segsize,'String'),'! !\r\n']);
fprintf(fileID,['invertopt=',num2str(get(handles.invertopt,'Value')),'! (1 means invert signal) !\r\n']);
fprintf(fileID,['BLopt=',num2str(get(handles.BLopt,'Value')),'! (',BLopt_string,') !\r\n']);
fprintf(fileID,['thlen=',get(handles.thlen,'String'),'! {ms}, (Length of top-hat filter) !\r\n']);
fprintf(fileID,['tfilt=',num2str(get(handles.tfilt,'Value')),'! (',tfilt_string,') !\r\n']);
fprintf(fileID,['minpeak=',get(handles.minpeak,'String'),'! {ms} !\r\n']);
fprintf(fileID,['peakhigh=',get(handles.peakhigh,'String'),'! !\r\n']);
fprintf(fileID,['minnum=',get(handles.minnum,'String'),'! !\r\n']);
fprintf(fileID,['minbound=',get(handles.minbound,'String'),'! {ms} !\r\n']);
fprintf(fileID,['beforeGUI=',get(handles.beforeGUI,'String'),'! {ms} !\r\n']);
fprintf(fileID,['afterGUI=',get(handles.afterGUI,'String'),'! {ms} !\r\n']);
fprintf(fileID,['apdbl=',num2str(get(handles.apdbl,'Value')),'! (',apdbl_string,') !\r\n']);
fprintf(fileID,['apdblnum=',get(handles.apdblnum,'String'),'! {ms} !\r\n']);
fprintf(fileID,['aptime1=',num2str(get(handles.aptime1,'Value')),'! (',aptime1_string,') !\r\n']);
fprintf(fileID,['taustart=',get(handles.taustart,'String'),'! {Percent} !\r\n']);
fprintf(fileID,['taufinish=',get(handles.taufinish,'String'),'! {Percent} !\r\n']);
fprintf(fileID,['r2cut=',get(handles.r2cut,'String'),'! !\r\n']);
fprintf(fileID,['velalgo=',num2str(get(handles.velalgo,'Value')),'! (',velalgo_string,') (Activation Measure) !\r\n']);
fprintf(fileID,['isoopt=',num2str(get(handles.isoopt,'Value')),'! !\r\n']);
fprintf(fileID,['isomin=',get(handles.isomin,'String'),'! {ms} !\r\n']);
fprintf(fileID,['isomax=',get(handles.isomax,'String'),'! {ms} !\r\n']);
fprintf(fileID,['actfittimes=',num2str(get(handles.actfittimes,'Value')),'! (',actfittimes_string,') !\r\n']);
fprintf(fileID,['MINt=',get(handles.MINt,'String'),'! {ms} (Minimum activation time for multi-vector fit) !\r\n']);
fprintf(fileID,['MAXt=',get(handles.MAXt,'String'),'! {ms} (Maximum activation time for multi-vector fit) !\r\n']);
fprintf(fileID,['velout=',num2str(get(handles.velout,'Value')),'! (',velout_string,') (Local Velocity outlier removal) !\r\n']);
fprintf(fileID,['minvel=',get(handles.minvel,'String'),'! {ms} (Minimum calcualted velocity that is not discarded) !\r\n']);
fprintf(fileID,['maxvel=',get(handles.maxvel,'String'),'! {ms} (Maximum calculated velocity that is not discarded) !\r\n']);
fprintf(fileID,['winsize=',get(handles.winsize,'String'),'!{Pixels} (Local window size) !\r\n']);
fprintf(fileID,['scal=',get(handles.scal,'String'),'! (Size of overlaid velocity vectors) !\r\n']);
fprintf(fileID,['wint=',get(handles.wint,'String'),'! (Maximum time diffrence allowed in local window fit) !\r\n']);
fprintf(fileID,['apdout=',num2str(get(handles.apdout,'Value')),'! (',apdout_string,') (APD/CaD outlier removal) !\r\n']);
fprintf(fileID,['apdscale=',num2str(get(handles.apdscale,'Value')),'! (',apdscale_string,') !\r\n']);
fprintf(fileID,['cmin=',get(handles.cmin,'String'),'! {ms} (manual colour map minimum)!\r\n']);
fprintf(fileID,['cmax=',get(handles.cmax,'String'),'! {ms} (manual colour map maximum)!\r\n']);
fprintf(fileID,['t=',get(handles.t,'String'),'! {Percent} (APD/CaD)!\r\n']);
fprintf(fileID,['checkbox8=',num2str(get(handles.checkbox8,'Value')),'! (1 means normalised to {ms/mm2}, 0 absoulte in {ms})!\r\n']);
fprintf(fileID,['actmin=',get(handles.actmin,'String'),'! {Percent}!\r\n']);
fprintf(fileID,['actmax=',get(handles.actmax,'String'),'! {Percent}!\r\n']);
fprintf(fileID,['binnumber=',get(handles.binnumber,'String'),'!!\r\n']);

% close file
fclose(fileID);

case 'Load Settings'
[fname,PathName] = uigetfile('*.txt');
filename=[PathName,fname];
%% open file for reading
fid=fopen(filename,'r','b')
fstr=fread(fid,'int8=>char')';
fclose(fid);

%% Update GUI
set(handles.framerate,'String',varEM(fstr,'framerate',0))
set(handles.pixelsize,'String',varEM(fstr,'pixelsize',0))
set(handles.threshopt,'Value',varEM(fstr,'threshopt',1))
set(handles.manthresh,'String',varEM(fstr,'manthresh',0))
set(handles.sfilt,'Value',varEM(fstr,'sfilt',1))
set(handles.sfiltsize,'String',varEM(fstr,'sfiltsize',0))
set(handles.segchoice,'Value',varEM(fstr,'segchoice',1))
set(handles.segsize,'String',varEM(fstr,'segsize',0))
set(handles.invertopt,'Value',varEM(fstr,'invertopt',1))
set(handles.BLopt,'Value',varEM(fstr,'BLopt',1))
set(handles.tfilt,'Value',varEM(fstr,'tfilt',1))
set(handles.minpeak,'String',varEM(fstr,'minpeak',0))
set(handles.peakhigh,'String',varEM(fstr,'peakhigh',0))
set(handles.minnum,'String',varEM(fstr,'minnum',0))
set(handles.minbound,'String',varEM(fstr,'minbound',0))
set(handles.beforeGUI,'String',varEM(fstr,'beforeGUI',0))
set(handles.afterGUI,'String',varEM(fstr,'afterGUI',0))
set(handles.apdbl,'Value',varEM(fstr,'apdbl',1))
set(handles.apdblnum,'String',varEM(fstr,'apdblnum',0))
set(handles.aptime1,'Value',varEM(fstr,'aptime1',1))
set(handles.taustart,'String',varEM(fstr,'taustart',0))
set(handles.taufinish,'String',varEM(fstr,'taufinish',0))
set(handles.r2cut,'String',varEM(fstr,'r2cut',0))
set(handles.velalgo,'Value',varEM(fstr,'velalgo',1))
set(handles.isoopt,'Value',varEM(fstr,'isoopt',1))
set(handles.isomin,'String',varEM(fstr,'isomin',0))
set(handles.isomax,'String',varEM(fstr,'isomax',0))
set(handles.actfittimes,'Value',varEM(fstr,'actfittimes',1))
set(handles.MINt,'String',varEM(fstr,'MINt',0))
set(handles.MAXt,'String',varEM(fstr,'MAXt',0))
set(handles.velout,'Value',varEM(fstr,'velout',1))
set(handles.minvel,'String',varEM(fstr,'minvel',0))
set(handles.maxvel,'String',varEM(fstr,'maxvel',0))
set(handles.winsize,'String',varEM(fstr,'winsize',0))
set(handles.scal,'String',varEM(fstr,'scal',0))
set(handles.wint,'String',varEM(fstr,'wint',0))
set(handles.apdout,'Value',varEM(fstr,'apdout',1))
set(handles.apdscale,'Value',varEM(fstr,'apdscale',1))
set(handles.cmin,'String',varEM(fstr,'cmin',0))
set(handles.cmax,'String',varEM(fstr,'cmax',0))
set(handles.t,'String',varEM(fstr,'t',0))
set(handles.checkbox8,'Value',varEM(fstr,'checkbox8',1))
set(handles.actmin,'String',varEM(fstr,'actmin',0))
set(handles.actmax,'String',varEM(fstr,'actmax',0))
set(handles.binnumber,'String',varEM(fstr,'binnumber',0))
set(handles.thlen,'String',varEM(fstr,'thlen',0))
end



function thlen_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function thlen_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colmap.
function colmap_Callback(hObject, eventdata, handles)

handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function colmap_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in drawcon.
function drawcon_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)





% --- Executes on button press in roibutton.
function roibutton_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
choice = questdlg('Save Current ROI or load previous?', ...
'ROI', ...
'Save ROI','Load ROI','Save ROI');
% Handle response
switch choice
case 'Save ROI'
[filename,pathname] = uiputfile({'*.txt'}, 'Save ROI in text file');
file=[pathname,filename];
savemask=handles.mask;
figure,
imshow(savemask,[])
dlmwrite(file,savemask);
lmask=maskload(file)
figure,
imshow(lmask,[])
size(savemask)
size(lmask)
case 'Load ROI'
[filename,pathname] = uigetfile('*.txt','Select the ROI File');
file=[pathname,filename];
handles.loadedmask=maskload(file)
handles.herefromroiload=1;
guidata(hObject, handles);
pushload_Callback(hObject, eventdata, handles)
handles = guidata(hObject); %%update handles after doing OMimload
handles.herefromroiload=0;
guidata(hObject, handles);
end
guidata(hObject, handles);


% --- Executes on button press in rawvid.
function rawvid_Callback(hObject, ~, ~)

handles = guidata(hObject);
axes(handles.imageaxes)
[rows, cols]=size(handles.im);
images=handles.images;
mask=handles.mask;
images=imcomplement(images);
im=handles.im;
im=double(im);
im=im-min(min(im));
im=im./max(max(im));
im=im*65535;
im=uint16(im);

savegif=1;
if savegif == 1
[a,b]=uiputfile('*.gif');
filename=[b,a];
end
prompt = {'Fluorescence threshold (0-1):','Video Start (s):','Video End (s)','Normalise? (0=no, 1=yes)'};
dlg_title = 'Raw Video Options';
num_lines = 1;
exposure=1/str2double(get(handles.framerate,'String'));
defaultans = {'0.2','0',num2str(size(images,3)/1000*exposure),'1'};
opts = inputdlg(prompt,dlg_title,num_lines,defaultans);
flthresh=str2double(opts{1});
istart=str2double(opts{2});

iend=str2double(opts{3});
normF=str2double(opts{4});
%change is to frame #
istart=round((istart/exposure)*1000);
iend=round((iend/exposure)*1000);
if istart == 0
istart = 1;
end

if iend > size(images,3)
iend = size(images,3);
end

for r=1:rows
for c=1:cols
sig=squeeze(images(r,c,:));
sig=sig-min(sig);
sig=double(sig);
if normF == 1
sig=(sig./max(sig))*65535;
end
sig=uint16(sig);
images(r,c,:)=sig;
end
end
wb=waitbar(0,'Saving Raw Video');
images=double(images);
mask=double(mask);
maxval=max(max(max(images)));
background = repmat(im, [1, 1, 3]);
for i =istart:iend
waitbar(i/(iend-istart),wb,'Saving Raw Video');
combinedImage = background;
foreground=images(:,:,i).*mask;
foreground=foreground./maxval;
foreground(foreground < flthresh) = 0;
foregroundColourised = colouriseData(foreground, 'j',flthresh,1);
c1 = combinedImage(:, :, 1);
c2 = combinedImage(:, :, 2);
c3 = combinedImage(:, :, 3);

f1 = foregroundColourised(:, :, 1);
f2 = foregroundColourised(:, :, 2);
f3 = foregroundColourised(:, :, 3);

c1(sum(foreground, 3) ~= 0) = f1(sum(foreground, 3) ~= 0);
c2(sum(foreground, 3) ~= 0) = f2(sum(foreground, 3) ~= 0);
c3(sum(foreground, 3) ~= 0) = f3(sum(foreground, 3) ~= 0);

combinedImage(:, :, 1) = c1;
combinedImage(:, :, 2) = c2;
combinedImage(:, :, 3) = c3;

hold on

axis image;
axis off;

if savegif == 1
delay=0.01;
[imind,cm] = rgb2ind(combinedImage,256);
% Write to the GIF File
if i == istart
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay);
else
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
end
end
end
delete(wb)

% --- Executes on button press in resegment.
function resegment_Callback(hObject, eventdata, handles)

handles.herefromsegmentpush=1;
guidata(hObject, handles);
pushprocess_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.herefromsegmentpush = 0;
guidata(hObject, handles);
producemaps_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
guidata(hObject, handles);

% --- Executes on button press in B2B.
function B2B_Callback(hObject, eventdata, handles)

handles.herefromsegmentpush=1;
set(handles.segsize,'String',1);
set(handles.segchoice,'Value',2);
guidata(hObject, handles);
pushprocess_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.herefromsegmentpush = 0;
guidata(hObject, handles);
producemaps_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
guidata(hObject, handles);


% --- Executes on selection change in winopt.
function winopt_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function winopt_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segsignal.
function segsignal_Callback(hObject, eventdata, handles)

handles.herefromsegmentpush=1;
guidata(hObject, handles);
pushprocess_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 



% --- Executes during object creation, after setting all properties.
function segsignal_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function sfiltsigma_Callback(~, ~, ~)


% --- Executes during object creation, after setting all properties.
function sfiltsigma_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in usespline.
function usespline_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function usespline_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function splineN_Callback(~, ~, ~)



% --- Executes during object creation, after setting all properties.
function splineN_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(~, ~, ~)


% --------------------------------------------------------------------
function freqmapopt_Callback(hObject, eventdata, ~)

handles = guidata(hObject); 
prompt = {'Minimum Frequency (Hz):','Maximum Frequency (Hz):','Frequency Bin Size (Hz)','Window? 0 = no, 1 = hann'};
dims = [1 35];
definput = {num2str(handles.fmin),num2str(handles.fmax),num2str(handles.fbin),num2str(handles.dfwin)};
answer = inputdlg(prompt,'Frequnecy Mapping Options',dims,definput)
handles.fmin=str2double(answer{1});
handles.fmax=str2double(answer{2});
handles.fbin=str2double(answer{3});
handles.dfwin=str2double(answer{4});
guidata(hObject, handles);
if get(handles.Mapchoice,'Value') == 5;
Mapchoice_Callback(hObject, eventdata, handles)
end
guidata(hObject, handles);




% --------------------------------------------------------------------
function ROInum_Callback(hObject, ~, ~)

handles = guidata(hObject); 
prompt = {'Number of ROIs:','Remove overlapping pixels? (0=Yes) (1=No)'};
dims = [1 35];
definput = {num2str(handles.roinum),num2str(handles.roisum)};
answer = inputdlg(prompt,'ROI options',dims,definput)
handles.roinum=str2double(answer{1});
handles.roisum=str2double(answer{2});
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(~, ~, ~)


% --------------------------------------------------------------------
function Untitled_5_Callback(~, ~, ~)


%% ColourMaps
function coljet_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',1);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);



% --------------------------------------------------------------------
function colhsv_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',2);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function colhot_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',3);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colcool_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',4);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colparula_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',5);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function colspring_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',6);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colsummer_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',7);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colautumn_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',8);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colwinter_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
set(handles.colmap,'Value',9);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_3_Callback(~, ~, ~)


% --------------------------------------------------------------------
function Untitled_4_Callback(~, ~, ~)


% --------------------------------------------------------------------
function Untitled_6_Callback(~, ~, ~)


% --------------------------------------------------------------------
function bgblack_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
handles.bgcol='k';
handles.bgon=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function bgwhite_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
handles.bgcol='w';
handles.bgon=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function bgtran_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
handles.bgon=0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function snrcalc_Callback(hObject, eventdata, ~)
handles = guidata(hObject); 
prompt = {'Signal before time (from peak) (ms):','Signal after time (from peak) (ms):'};
dims = [1 35];
definput = {num2str(handles.snrt1),num2str(handles.snrt2)};
answer = inputdlg(prompt,'Frequnecy Mapping Options',dims,definput);
handles.snrt1=str2double(answer{1});
handles.snrt2=str2double(answer{2});
guidata(hObject, handles);
if get(handles.Mapchoice,'Value') == 10 || get(handles.Mapchoice,'Value') == 1
Mapchoice_Callback(hObject, eventdata, handles)
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function ttpset_Callback(hObject, eventdata, ~)

handles = guidata(hObject); 
prompt = {'Start Point (%):','End Point (%):'};
dims = [1 35];
definput = {num2str(handles.ttpstart),num2str(handles.ttpend)};
answer = inputdlg(prompt,'Frequnecy Mapping Options',dims,definput);
handles.ttpstart=str2double(answer{1});
handles.ttpend=str2double(answer{2});
guidata(hObject, handles);
if get(handles.Mapchoice,'Value') == 7
Mapchoice_Callback(hObject, eventdata, handles)
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_7_Callback(~, ~, ~)



% --------------------------------------------------------------------
function connnnnnnnnnn_Callback(~, ~, ~)


% --------------------------------------------------------------------
function conoff_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
handles.drawcon=0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function conon_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
handles.drawcon=1;
prompt = {'Contour spacing (map units):'};
dims = [1 35];
if isempty(handles.conbon) == 1
handles.framerate
conbount=1/str2double(get(handles.framerate,'String'));
handles.conbon=num2str(conbount);
end
definput = {num2str(handles.conbon)};
answer = inputdlg(prompt,'Contour Setting',dims,definput);
handles.conbon=answer{1};
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function nomedifilt_Callback(hObject, eventdata, ~)

handles = guidata(hObject);
handles.medifilt=0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function yesmedifilt_Callback(hObject, eventdata, ~)
handles = guidata(hObject);
handles.medifilt=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)

function squareROI_Callback(~,~,~)


% --------------------------------------------------------------------
function frameremoval_Callback(hObject, eventdata, handles)
% hObject    handle to frameremoval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.drawcon=1;
prompt = {'Frames before pulse to remove:','Frames after pulse to remove:'};
dims = [1 35];
definput = {num2str(handles.pbefore),num2str(handles.pafter)};
answer = inputdlg(prompt,'Contour Setting',dims,definput);
handles.pbefore=str2double(answer{1});
handles.pafter=str2double(answer{2});
guidata(hObject, handles);
