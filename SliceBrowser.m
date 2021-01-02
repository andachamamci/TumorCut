function varargout = SliceBrowser(varargin)
%SLICEBROWSER M-file for SliceBrowser.fig
%       SliceBrowser is an interactive viewer of 3D volumes, 
%       it shows 3 perpendicular slices (XY, YZ, ZX) with 3D pointer.
%   Input:  a) VOLUME - a 3D matrix with volume data
%           b) VOLUME - a 4D matrix with volume data over time
%   Control:
%       - Clicking into the window changes the location of 3D pointer.
%       - 3D pointer can be moved also by keyboard arrows.
%       - Pressing +/- will switch to next/previous volume.
%       - Pressing 1,2,3 will change the focus of current axis.
%       - Pressing 'e' will print the location of 3D pointer.
%   Example of usage:
%       load mri.dat
%       volume = squeeze(D);
%       SliceBrowser(volume);
%
% Author: Marian Uhercik, CMP, CTU in Prague
% Web: http://cmp.felk.cvut.cz/~uhercik/3DSliceViewer/3DSliceViewer.htm
% Last Modified by 05-Aug-2008

% Documentation generated GUIDE:
%
%SLICEBROWSER M-file for SliceBrowser.fig
%      SLICEBROWSER, by itself, creates a new SLICEBROWSER or raises the existing
%      singleton*.
%
%      H = SLICEBROWSER returns the handle to a new SLICEBROWSER or the handle to
%      the existing singleton*.
%
%      SLICEBROWSER('Property','Value',...) creates a new SLICEBROWSER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SliceBrowser_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SLICEBROWSER('CALLBACK') and SLICEBROWSER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SLICEBROWSER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SliceBrowser

% Last Modified by GUIDE v2.5 20-Feb-2015 14:47:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SliceBrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @SliceBrowser_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%disp('gui_mainfcn called...');

% End initialization code - DO NOT EDIT


% --- Executes just before SliceBrowser is made visible.
function SliceBrowser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SliceBrowser,

addpath('./NIFTI_20100413');
handles.output = [];%hObject;

[handles.fname,handles.pname,filterindex] = uigetfile('*.*');

disp([handles.pname handles.fname]);
%nii = load_nii([pname fname]);
handles.nii =load_untouch_nii([handles.pname handles.fname]);
handles.Spacing = [handles.nii.hdr.dime.pixdim(2) handles.nii.hdr.dime.pixdim(3) handles.nii.hdr.dime.pixdim(4)];
handles.MRI{1} = uint16(handles.nii.img);

%path('C:\Users\andachamamci\Desktop\MICCAI 2012-Tumor Seg Challenge\TumorSegmentationCode\ReadData3D\',path);
%MRI = ReadData3D;
%MRI = uint16(MRI);
handles.Spacing = [1 1 1];
%delete nii;
%disp('Default Output assigned...');
% Update handles structure
%guidata(hObject, handles);

% UIWAIT makes SliceBrowser wait for user response (see UIRESUME)
%uiwait(handles.figure1);

if (length(handles.MRI) <=0)
    error('Input volume has not been specified.');
end;
volume = handles.MRI{1};
if (ndims(volume) ~= 3 && ndims(volume) ~= 4)
    error('Input volume must have 3 or 4 dimensions.');
end;
handles.volume = volume;

contour = [];
if (length(handles.MRI) > 1)
    contour = handles.MRI{2};
    if (ndims(contour) ~= ndims(volume))
        error('Input volume and contour must have same dimensions.');
    end
    if (size(contour) ~= size(volume))
        error('Input volume and contour must be in same size.');
    end
end
handles.contour = not(contour == 0);

% set main wnd title
set(gcf, 'Name', 'Slice Viewer')

% init 3D pointer
vol_sz = size(volume); 
if (ndims(volume) == 3)
    vol_sz(4) = 1;
end;
pointer3dt = floor(vol_sz/2)+1;
handles.pointer3dt = pointer3dt;
handles.vol_sz = vol_sz;

handles.IsSelected = false;
handles.HasLine = false;

handles.FirstP = [];
handles.SecondP = [];
plot3slices(hObject, handles);

% stores ID of last axis window 
% (0 means that no axis was clicked yet)
handles.last_axis_id = 0;
handles.figure1 = hObject;
% Update handles structure
guidata(hObject, handles);
%uiwait(gcf);
% uiwait(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = SliceBrowser_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%varargout{1} = {handles.FirstP(1)};
% Get default command line output from handles structure
% varargout{1} = [handles.FirstP(2) handles.FirstP(1) handles.FirstP(3); handles.SecondP(2) handles.SecondP(1) handles.SecondP(3)];


% --- Executes on mouse press over axes background.
function Subplot1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the XY slice

%disp('Subplot1:BtnDown');

pt=get(gca,'currentpoint');
xpos=round(pt(1,2)); ypos=round(pt(1,1));
zpos = handles.pointer3dt(3);
tpos = handles.pointer3dt(4);

if (handles.IsSelected == false)
    handles.FirstP = [ypos xpos zpos];
    handles.IsSelected = true;
    handles.HasLine = false;
else
    handles.SecondP = [ypos xpos zpos];
    handles.IsSelected = false;
    handles.HasLine = true;
end

handles.pointer3dt = [xpos ypos zpos tpos];
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
plot3slices(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% --- Plots all 3 slices XY, YZ, XZ into 3 subplots
function [sp1] = plot3slices(hObject, handles)
% pointer3d     3D coordinates in volume matrix (integers)

handles.pointer3dt;
size(handles.volume);
value3dt = handles.volume(handles.pointer3dt(1), handles.pointer3dt(2), handles.pointer3dt(3), handles.pointer3dt(4));

text_str = ['[X:' int2str(handles.pointer3dt(1)) ...
           ', Y:' int2str(handles.pointer3dt(2)) ...
           ', Z:' int2str(handles.pointer3dt(3)) ...
           ', Time:' int2str(handles.pointer3dt(4)) '/' int2str(handles.vol_sz(4)) ...
           '], value:' num2str(value3dt)];
set(handles.pointer3d_info, 'String', text_str);
guidata(hObject, handles);

sliceXY = squeeze(handles.volume(:,:,handles.pointer3dt(3),handles.pointer3dt(4)));

sp1 = subplot(1,1,1);

%colorbar;
imagesc(sliceXY);
hold on;
axis image;
title('Slice XY');
ylabel('X');xlabel('Y');
if (not(isempty(handles.contour)))
    b = bwboundaries(handles.contour(:,:,handles.pointer3dt(3),handles.pointer3dt(4)));

    for k = 1:numel(b)
         plot(b{k}(:,2), b{k}(:,1), 'r', 'Linewidth', 2)
    end
end
if (handles.HasLine) 
    if (handles.pointer3dt(3) == handles.FirstP(3))
        line([handles.FirstP(1) handles.SecondP(1)], [handles.FirstP(2) handles.SecondP(2)]);
    end
end
%line([handles.pointer3dt(2) handles.pointer3dt(2)], [0 size(handles.volume,1)]);
%line([0 size(handles.volume,2)], [handles.pointer3dt(1) handles.pointer3dt(1)]);
%set(allchild(gca),'ButtonDownFcn',@Subplot1_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','SliceBrowser(''Subplot1_ButtonDownFcn'',gca,[],guidata(gcbo))');
guidata(hObject, handles);

function pointer3d_out = clipointer3d(pointer3d_in,vol_size)
pointer3d_out = pointer3d_in;
for p_id=1:4
    if (pointer3d_in(p_id) > vol_size(p_id))
        pointer3d_out(p_id) = vol_size(p_id);
    end;
    if (pointer3d_in(p_id) < 1)
        pointer3d_out(p_id) = 1;
    end;
end;

% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiresume(handles.figure1);
% delete(gcf);
% 
% delete(hObject);
% R = SliceBrowser(MRI);
R=[handles.FirstP(2) handles.FirstP(1) handles.FirstP(3); handles.SecondP(2) handles.SecondP(1) handles.SecondP(3)];
VOI= CreateVOI(handles.volume, handles.Spacing, R, single(1.35));
[LSL1Result,CAL1Result] = tumorsegmentation_L1_new(VOI.MRI, VOI.Seeds,handles.Spacing,single(1.0), single(0.2));

% 5% of average radius Offset invards
volume = sum(LSL1Result(:)) * handles.Spacing(1) * handles.Spacing(2) * handles.Spacing(3);
% assignin('base','volume',volume);
avgrad = (3*volume / 4*pi) ^ (1/3);
offset = avgrad * 0.05;%offset_ratio
Contour = false(size(LSL1Result));
for i=1:size(Contour,3)
    b = bwboundaries(LSL1Result(:,:,i));
    for k = 1:numel(b)
        for n=1:size(b{k},1);
            Contour(b{k}(n,1), b{k}(n,2),i) = true;
        end
    end
end
CI_SDM = bwdistsc(Contour, handles.Spacing);
CI_SDM(LSL1Result==1) = -CI_SDM(LSL1Result==1);
CI_SDM = CI_SDM + offset;
OffsetResult = CI_SDM < 0;

%  Calculate OTSU Treshold
NI = double(VOI.MRI(OffsetResult));
[nf,nxi] = ksdensity(NI);
minNI = min(NI);
widthNI = max(NI-min(NI));
NNI = (NI - min(NI)) / widthNI;
%level = otsu(NNI);
level = graythresh(NNI);
otsu_treshold = (level * widthNI) + minNI;

nec_area = single(sum(NI<otsu_treshold));
enh_area = single(sum(NI>otsu_treshold));

area_ratio = single(0.25);

for k=min(NI):max(NI)
    if single(sum(NI < k)) >= single(nec_area)*area_ratio
        break
    end
end
nec_treshold = k;
for k=min(NI):max(NI)
    if single(sum(NI > k)) <= single(enh_area)*area_ratio
        break
    end
end
enh_treshold = k;


Seeds = int8(zeros(size(VOI.MRI)));
Seeds((VOI.MRI < uint16(nec_treshold)) & OffsetResult) = 2;
Seeds((VOI.MRI > uint16(enh_treshold)) & OffsetResult) = 1;
%Seeds(not(VOI.CI)) = -1; %Mask healthy tissue
Seeds(not(OffsetResult)) = -1; %Mask healthy tissue
%Seeds(not(VOI.CI)) = 0;
%Seeds(VOI.Seeds==2) = 3; %3 label segmentation

[CAMap,CAStr] = castrength(VOI.MRI, Seeds, handles.Spacing, single(1.0));

LSL1Result(CAMap==2) = 2;

fullLSL1Result = uint8(zeros(size(handles.volume)));
%fullLSL1Result(VOI.Start(1):VOI.End(1),VOI.Start(2):VOI.End(2),VOI.Start(3):VOI.End(3)) = (LSL1Result==1);
fullLSL1Result(VOI.Start(1):VOI.End(1),VOI.Start(2):VOI.End(2),VOI.Start(3):VOI.End(3)) = LSL1Result;
nii = handles.nii;
nii.img = uint16(fullLSL1Result);
nii.hdr.dime.datatype = 512;
nii.hdr.dime.bitpix = 16;
% assignin('base','fullLSL1Result',fullLSL1Result);
% assignin('base','volume',handles.volume);
if (exist('LSL1Result')&&~isempty(LSL1Result))
    contour =fullLSL1Result ;
    if (ndims(contour) ~= ndims(handles.volume))
        error('Input volume and contour must have same dimensions.');
    end
    if (size(contour) ~= size(handles.volume))
        error('Input volume and contour must be in same size.');
    end
end
handles.contour = not(contour == 0);

% set main wnd title
set(gcf, 'Name', 'Slice Viewer')

% init 3D pointer
vol_sz = size(handles.volume); 
if (ndims(handles.volume) == 3)
    vol_sz(4) = 1;
end;
pointer3dt = floor(vol_sz/2)+1;
handles.pointer3dt = pointer3dt;
handles.vol_sz = vol_sz;

% handles.IsSelected = false;
% handles.HasLine = false;

% handles.FirstP = [];
% handles.SecondP = [];
assignin('base','TumorContour',handles.contour);
plot3slices(hObject, handles);

% stores ID of last axis window 
% (0 means that no axis was clicked yet)
handles.last_axis_id = 0;
handles.figure1 = hObject;
% Update handles structure
guidata(hObject, handles);
% R = SliceBrowser(MRI,LSL1Result);
%Result(CAMap==1) = 1;
%Result(CAMap==2) = 2;
%nii.img = single(fullLSL1Result);
save_untouch_nii(nii, [handles.pname 's' handles.fname]);
% SliceBrowser(MRI,fullLSL1Result);


% --- Executes on button press in btnPrevSlice.
function btnPrevSlice_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrevSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.pointer3dt(3) = handles.pointer3dt(3) - 1;
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
handles.IsSelected = false;
handles.HasLine = false;
plot3slices(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btnNextSlice.
function btnNextSlice_Callback(hObject, eventdata, handles)
% hObject    handle to btnNextSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.pointer3dt(3) = handles.pointer3dt(3) + 1;
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
handles.IsSelected = false;
handles.HasLine = false;

plot3slices(hObject, handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes when user attempts to close SliceBrowserFigure.
function SliceBrowserFigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to SliceBrowserFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object creation, after setting all properties.
function Subplot1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Subplot1


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
