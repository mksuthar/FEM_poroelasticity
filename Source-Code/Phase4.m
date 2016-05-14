function varargout = Phase4(varargin)
% PHASE4 MATLAB code for Phase4.fig
%      PHASE4, by itself, creates a new PHASE4 or raises the existing
%      singleton*.
%
%      H = PHASE4 returns the handle to a new PHASE4 or the handle to
%      the existing singleton*.
%
%      PHASE4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE4.M with the given input arguments.
%
%      PHASE4('Property','Value',...) creates a new PHASE4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Phase4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Phase4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Phase4

% Last Modified by GUIDE v2.5 06-Apr-2016 15:21:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Phase4_OpeningFcn, ...
                   'gui_OutputFcn',  @Phase4_OutputFcn, ...
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


% --- Executes just before Phase4 is made visible.
function Phase4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Phase4 (see VARARGIN)

% Choose default command line output for Phase4
handles.output = hObject;
load nodes.mat
load elements.mat
% Update handles structure
guidata(hObject, handles);
Bnames = fieldnames(bM)
ir = 0;
maxc= 0;
for i=1:length(Bnames)
     if maxc < length(bM.(Bnames{i}))
        ir = i;
        maxc = length(bM.(Bnames{i}));
     end
end
bb = Bnames;
bb{ir} = [];
bbc = {};
for i = 1: length(bb)
    if length(bb{i}) > 0
        bbc{i} = bb{i};
    end
end
        
Bnames{ir} = [Bnames{ir},' - Surface'];
set(handles.uitable5,'data', zeros(length(bbc), 3));
set(handles.uitable5,'RowName',bbc, 'ColumnEditable', true(1, 3));
set(handles.listbox1,'String',Bnames,'Value',1)
set(handles.uitable5,'ColumnEditable', true(1, 3));
C = bM.(Bnames{1});
N = nodess;
%handles.(item_selected)
%guidata(hObject,handles)
ix = 0;
iy = 0;
if sum(N(:,2)) > 0 && sum(N(:,3)) > 0
    ix = 2;
    iy = 3;
end

if sum(N(:,3)) > 0 && sum(N(:,4)) > 0
    ix = 3;
    iy = 4;
end

if sum(N(:,2)) > 0 && sum(N(:,4)) > 0
    ix = 2;
    iy = 4;
end

plot(handles.axes2,N(C(1,[2 3 4]),ix),N(C(1,[2 3 4]),iy), '-ok')
hold on
for i=2:length(C(2:end,1))
    plot(handles.axes2,N(C(i,[2 3 4]),ix),N(C(i,[2 3 4]),iy),'ok')
end
hold off
xlabel('\bfX')
xlim([min(N(C(:,:),ix))-10 max(N(C(:,:),ix))+10])
ylim([min(N(C(:,:),iy))-10 max(N(C(:,:),iy))+10])
ylabel('\bfY')
title('Geometry of the problem (circle = displacment node, triangle = pressure node)')
grid on

% UIWAIT makes Phase4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Phase4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load nodes.mat
load elements.mat
items = get(hObject,'String');
index_selected = get(hObject,'Value');
item_selected = items{index_selected};

if length(strfind(item_selected,'Surface')) > 0
    item_selected = strrep(item_selected, ' - Surface', '')
end
set(handles.uitable5,'ColumnEditable', true(1, 3));
C = bM.(item_selected);
N = nodess;
%handles.(item_selected)
%guidata(hObject,handles)
ix = 0;
iy = 0;
if sum(N(:,2)) > 0 && sum(N(:,3)) > 0
    ix = 2;
    iy = 3;
end

if sum(N(:,3)) > 0 && sum(N(:,4)) > 0
    ix = 3;
    iy = 4;
end

if sum(N(:,2)) > 0 && sum(N(:,4)) > 0
    ix = 2;
    iy = 4;
end

plot(handles.axes2,N(C(1,[2 3 4]),ix),N(C(1,[2 3 4]),iy), '-ok')
hold on
for i=2:length(C(2:end,1))
    plot(handles.axes2,N(C(i,[2 3 4]),ix),N(C(i,[2 3 4]),iy),'ok')
end
hold off
xlabel('\bfX')
xlim([min(N(C(:,:),ix))-10 max(N(C(:,:),ix))+10])
ylim([min(N(C(:,:),iy))-10 max(N(C(:,:),iy))+10])
ylabel('\bfY')
title('Geometry of the problem (circle = displacment node, triangle = pressure node)')
grid on



% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load nodes.mat
load elements.mat

N = nodess;

% Shear modulus
shm = handles.edit7.String;
shm = str2double (shm);
handles.edit1.String = 0.00001;
% poission
ps = handles.edit5.String;
ps = str2double (ps);

gravity = handles.edit4.String;
gravityn = str2double (gravity);

% density
density = handles.edit3.String;
densityn = str2double (density);

% biot
bc = handles.edit2.String;
bcn = str2double (bc);

% 1/M
M_1 = handles.edit1.String;
M_1= str2double (M_1);

valid_biot = (bcn > 0) & (bcn < 1);
valid_M = (M_1 > 0) & (M_1 < 1);
valid_shm = (shm > 0) & (shm < 100);
valid_d = (densityn > 0) & (densityn < 5);
valid_ps = (ps > -0.5) & (ps < 0.5);
tline = [valid_biot, valid_M, valid_shm,valid_d,  valid_ps];
if sum(tline) < length(tline)
run 'Phase3a.m'
else
    parameters = [shm, ps, gravityn, densityn, bcn, M_1]
    save ('parameters.mat','parameters')

row_names = handles.uitable5.RowName;
for i=1:length(handles.uitable5.RowName)
    C = bM.(row_names{i});
    bnodes = N(C(:,[2,3, 4]));
    bpressure = zeros(1,length(bnodes(:))) + handles.uitable5.Data(i,1);
    bdispx = zeros(1,length(bnodes(:))) + handles.uitable5.Data(i,2);
    bdispy = zeros(1,length(bnodes(:))) + handles.uitable5.Data(i,3);
    csvwrite(['bc_',row_names{i},'.csv'],[bnodes(:), bpressure', bdispx', bdispy']);
end
close (gcf)
run 'Phase4a.m'
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
