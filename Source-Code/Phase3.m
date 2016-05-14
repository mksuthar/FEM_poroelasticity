function varargout = Phase3(varargin)
% PHASE3 MATLAB code for Phase3.fig
%      PHASE3, by itself, creates a new PHASE3 or raises the existing
%      singleton*.
%
%      H = PHASE3 returns the handle to a new PHASE3 or the handle to
%      the existing singleton*.
%
%      PHASE3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE3.M with the given input arguments.
%
%      PHASE3('Property','Value',...) creates a new PHASE3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Phase3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Phase3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Phase3

% Last Modified by GUIDE v2.5 27-Mar-2016 21:29:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Phase3_OpeningFcn, ...
                   'gui_OutputFcn',  @Phase3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
load nodes.mat;
load elements.mat;

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Phase3 is made visible.
function Phase3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Phase3 (see VARARGIN)

% Choose default command line output for Phase3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Phase3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Phase3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_stress_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stress as text
%        str2double(get(hObject,'String')) returns contents of edit_stress as a double


% --- Executes during object creation, after setting all properties.
function edit_stress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sbm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sbm as text
%        str2double(get(hObject,'String')) returns contents of edit_sbm as a double



% --- Executes during object creation, after setting all properties.
function edit_sbm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fbm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fbm as text
%        str2double(get(hObject,'String')) returns contents of edit_fbm as a double


% --- Executes during object creation, after setting all properties.
function edit_fbm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_porosity_Callback(hObject, eventdata, handles)
% hObject    handle to edit_porosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_porosity as text
%        str2double(get(hObject,'String')) returns contents of edit_porosity as a double


% --- Executes during object creation, after setting all properties.
function edit_porosity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_porosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dv as text
%        str2double(get(hObject,'String')) returns contents of edit_dv as a double


% --- Executes during object creation, after setting all properties.
function edit_dv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gravity_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gravity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gravity as text
%        str2double(get(hObject,'String')) returns contents of edit_gravity as a double

% --- Executes during object creation, after setting all properties.
function edit_gravity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gravity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_density_Callback(hObject, eventdata, handles)
% hObject    handle to edit_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_density as text
%        str2double(get(hObject,'String')) returns contents of edit_density as a double


% --- Executes during object creation, after setting all properties.
function edit_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bc as text
%        str2double(get(hObject,'String')) returns contents of edit_bc as a double

% --- Executes during object creation, after setting all properties.
function edit_bc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_permeability_Callback(hObject, eventdata, handles)
% hObject    handle to edit_permeability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_permeability as text
%        str2double(get(hObject,'String')) returns contents of edit_permeability as a double


% --- Executes during object creation, after setting all properties.
function edit_permeability_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_permeability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stress = handles.edit_stress.String;
stressn = str2double (stress);
sbm = handles.edit_sbm.String;
sbmn = str2double(sbm);
fbm = handles.edit_fbm.String;
fbmn = str2double (fbm);
porosity = handles.edit_porosity.String;
porosityn = str2double (porosity);
dv = handles.edit_dv.String;
dvn = str2double (dv);
gravity = handles.edit_gravity.String;
gravityn = str2double (gravity);
density = handles.edit_density.String;
densityn = str2double (density);
bc = handles.edit_bc.String;
bcn = str2double (bc);
permeability = handles.edit_permeability.String;
permeabilityn = str2double (permeability);
if stressn <= 0 | fbmn <= 0 | sbmn <= 0 | porosityn <= 0 | dvn <= 0....
    | gravityn <= 0 | densityn <= 0 | bcn <= 0 | permeability <= 0  
run 'Phase3a.m'
else
    close (gcf)
    run 'Phase4.m'
end
