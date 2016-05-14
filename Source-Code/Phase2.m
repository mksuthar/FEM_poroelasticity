function varargout = Phase2(varargin)
% PHASE2 MATLAB code for Phase2.fig
%      PHASE2, by itself, creates a new PHASE2 or raises the existing
%      singleton*.
%
%      H = PHASE2 returns the handle to a new PHASE2 or the handle to
%      the existing singleton*.
%
%      PHASE2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE2.M with the given input arguments.
%
%      PHASE2('Property','Value',...) creates a new PHASE2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Phase2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Phase2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Phase2

% Last Modified by GUIDE v2.5 03-Apr-2016 19:47:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Phase2_OpeningFcn, ...
                   'gui_OutputFcn',  @Phase2_OutputFcn, ...
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


% --- Executes just before Phase2 is made visible.
function Phase2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Phase2 (see VARARGIN)

% Choose default command line output for Phase2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Phase2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Phase2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_textdisp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_textdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_textdisp as text
%        str2double(get(hObject,'String')) returns contents of edit_textdisp as a double


% --- Executes during object creation, after setting all properties.
function edit_textdisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_textdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_browse.
function pushbutton_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.text3.String = 'Reading and Organizing (>1MB File)';
handles.text3.ForegroundColor = 'r';
[filename, pathname] = uigetfile({'*.msh';'*.*'},'File Selector');
handles.edit_textdisp.String = [pathname,filename]
fid = fopen([pathname,filename]);
C = textscan(fid,'%s','Delimiter','\n');
rowcount=str2num(C{1}{5});
rowcountnum = rowcount + 4;
nodess = dlmread([pathname,filename], ' ', [5 0 rowcountnum 3]);
fclose(fid);
save('nodes.mat','nodess')
line_begin = find(not(cellfun('isempty', strfind(C{1}, '$Elements'))));
line_end = find(not(cellfun('isempty', strfind(C{1}, '$EndElements'))));
delete C

g = extract_mesh2_elem(filename,line_begin, line_end);
M = [g{1}, g{2}, g{3}, g{4}, g{5}, g{6}, g{7}, g{8}, g{9}, g{10}, g{11}, g{12}, g{13}, g{14}];
og = unique(g{4});
delete g
bM_name = cell(1,length(og));
bM = struct ;
for i=1:length(og)
    bM.(['b', num2str(og(i))]) = [];
    bM_name{i} = ['b', num2str(og(i))];
end

for i = 1:length(M(:,1))
 nodes=M(i,4+M(1,3):end);
 nodes(find(isnan(nodes) == 1)) = [];
 phys_num = M(i,4);
 bM.(['b', num2str(phys_num)]) = [bM.(['b', num2str(phys_num)]); M(i,1), nodes ];
end

maxS = '';
maxN = 0;
for i=1:length(og)
    k = bM.(['b', num2str(og(i))]);
    if length(k(:)) > maxN
        maxS = ['b', num2str(og(i))];
        maxN = length(k(:));
    end
end

save('elements.mat','bM')
N = nodess;
C = bM.(maxS);
delete nodess 
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

order =[ 2 6 3 7 4 8 5 9 2];
plot (handles.axes1,N(C(1,order),ix),N(C(1,order),iy), '-ok')
hold on
for i=2:length(C(2:end,1))
    plot(handles.axes1,N(C(i,order),ix),N(C(i,order),iy),'-ok')
    plot(handles.axes1,N(C(i,10),ix),N(C(i,10),iy),'ok')
    plot(handles.axes1,N(C(i,[2 3 4 5]),ix),N(C(i,[2 3 4 5]),iy),'^k')
end
hold off
xlabel('\bfX')
xlim([min(N(C(:,2:8),ix))-10 max(N(C(:,2:8),ix))+10])
ylim([min(N(C(:,2:8),iy))-10 max(N(C(:,2:8),iy))+10])
ylabel('\bfY')
title('Geometry of the problem (circle = displacment node, triangle = pressure node)')
grid on

handles.text3.String = 'Ready !!';
handles.text5.String = ['Total Nodes: ', num2str(length(N(:,1)))];
handles.text6.String = ['Total Elements: ', num2str(length(C(:,1)))];
handles.pushbutton_continue.Enable = 'on';

% figure
% C = bM.(bM_name{5});
% plot (N(C(1,:),ix),N(C(1,:),iy), 'or','MarkerSize',2)
% hold on
% for i=2:length(C(1:end,1))
%     plot(N(C(i,:),ix),N(C(i,:),iy),'or','MarkerSize',2)
% end
% hold off
% xlabel('\bfX')
% ylabel('\bfY')
% title('Geometry of the problem (circle = displacment node, triangle = pressure node)')
% grid on







% --- Executes on button press in pushbutton_continue.
function pushbutton_continue_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close (gcf)
run 'Phase4.m'
