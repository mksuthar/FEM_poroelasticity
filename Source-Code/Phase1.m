function varargout = Phase1(varargin)
% PHASE1 MATLAB code for Phase1.fig
%      PHASE1, by itself, creates a new PHASE1 or raises the existing
%      singleton*.
%
%      H = PHASE1 returns the handle to a new PHASE1 or the handle to
%      the existing singleton*.
%
%      PHASE1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE1.M with the given input arguments.
%
%      PHASE1('Property','Value',...) creates a new PHASE1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Phase1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Phase1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Phase1

% Last Modified by GUIDE v2.5 25-Mar-2016 22:40:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Phase1_OpeningFcn, ...
                   'gui_OutputFcn',  @Phase1_OutputFcn, ...
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


% --- Executes just before Phase1 is made visible.
function Phase1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Phase1 (see VARARGIN)

% Choose default command line output for Phase1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Phase1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Phase1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_agree.
function pushbutton_agree_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_agree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete('nodes.mat')
delete('parameters.mat')
delete('elements.mat')
k = {};
c = 1;
listing = dir
for i=1:length(listing)
 if length(findstr(listing(i).name, 'bc_')) > 0
     if length(findstr(listing(i).name, '.csv')) > 0
        k{c} = listing(i).name;
        c= c+1;
     end
 end
end
for i=1:length(k)
    delete (k{i})
end
close(gcf);
run 'Phase2.m'


% --- Executes on button press in pushbutton_disagree.
function pushbutton_disagree_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_disagree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
