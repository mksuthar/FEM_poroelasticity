function varargout = Phase4a(varargin)
% PHASE4A MATLAB code for Phase4a.fig
%      PHASE4A, by itself, creates a new PHASE4A or raises the existing
%      singleton*.
%
%      H = PHASE4A returns the handle to a new PHASE4A or the handle to
%      the existing singleton*.
%
%      PHASE4A('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE4A.M with the given input arguments.
%
%      PHASE4A('Property','Value',...) creates a new PHASE4A or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Phase4a_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Phase4a_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Phase4a

% Last Modified by GUIDE v2.5 05-Apr-2016 23:29:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Phase4a_OpeningFcn, ...
                   'gui_OutputFcn',  @Phase4a_OutputFcn, ...
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


% --- Executes just before Phase4a is made visible.
function Phase4a_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Phase4a (see VARARGIN)

% Choose default command line output for Phase4a
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% ASSEMBLY PROCESS ------------------------------------------------------
load nodes.mat
load elements.mat
load parameters.mat
shm = parameters(1);
ps = parameters(2);
gravityn = parameters(3);
densityn = parameters(4);
bcn = parameters(5);
M_1 = parameters(6);

% find the surface
ir = 0;
maxc= 0;
Bnames = fieldnames(bM);
for i=1:length(Bnames)
     if maxc < length(bM.(Bnames{i}))
        ir = i;
        maxc = length(bM.(Bnames{i}));
     end
end
C = bM.(Bnames{ir});
N = nodess;
C(:,1) = [];

num_elements = length(C(:,1));
num_nodes = length(N(:,1));

BCN = struct;

listing = dir
names_csv = {};
c = 1;
for i=1:length(listing)
 if length(findstr(listing(i).name, 'bc_')) > 0
     if length(findstr(listing(i).name, '.csv')) > 0
        names_csv{c} = listing(i).name;
        c= c+1;
     end
 end
end

% loop to read the nodes
for i=1:length(names_csv)
    M = csvread(names_csv{i});
    k = strrep(names_csv{i}, 'bc_', '');
    k = strrep(k, '.csv', '');
    BCN.(k) = M;
end


% ASSEMBLY
% -------------------------------------------------------------------------

dof_map = zeros(num_nodes, 3);
disp_only_nodes = C(:,5:end);
    
dof_map(disp_only_nodes(:),3) = -1;
df_total = sum((dof_map(:) == 0));

GLOBAL_K = zeros(df_total,df_total);
GLOBAL_F = zeros(df_total,1);

k = 1;    
for i=1:length(dof_map(:,1))
    for j=1:length(dof_map(i,:))
        if dof_map(i,j) == 0
            dof_map(i,j) = k;
            k = k+1;
        else
            dof_map(i,j) = nan;
        end
    end    
end

mu = shm .* 1e9;%1; % shera modulus
nu = ps;%0.3; % 
alpha = bcn; %1;

Ey = 2 * mu * (1 + nu) ;
%Ey = 2e10;
c11 = Ey * (1 - nu * nu) / ((1 + nu) * (1 - nu - 2 * nu * nu))
c12 = Ey * nu / (1 - nu - 2 * nu * nu)
c66 = Ey / (2 * (1 + nu))

Cmat = [c11, c12, 0; c12, c11, 0; 0, 0, c66];

% Cmat = [1.12221e7, 3.70328e6, 0; ...
%      3.70329e6, 1.12221e7, 0; ...
%      0,0,3.7594e6];
 
alpha = 1;

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

for i=1:num_elements
    nodes = C(i,1:end);
    x_vector = N(nodes,ix);
    y_vector = N(nodes,iy);
    
    [ K, Q, H ] = FEM_poroelasticity(x_vector', y_vector', Cmat, alpha,mu,ps);
    
    disp_dofs = dof_map(nodes,1:2);
    pres_dofs = dof_map(nodes,3);
    pres_dofs(find(isnan(pres_dofs) == 1)) = [];
    
    disp_dofs = reshape(disp_dofs',[],1);
    pres_dofs = reshape(pres_dofs',[],1);
    GLOBAL_K(disp_dofs,disp_dofs) = GLOBAL_K(disp_dofs,disp_dofs) + K;
    GLOBAL_K(disp_dofs,pres_dofs) = GLOBAL_K(disp_dofs,pres_dofs) + Q;
    GLOBAL_K(pres_dofs,disp_dofs) = GLOBAL_K(pres_dofs,disp_dofs)+ transpose(Q);
    
end
% problem.apply_essential_bc(ns3,np.zeros(len(ns3)),dof")
% [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'uy', node_3s, 0.*node_3s);
% [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'ux', node_4s, 0.*node_4s);
% 
% node_1s_p = intersect(node_1s,C(:,1:4));
% node_2s_p = intersect(node_2s,C(:,1:4));
% 
% [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'p', node_1s_p, 5 + 0.*node_1s);
% [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'p', node_2s_p, 1 + 0.*node_2s);

%%% LOOOP TO COVER
for i=1:length(names_csv)
    k = strrep(names_csv{i}, 'bc_', '');
    k = strrep(k, '.csv', '');
    M=BCN.(k);
    bc_nodes = M(:,1);
    pv = M(:,2);
    ux = M(:,3);
    uy = M(:,4);
    
    if ux(1) > -100
        [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'uy', bc_nodes, uy);
    end
    if uy(1) > -100
        [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'ux', bc_nodes, ux);
    end
    node_p = intersect(bc_nodes,C(:,1:4));
    [GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'p', node_p, pv);
end


%%%


U = GLOBAL_K \ GLOBAL_F;
x = gmres(GLOBAL_K,GLOBAL_F);

disp_dofs = dof_map(:,1:2);
disp_dofs=reshape(disp_dofs',length(disp_dofs(:)),[]);
disp_dofs(find(isnan(disp_dofs) == 1)) = -1;   
displacment = U(disp_dofs);
displacment = reshape(displacment(:),2,[])';
deformed_position = N(:,[ix, iy]) + displacment;


pres_dof = dof_map(:,3);
pres_dof = pres_dof(find(isnan(pres_dof) == 0));
pressure = U(pres_dof);
p_index = find(isnan(dof_map(:,3)) == 0);

X = N(:,ix);
Y = N(:,iy);
j = 0:0.05:11;
[xq,yq] = meshgrid(j, j);

xdis = griddata(X,Y,displacment(:,1),xq,yq,'cubic');
ydis = griddata(X,Y,displacment(:,2),xq,yq,'cubic');
dis_mag = griddata(X,Y,sqrt(displacment(:,1).^2 + displacment(:,2).^2),xq,yq,'cubic');
pres_all = griddata(N(p_index,ix),N(p_index,iy),pressure,xq,yq,'cubic');

Nz = [N(:,ix),N(:,iy)];

% Compute stress
[ s, x_loc, y_loc ] = compute_stress_dispacment(Nz,C,Cmat, displacment,deformed_position)

Sigma_x = griddata(x_loc', y_loc',s(1,:)',xq,yq,'cubic');
Sigma_y = griddata(x_loc', y_loc',s(2,:)',xq,yq,'cubic');
Sigma_xy = griddata(x_loc', y_loc',s(3,:)',xq,yq,'cubic');

plot_patches(N,C,deformed_position,xq,yq,xdis,ydis,dis_mag, pres_all,Sigma_x,Sigma_y, Sigma_xy)

% ASSEMBLY
% -------------------------------------------------------------------------
handles.text2.String = 'Finished! Results are available in PDF format, check the program directory';



% UIWAIT makes Phase4a wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Phase4a_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




