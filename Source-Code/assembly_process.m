clc, clear, close all

load 'matlab_kk.mat';

tic

% Logic for assembly 
filename = 'connect.csv';
filename2 = 'coords.csv';
C = csvread(filename);
N = csvread(filename2);

num_elements = length(C(:,1));
num_nodes = length(N(:,1));

node_1s = csvread('nodeset1.csv');
node_2s = csvread('nodeset2.csv');
node_3s = csvread('nodeset3.csv');
node_4s = csvread('nodeset4.csv');

% ASSEMBLY 

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

mu = 1;
nu = 0.3;
alpha = 1;

Ey = 2 * mu * (1 + nu);
Ey = 2e10;
c11 = Ey * (1 - nu * nu) / ((1 + nu) * (1 - nu - 2 * nu * nu))
c12 = Ey * nu / (1 - nu - 2 * nu * nu)
c66 = Ey / (2 * (1 + nu))

Cmat = [c11, c12, 0; c12, c11, 0; 0, 0, c66];

% Cmat = [1.12221e7, 3.70328e6, 0; ...
%      3.70329e6, 1.12221e7, 0; ...
%      0,0,3.7594e6];
 
alpha = 1;

for i=1:num_elements
    nodes = C(i,:);
    x_vector = N(nodes,1);
    y_vector = N(nodes,2);
    
    [ K, Q, H ] = FEM_poroelasticity(x_vector', y_vector', Cmat, alpha,1,0.3);
    
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
[GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'uy', node_3s, 0.*node_3s);
[GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'ux', node_4s, 0.*node_4s);

node_1s_p = intersect(node_1s,C(:,1:4));
node_2s_p = intersect(node_2s,C(:,1:4));

[GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'p', node_1s_p, 5 + 0.*node_1s);
[GLOBAL_K, GLOBAL_F] = FEM_poroelasticity_BC(GLOBAL_K, GLOBAL_F,dof_map, 'd', 'p', node_2s_p, 1 + 0.*node_2s);

U = GLOBAL_K \ GLOBAL_F;
x = gmres(GLOBAL_K,GLOBAL_F);

disp_dofs = dof_map(:,1:2);
disp_dofs=reshape(disp_dofs',length(disp_dofs(:)),[]);
disp_dofs(find(isnan(disp_dofs) == 1)) = -1;   
displacment = U(disp_dofs);
displacment = reshape(displacment(:),2,[])';
deformed_position = N + displacment;


pres_dof = dof_map(:,3);
pres_dof = pres_dof(find(isnan(pres_dof) == 0));
pressure = U(pres_dof);
p_index = find(isnan(dof_map(:,3)) == 0);

X = N(:,1);
Y = N(:,2);
j = 0:0.05:11;
[xq,yq] = meshgrid(j, j);

xdis = griddata(X,Y,displacment(:,1),xq,yq,'cubic');
ydis = griddata(X,Y,displacment(:,2),xq,yq,'cubic');
dis_mag = griddata(X,Y,sqrt(displacment(:,1).^2 + displacment(:,2).^2),xq,yq,'cubic');
pres_all = griddata(N(p_index,1),N(p_index,2),pressure,xq,yq,'cubic');


% Compute stress
[ s, x_loc, y_loc ] = compute_stress_dispacment(N,C,Cmat, displacment,deformed_position)

Sigma_x = griddata(x_loc', y_loc',s(1,:)',xq,yq,'cubic');
Sigma_y = griddata(x_loc', y_loc',s(2,:)',xq,yq,'cubic');
Sigma_xy = griddata(x_loc', y_loc',s(3,:)',xq,yq,'cubic');

plot_patches(N,C,deformed_position,xq,yq,xdis,ydis,dis_mag, pres_all,Sigma_x,Sigma_y, Sigma_xy)

