function [ K, Q, H ] = FEM_poroelasticity(x_vector, y_vector, C, alpha,k,mu)
% Magic !
tic
% x_vect = [0, 4, 4, 1.2];
% y_vect = [0 0.2 1 0.78];
% x_vect =[0.1000    0.1248    0.1684    0.1827    0.1124    0.1466    0.1755    0.1413    0.1440];
% y_vect = [-0.3500   -0.2255   -0.2101   -0.2541   -0.2878   -0.2178   -0.2321   -0.3020   -0.2599];
ip = numerically_int('Q93x3');
m2d = [1 1 0]';
K = compute_K_Matrix(x_vector,y_vector,ip,C,'Q9',9);
Q = compute_Q_Matrix(x_vector,y_vector,ip,'Q9',9,'Q4',4,m2d, alpha);
H = compute_H_Matrix(x_vector, y_vector, ip, 4, 'Q4',k,mu);
end

function [K] = compute_K_Matrix(x_vect,y_vect,int_points,C,type_elem,df)
% This works, dont fuck with this
    K = 0;
    for i=1:length(int_points(:,1))
        [N, dN, Nou] = shape_function(type_elem,int_points(i,1),int_points(i,2));
        [dN_xy, detJ] = compute_inv_jac(x_vect,y_vect, dN,df);
        B = compute_B_2d_matrix(df,dN_xy(:,1),dN_xy(:,2));
        K = K + int_points(i,3)*abs(detJ)*transpose(B)*C*B;
    end
end

function [Q] = compute_Q_Matrix(x_vect, y_vect, int_points, telem_u, dfu, telem_p, dfp, m, alpha)
    % Q = int(B*m*a*Np)
    Q = 0;
    for i=1:length(int_points(:,1))
        [Nu, dNu,Nou] = shape_function(telem_u,int_points(i,1),int_points(i,2));
        [Np, dNp,Nop] = shape_function(telem_p,int_points(i,1),int_points(i,2));
        [dNu_xy, detJu] = compute_inv_jac(x_vect,y_vect, dNu,dfu);
        B = compute_B_2d_matrix(dfu,dNu_xy(:,1),dNu_xy(:,2));
        Q = Q + int_points(i,3)*abs(detJu)*transpose(B)*m*Np';
    end 
end

function [H] = compute_H_Matrix(x_vect, y_vect, int_points, dfp, telem_p,k,mu)
    H = 0;
    for i=1:length(int_points(:,1))
        [Np, dNp,Nop] = shape_function(telem_p,int_points(i,1),int_points(i,2));
        [dNp_xy, detJp] = compute_inv_jac(x_vect(1:4),y_vect(1:4),dNp,dfp);
        H = H + int_points(i,3)*abs(detJp)*transpose(dNp_xy')* (k ./ mu) * dNp_xy'; 
    end    
end

function [S] = compute_S_Matrix(x_vect, y_vect, int_points, dfp, telem_p,a,n,ks,kf)
    S = 0;
    for i=1:length(int_points(:,1))
        [Np, dNp,Nop] = shape_function(telem_p,int_points(i,1),int_points(i,2));
        [dNp_xy, detJp] = compute_inv_jac(x_vect(1:4),y_vect(1:4), dNp,dfp);
        S = S + int_points(i,3)*abs(detJp)*transpose(Np')*((a-n)/ks + n/Kf) * Np'; 
    end
end

function [Fu] = compute_Fu_Matrix()

end

function [Fp] = compute_Fp_Matrix()

end

% function [K,F] = apply_D_bc(nodes, values, bc_type,K,F)
% 
%     switch bc_type
%         case 'u_x'
%             df_index = 0;
%         case 'u_y'
%             df_index = 1;
%         case 'p'
%             df_index = 2;
%         otherwise
%             disp('Error')
%     end
%     
%     for i=1:length(nodes)
%         rows = dof_map(nodes, df_index);
%         K(rows) = 0;
%         K(rows, rows) = 1;
%         F(rows) = values[i];
%     end
%     
% end

function [p] = numerically_int(type)
    switch type

        case 'Q42x2'
            p = [-0.577350269189626, -0.577350269189626, 1;...
                 0.577350269189626, -0.577350269189626, 1; ...
                 -0.577350269189626, 0.577350269189626, 1;...
                 0.577350269189626, 0.577350269189626, 1];
        case 'Q93x3'
            p = [-0.774596669241483, -0.774596669241483, 0.308641975308642; ...
                                  0, -0.774596669241483, 0.493827160493827; ...
                  0.774596669241483, -0.774596669241483, 0.3086419753086423; ...
                 -0.774596669241483,                  0, 0.493827160493827; ...
                                  0,                  0, 0.790123456790123; ...
                  0.774596669241483,                  0, 0.493827160493827; ...
                 -0.774596669241483,  0.774596669241483, 0.308641975308642; ...
                                  0,  0.774596669241483, 0.493827160493827; ...
                  0.774596669241483,  0.774596669241483, 0.308641975308642];
        otherwise
            disp('Integration Method')
    end   
end


function [dN_xy, det_J] = compute_inv_jac(x_vect, y_vect, dN,df)
    % Compute inverse jacobian 
    dx_ds = sum(x_vect .* dN(:,1)');
    dy_dn = sum(y_vect .* dN(:,2)');
    dy_ds = sum(y_vect .* dN(:,1)');
    dx_dn = sum(x_vect .* dN(:,2)');
    det_J = dx_ds*dy_dn - dy_ds*dx_dn;
    J_inv = (1/det_J).*[dy_dn, -dy_ds;-dx_dn,dx_ds];
    dN_xy = zeros(df,2);
    for i=1:df
        dN_xy(i,:) = (J_inv * [dN(i,1); dN(i,2)])';
    end
end


function [B] = compute_B_2d_matrix(df, N_xs, N_ys)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE B MATRIX
%                                              
%        |   N_1,x      0   N_2,x       0   ... |
%  B  =  |       0  N_1,y       0   N_2,y   ... |
%        |   N_1,y  N_2,x   N_2,y   N_2,x   ... |
%         -                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%B matrix
    
B = zeros(3,df*2);
B(1,1:2:2*df) = N_xs;
B(2,2:2:2*df) = N_ys;
B(3,1:2:2*df) = N_ys;
B(3,2:2:2*df) = N_xs;

end

function [N,dN,Norg] = shape_function(type,s,n)

% N = [N1; N2; ... N9]
%
% dN = [ dN1/ds dN1/dn
%        dN2/ds dN2/dn
%        dN3/ds dN3/dn .....]

    switch type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Quadrangle9:
%
%        n
%        |
%  4-----7-----3 
%  |     |     | 
%  |     |     | 
%  8     9     6 -----> s
%  |           | 
%  |           | 
%  1-----5-----2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Q9'
        % N = [N1; N2; ... N9]
        N = [  1/4 * (s^2 - s) * (n^2 - n); ...
               1/4 * (s^2 + s) * (n^2 - n); ...
               1/4 * (s^2 + s) * (n^2 + n); ...
               1/4 * (s^2 - s) * (n^2 + n); ...
               1/2 * (1 - s^2) * (n^2 - n); ...
               1/2 * (s^2 + s) * (1 - n^2); ...
               1/2 * (1 - s^2) * (n^2 + n); ...
               1/2 * (s^2 - s) * (1 - n^2); ...
               (1 - s^2) * (1 - n^2)];
           
        Norg = [N(1)*eye(2), N(2)*eye(2), N(3)*eye(2), N(4)*eye(2), ...
                N(5)*eye(2), N(6)*eye(2), N(7)*eye(2), N(8)*eye(2),...
                N(9)*eye(2)];

        dN = [ -(s/2 - 1/4)*(- n^2 + n), -(2*n - 1)*(- s^2/4 + s/4); ...
               -(s/2 + 1/4)*(- n^2 + n),    (2*n - 1)*(s^2/4 + s/4); ...
                  (s/2 + 1/4)*(n^2 + n),    (2*n + 1)*(s^2/4 + s/4); ...
                  (s/2 - 1/4)*(n^2 + n), -(2*n + 1)*(- s^2/4 + s/4); ...
                          s*(- n^2 + n),   -(2*n - 1)*(s^2/2 - 1/2); ...
                   -(n^2 - 1)*(s + 1/2),         -2*n*(s^2/2 + s/2); ...
                           -s*(n^2 + n),   -(2*n + 1)*(s^2/2 - 1/2); ...
                   -(n^2 - 1)*(s - 1/2),        2*n*(- s^2/2 + s/2); ...
                          2*s*(n^2 - 1),              2*n*(s^2 - 1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Quadrangle4:
%
%        n
%        |
%  4-----------3 
%  |           | 
%  |           | 
%  |           | -----> s
%  |           | 
%  |           | 
%  1-----------2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Q4'
        N = [1/4 * (1-s) * (1-n);...
             1/4 * (1+s) * (1-n);...
             1/4 * (1+s) * (1+n);...
             1/4 * (1-s) * (1+n)];
        
        Norg = [N(1)*eye(2), N(2)*eye(2), N(3)*eye(2), N(4)*eye(2)];
        
        dN = [   n/4 - 1/4,   s/4 - 1/4; ...
                 1/4 - n/4, - s/4 - 1/4; ...
                 n/4 + 1/4,   s/4 + 1/4; ...
               - n/4 - 1/4,   1/4 - s/4];
    otherwise
        warning('Unexpected Element Type')
    end
end


% 
% % Some clever shit from AMATH 442 HW to impose Dirchlet BC
% %
% % 	% modify global matrix and vector
% %     K = Aglobal(:,bnode);
% %     Aglobal(bnode,:) = 0;   
% %     Aglobal(:,bnode) = 0;
% %     Aglobal(bnode,bnode) = 1;
% %     Fglobal = Fglobal -gD*K;
% %     Fglobal(bnode) = gD;





