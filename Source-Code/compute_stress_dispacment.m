function [s, x_loc, y_loc] = compute_stress_dispacment(N,C,Cmat,disp,deformed)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    num_elements = length(C(:,1));
    ip=numerically_int('Q93x3');
    [x_loc,y_loc] = find_xy_gauss_points(N,C,num_elements,deformed)
    s = []
    for e=1:num_elements       
        nodes = C(e,:);
        x_vect = N(nodes,1);
        y_vect = N(nodes,2);
        curr_disp = disp(C(e,:),:);
        s = [s,compute_stress(ip,x_vect,y_vect,curr_disp,Cmat)];
    end
    
end

function [x_loc,y_loc] = find_xy_gauss_points(N,C,num_elements,deformed)
    x_loc = [];
    y_loc = [];
    ip=numerically_int('Q93x3');
    for e=1:num_elements
        disp(e)
        nodes = C(e,:);
        x_vect = deformed(nodes,1);
        y_vect = deformed(nodes,2);
        for j = 1:length(ip(:,1))
            [Ns,dN,Norg] = shape_function('Q9',ip(j,1),ip(j,2));
            x_loc = [x_loc, sum(x_vect .* Ns)]; 
            y_loc = [y_loc, sum(y_vect .* Ns)];
        end       
    end
end



function S = compute_stress(ip,x_vect,y_vect,disp,Cmat)
    S = [];
    for i=1:length(ip(:,1))
            [N, dN, Nou] = shape_function('Q9',ip(i,1),ip(i,2));
            [dN_xy, detJ] = compute_inv_jac(x_vect',y_vect', dN,9);
            Bmat = compute_B_2d_matrix(9,dN_xy(:,1),dN_xy(:,2));
            disp_array = reshape(disp',[],length(disp(:)))';
            S = [S,Cmat*Bmat*disp_array];
    end
end

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


