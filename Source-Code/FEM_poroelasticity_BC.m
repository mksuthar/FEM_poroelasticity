function [ K,F ] = FEM_poroelasticity_BC(K,F,dof_map,type, dof_type, nodes, values)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    switch type
        case 'n'
            [K,F] =neuman_BC(K,F,dof_type, nodes, values,dof_map);
        case 'd'
            [K,F] = dirchlet_BC(K,F,dof_type, nodes, values,dof_map);
    end

end

function [K,F] = neuman_BC(K,F,dof_type, nodes, values,dof_map)
    
end

function [K,F] = dirchlet_BC(K,F,dof_type, nodes, values,dof_map)
        
        num_nodes = length(nodes);

        switch dof_type
            case 'ux'
                col = 1;
            case 'uy'
                col = 2;
            case 'p'
                col = 3;
        end
        
        for i=1:num_nodes
            row_index = dof_map(nodes(i), col);
            K(row_index, :) = 0;
            K(row_index, row_index) = 1;
            F(row_index) = values(i);
        end
        
end

