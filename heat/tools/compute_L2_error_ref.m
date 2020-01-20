function err = compute_L2_error_ref(fespace,values,values_ref)
% Compute L2 error with respect to reference solution
% input=
%           fespace: finite element space
%           values: vector of degrees of freedom
%           values_ref: vector of degrees of freedom of reference solution
% output=
%           err: L2 error
%

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;

n_elements = size(connectivity,1);

n_gauss = 7;
[gp,weights,~] = gauss_points2D(n_gauss);

err = 0;

if (~strcmp(fespace.mesh.type,'structured'))
    for i = 1:n_elements
        indices = connectivity(i,1:end-1);
        x1 = vertices(indices(1),1:2)';
        x2 = vertices(indices(2),1:2)';
        x3 = vertices(indices(3),1:2)';
        
        mattransf = [x2-x1 x3-x1];
        
        % transformation from parametric to physical
        transf = @(x) mattransf*x + x1;
        dettransf = abs(det(mattransf));
        
        for j = 1:n_gauss
            functions = fespace.functions(gp(:,j));
            value_in_gp = functions'*values(indices);
            value_ref_in_gp = functions'*values_ref(indices);
            err = err + dettransf*(value_in_gp-value_ref_in_gp)^2*weights(j)/2;
        end
    end
else
    [fespace,gp,weights] = add_members_structured_meshes(fespace,n_gauss);
    for i = 1:n_elements
        indices = connectivity(i,1:end-1);
        x1 = vertices(indices(1),1:2)';

        % then the triangle is in this configuration /|
        if (indices(2) == indices(1) + 1)
            for j = 1:n_gauss
                value_in_gp = fespace.transffuns1{j}*values(indices);
                value_ref_in_gp = fespace.transffuns1{j}*values_ref(indices);
                err = err + fespace.dettransf1*...
                    (value_in_gp-value_ref_in_gp)^2*weights(j)/2;
            end
        else
            for j = 1:n_gauss
                value_in_gp = fespace.transffuns2{j}*values(indices);
                value_ref_in_gp = fespace.transffuns2{j}*values_ref(indices);
                err = err + fespace.dettransf2*...
                    (value_in_gp-value_ref_in_gp)^2*weights(j)/2;
            end
        end
    end
end
% err = sqrt(err);
