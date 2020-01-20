function V = POD_Xh(S, Xh, err_POD)
% Function that returns an Xh-orthonormal reduced basis for a given set of snapshots and a
% given POD error
%
% Inputs : 
%           S : Set of snapshots
%           Xh : Matrix norm
%           err_POD : POD error
%
% Outputs : 
%           V : Reduced basis
%           

% Form Correlation matrix
C = S' * Xh * S;

% Get eigen values and eigen vectors of the correlation matrix
[eig_vectors, eig_values] = eig(C);

% Sort the eigen values and eigenvectors in decreasing order
eig_values = diag(eig_values);
[eig_values, idx] = sort(eig_values, 'descend');
eig_vectors = eig_vectors(:, idx);

% Take the sum of the eigenvalues
sum_eig_values = sum(eig_values);

sum_reduced_eig_values = 0;
N = 0;
V = [];
while sum_reduced_eig_values / sum_eig_values < 1 - err_POD^2
    N = N + 1;
    sum_reduced_eig_values = sum_reduced_eig_values + eig_values(N);
    
    V = horzcat(V, (1. / sqrt(eig_values(N))) * S * eig_vectors(:, N));
end

return


