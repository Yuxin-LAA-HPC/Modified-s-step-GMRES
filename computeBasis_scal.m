% computeBasis_scal.m
% This file computes a Krylov subspace basis of size s with starting vector
% v0 and matrix A, using the polynomial basis defined by alp, beta, and
% gam.
%
% Input:
%   A : a square matrix for which we will compute a Krylov basis,
%   v0 : the vector defining the Krylov subspace vector
%   s : the dimension of the resulting basis
%   basis_info: a struct which has vectors 'alp', 'beta', and 'gam, defining the polynomial basis parameters.
%
% Output:
%   V : a matrix whose columns form a basis for the defined Krylov subspace
%
% This file is based on computeBasis.m copyed from
% https://github.com/eccarson/mixedsstep.git.
%
% Last edited by: Yuxin Ma, 2024
%



function V = computeBasis_scal(A, v0, s, basis_info)

% Get basis parameters stored in struct
alp = basis_info.alp;
bet = basis_info.bet;
gam = basis_info.gam;

% Get dimension of A
n = size(A, 1);

% Set first basis vector
V(:, 1) = v0;

% Construct the remaining basis vectors iteratively
if(s>1)
    V(:,2) = V(:,1)./gam(1);
    V(:,2) = A*V(:,2) - alp(1).*V(:,2);
    V(:, 2) = V(:, 2)/norm(V(:, 2));
    % V(:,2) = (A-alp(1).*speye(n))*V(:,1)./gam(1);
    if(s>2)
        for qq = 2:(s-1)  
            V(:,qq+1) = (A*V(:,qq) - alp(qq).*V(:,qq) - bet(qq-1)*V(:,qq-1))./gam(qq);
            V(:, qq+1) = V(:, qq+1)/norm(V(:, qq+1));
            % V(:,qq+1) = ((A-alp(qq).*speye(n))*V(:,qq) - bet(qq-1)*V(:,qq-1))./gam(qq);
        end
    end
end

end