function [x, W, its, flag, error_res, error_orth, error_innerorth] = gmres_sstep(A, x0, b, s, m, maxiter, basis_info, tolres, tolH, output_res, usetolH)

% Solves the linear system Ax = b
% using the restarted s-step Generalized Minimal residual ( GMRES ) method
% with monomial/Newton/Chebyshev basis.
% Currently uses ||Ax-b|| <= tolres*(||A||*||x|| + ||b||) to check for
% convergence in every s iteratiions.
%
% input   A            REAL nonsymmetric positive definite matrix
%         x0           REAL initial guess vector
%         b            REAL right hand side vector
%         s            INTEGER the block size
%         m            INTEGER: the number of restarted step
%         maxiter      INTEGER: the number of the maximal iterations
%         basis_info   Structure basis_info.type = 'monomial' or 'newton' or 'chebyshev'
%         tol          REAL error tolerance
%         tolres       REAL the backward error tolerance
%         tolH         REAL error tolerance
%         output_res   INTEGER: 1 = outputing error_res, error_orth, and error_innerorth
%         usetolH      INTEGER: 1 = using tolH as a additional stopping criterion 
%        
%
% output  x            REAL solution vector
%         W            REAL the basis
%         its          INTEGER number of (inner) iterations performed
%         flag         INTEGER: 0 = solution found to tolres
%                               1 = the tolres cannot be reached at the key dimension
%                               2 = error: tolH is too small
%         error_res    REAL norm of relative backward error
%         error_orth   REAL condition number of the basis for each iteration
%         error_inner  REAL condition number of the sub basis for each iteration

x = x0;
flag = 0;
its = 0;

normestA = norm(A, 'fro');
normb = norm(b);

% Compute initial residual 
r = b-A*x;

% initialize workspace
[n, ~] = size(A);
ms = m/s;
W = zeros(n, m+1);
B = zeros(n, s+1);
V = zeros(n, m+1);
H = zeros(m+1, m);
R = zeros(m+1, m+1);
cs = zeros(m,1);
sn = zeros(m,1);
e1 = zeros(n,1);
e1(1) = 1.0;

% Output the norm of backward error of each iteration.
error_res = zeros(n, 1);
error_orth = zeros(n, 1);
error_innerorth = zeros(n, 1);
if output_res == 1
    error_res(1) = norm(r)/(normb + normestA*norm(x));
end

% Test convergence
if (norm(r) < tolres*(normb + normestA*norm(x)))
    return;
end

% Compute/set basis parameters
[alp, bet, gam, ~] = basisparams(s, A, basis_info);

% Store the basis parameters used for output
basis_info.alp = alp;
basis_info.bet = bet;
basis_info.gam = gam;



for iter = 1:maxiter
    % Begain the GMRES iteration.
    r = b - A*x;
    V(:,1) = r/norm( r );
    svec = norm( r )*e1;
    for i = 1:ceil(ms)
        sres = min(m - (i-1)*s, s);
    	is = min(i*s, m);
        its = its + sres;
    
        % Build the Krylov basis
        B = computeBasis(A, V(:, (i-1)*s+1), sres, basis_info);
        W(:, (i-1)*s+1:is) = B(:, 1:sres);
       
        % Orthogonalize AB(:, 1:s) against V(:, 1:(i-1)*s+1).
        B(:, 1:sres) = A*B(:, 1:sres);
        [V(:, (i-1)*s+2:is+1), R(1:is+1, (i-1)*s+2:is+1)] = ortho_against(sres, (i-1)*s+1, B(:, 1:sres), V(:, 1:(i-1)*s+1));
    
        % Generate H from R.
        H(1:is+1, (i-1)*s+1:is) = R(1:is+1, (i-1)*s+2:is+1);
    
        % Compute QR factorization of H by Givens rotations
        % for solving the least square problem.
        for j = (i-1)*s+1:is
            H(:, j) = apply_givens(cs(1:j-1), sn(1:j-1), H(:, j));
            [cs(j), sn(j)] = rotmat(H(j, j), H(j+1, j));
            temp = cs(j)*svec(j);
            svec(j+1) = -sn(j)*svec(j);
            svec(j) = temp;
            H(j, j) = cs(j)*H(j, j) + sn(j)*H(j+1, j);
            H(j+1, j) = 0.0; 
        end
    
        % Output the norm of backward error of each iteration.
        if output_res == 1
            for j = (i-1)*s+1:is
                y = H(1:j, 1:j)\svec(1:j);
                addvec = W(:, 1:j)*y;
                xtemp = x + addvec;
                error_res((iter-1)*m + j) = norm(b - A*xtemp)/(normb + normestA*norm(xtemp));
            end
            Wtemp = W;
            for j = 1:is
                Wtemp(:, j) = Wtemp(:, j)/norm(Wtemp(:, j), 2);
            end
            for j = (i-1)*s+1:is
                error_orth((iter-1)*m + j) = cond(Wtemp(:, 1:j));
            end
            error_innerorth((iter-1)*ceil(ms) + i) = cond(Wtemp(:, (i-1)*s+1:is));
        
        end
    
        % Test if it reaches the "key dimension".
        if usetolH == 1
            normAB = 0;
            for j = (i-1)*s+1:is
                normAB = sqrt(normAB^2 + norm(B(:, j-(i-1)*s))^2);
                if abs(H(j, j)) <= tolH*normAB
                    y = H(1:j-1, 1:j-1) \ svec(1:j-1);
                    addvec = W(:, 1:j-1)*y(1:j-1);
                    xtemp = x + addvec;
                    error_resold = norm(b - A*xtemp)/(normb + normestA*norm(xtemp));
                    y = H(1:j, 1:j) \ svec(1:j);
                    addvec = W(:, 1:j)*y;
                    x = x + addvec;
                    flag = 1;
                    its = j;
                    error_resnew = norm(b - A*x)/(normb + normestA*norm(x));
                    if error_resold <= error_resnew 
                        its = j-1;
                        x = xtemp;
                    end            
                    return;
                end
            end
        end
    
        % Solve the least square problem.
        y = H(1:is, 1:is)\svec(1:is);
        addvec = W(:, 1:is)*y;
        xtemp = x + addvec; 
    
        % Test convergence.
        if norm(A*xtemp - b) <= tolres*(normestA*norm(xtemp)+normb)
            x = xtemp;
            return;
        end
        
    end
    x = xtemp;
end

x = xtemp;
flag = 2;
end



function [Q, R] = ortho_against(m, mQ, A, Q0)
% BCGSI+ with Householder qr.
Rtemp = Q0'*A;
A = A - Q0*Rtemp;
[A, R1] = qr(A, 0);
R(1:mQ, :) = Q0'*A;
A = A - Q0*R(1:mQ, :);
[Q, R(mQ+1:mQ+m, :)] = qr(A, 0);
R(1:mQ, :) = Rtemp + R(1:mQ, :)*R1;
R(mQ+1:mQ+m, :) = R(mQ+1:mQ+m, :)*R1;
end


function h = apply_givens(cs, sn, h)
n = length(cs);
for k = 1:n
    temp = cs(k)*h(k) + sn(k)*h(k+1);
    h(k+1) = -sn(k)*h(k) + cs(k)*h(k+1);
    h(k) = temp;
end
end


function [c, s] = rotmat(a, b)
% Compute the Givens rotation matrix parameters for a and b.
if (b == 0.0)
    c = 1.0;
    s = 0.0;
elseif (abs(b) > abs(a))
    temp = a/b;
    s = 1.0/sqrt(1.0 + temp^2);
    c = temp*s;
else
    temp = b/a;
    c = 1.0/sqrt(1.0 + temp^2);
    s = temp*c;
end
end
