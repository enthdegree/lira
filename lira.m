function mtx_C = lira(v_x)
%% LIRA Lovasz Integer Relation Algorithm
% Lovasz Integer Relation Algorithm from "Finding Integer Relations in
% Polynomial Time" by Hastad, Just, Lagrias and Schnorr.
%
% Find a basis for the collection of integer vectors orthogonal to v_x.
% Paraphrasing the paper:
%
% > On input v_x in Z^n, this algorithm finds a basis of the
% > (n-1)-dimensional integer lattice L whose complement is the
% > integer-span of v_x. The algorithm and output have nice properties by
% > Theorem 6.1.
%
% Input: 
%   v_x = integer n-vector
%
% Output:
%  mtx_C = n*(n-1) integer matrix, each column an integer vector orthogonal
%          to v_x, with properties according to 
%
% Usage Example:
%   >> v_x = [1;4;2;1];
%   >> mtx_C = lira(v_x)
%
%   ans =
%
%       -1    -1     1
%        1     0     0
%       -1     1     0
%       -1    -1    -1
%
%   >> v_x'*mtx_C 
%
%   ans = 
%
%       0     0     0
%

    n_n = size(v_x,1);
    
    %% 1. Initiation.
    mtx_b = [v_x/vgcd(v_x), eye(n_n)];
    [v_bsn, mtx_mu] = gsn(mtx_b(:,1), mtx_b(:,2:end));
    n_i = 1;
    
    while(true)
        %% 2. Termination test.
        if(n_i == n_n)
            mtx_C = ((mtx_b(:, [1, 3:(n_n+1)]))^(-1))';
            mtx_C = mtx_C(:,2:end);
            return;
        end

        while(true)
            %% 3. Reduction in size.
            for n_j=n_i:(-1):0
                % (paper index "i")+1 = matrix index
                mtx_b(:,(n_i+1)+1) = mtx_b(:,(n_i+1)+1) - ...
                    round(mtx_mu((n_i+1)+1, n_j+1))*mtx_b(:,n_j+1);
                for n_nu=0:n_j
                    mtx_mu((n_i+1)+1,n_nu+1) = mtx_mu((n_i+1)+1, n_nu+1) - ...
                        round(mtx_mu((n_i+1)+1, n_j+1)) * mtx_mu(n_j+1, n_nu+1);
                end
            end
            if(v_bsn(n_i+1) <= 2*v_bsn((n_i+1)+1))
                n_i = n_i + 1;
                % GOTO 2
                break;
            end

            %% 4. Exchange step.
            v_idx = 1:(n_n+1);
            v_idx(n_i+1) = (n_i+1)+1;
            v_idx((n_i+1)+1) = n_i+1;
            mtx_b = mtx_b(:,v_idx);
            [v_bsn, mtx_mu] = gsn(mtx_b(:,1), mtx_b(:,2:end));
            if(n_i> 1)
                n_i = n_i-1;
            end
            % GOTO 3
        end
    end
end
function d = vgcd(v)
%% VGCD Greatest common divisor of a vector of integers
% Input: 
%   v = real integer vector
% Output:
%   d = greatest common divisor of all components of v
    if(isempty(v))
        d = 1;
        return;
    else
        d = abs(v(1));
        for ii=2:length(v)
            d = gcd(d, v(ii));
        end
        return;
    end
end

function [v_bsn, mtx_mu, mtx_bs] = gsn(v_b0, mtx_b)
%% GSN Compute the Gram-Schmidt numbers of a list of vectors
% Input: 
%   v_b0 = n*1 initial vector (0 vector for the usual Gram-Schmidt)
%   mtx_b = n*m Matrix, each i-th column a vector b_i
% Outputs:
%   v_bsn = 1*(m+1) vector, norms of orthogonalized vectors:
%           |b_0^\ast|^2, ..., |b_m^\ast|^2
%   mtx_mu = (m+1)*m matrix, lower-triangular with zero-diagonals,
%            coefficients to obtain b_i^\ast from v_x,mtx_b. That is:
%              b_0^\ast = v_x
%              b_i^\ast = mtx_b(:,i) - [b_0^\ast, ..., b_{i-1}^\ast]*mtx_mu(i,1:i)'
%   mtx_mu(i,j) = mu_{i,j-1}
%   mtx_bs = n*(m+1) matrix of Gram-Schmidt orthogonalized vectors:
%              b_{i-1} = mtx_bs(:,), i=1,...,m+1

    mtx_b = [v_b0, mtx_b];
    
    n_dim = size(mtx_b,1);
    n_m = size(mtx_b,2);
    
    v_bsn = zeros(1,n_m);
    mtx_bs = zeros(n_dim, n_m);
    mtx_mu = zeros(n_m, n_m);
    
    mtx_bs(:,1) = v_b0;
    v_bsn(1) = norm(mtx_bs(:,1))^2;
    for ii=2:n_m
    for jj=1:ii
        if(norm(mtx_bs(:,jj))==0)
            mtx_mu(ii,jj) = 0;
        else 
            mtx_mu(ii,jj) = ...
                (mtx_b(:,ii)'*mtx_bs(:,jj))/norm(mtx_bs(:,jj))^2;
        end
    end
    mtx_bs(:,ii) = mtx_b(:,ii) - mtx_bs*mtx_mu(ii,:)';
    v_bsn(ii) = norm(mtx_bs(:,ii))^2;
    end
    
end
