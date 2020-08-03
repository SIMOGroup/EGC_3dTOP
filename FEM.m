function [c,dc]=FEM(ke,iK,jK,x,nel,elsize,edofMat,freedofs,penal,U,F,posD,iD)
%% modifiled SIMP
Emin=1e-4; 
sK = reshape(ke(:)*(Emin+(1-Emin)*x(:)'.^penal),24*24*nel,1); 
K = sparse(iK,jK,sK); K = 0.5*(K+K');
dd = full(sparse(iD,1,sK(posD))); % Diagonal
clearvars sK
%% delete DOFs of void element in K
skip = find(dd ==0);                   % Skipped DOFs
freedofs = setdiff(freedofs,skip);     % Remaining DOFs
% Preconditioner
alfa = 0;
while 1
    try
        P = ichol(K(freedofs,freedofs),struct('type','ict','droptol',0.001,'diagcomp',alfa));
        break;
    catch
        alfa = (alfa+0.0001)*2;
    end
end
U(freedofs)= pcg(K(freedofs,freedofs),F(freedofs),[],10000,P,P',U(freedofs));
% U(freedofs) = K(freedofs,freedofs)\F(freedofs,:);       % direct solver
c = F'*U; 
ce = reshape(sum((U(edofMat)*ke).*U(edofMat),2),elsize);
dc = -penal*(1-Emin)*x.^(penal-1).*ce;
dc=dc(:);