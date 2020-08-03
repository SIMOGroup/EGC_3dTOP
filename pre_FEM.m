function [U,F,edofMat,freedofs,Xe,Ye,Ze,iK,jK,posD,iD]=pre_FEM(nelx,nely,nelz,dx,dy,dz,nel)
db=0; % distance from boundary
Xe=dx*[0.5:1:nelx]'; Ye=dy*[nely-0.5:-1:0]'; Ze=dz*[0.5:1:nelz]';
[xx,yy,zz]=meshgrid(Xe,Ye,Ze);
Xe=xx(:);
Ye=yy(:);
Ze=zz(:);   % centroid of element cordinates
%% Load dofs  
loadid=1+(nelx+1)*(nely+1)*nelz/2;   % upper midpoint
loaddof = 3*loadid-1; % DOFs
loaddof=loaddof(:);
% PREPARE FINITE ELEMENT ANALYSIS
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nel,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nel,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nel,1);
elnodeID=edofMat(:,[3:3:24])/3;     % element node ID
posD = find(jK == iK);              % Positions of diagonal of matrix
iD = iK(posD);
%% Fixed dofs
[iif,jf,kf] = meshgrid(db,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs on the left side
fixeddof1 = [3*fixednid-2]; % DOFs
fixednid2=(nelx+1-db)*(nely+1):(nelx+1)*(nely+1):(nelx+1)*(nely+1)*(nelz+1);   % right side
fixeddof2=[3*fixednid2-1]; % DOFs
fixeddof=[fixeddof1(:);fixeddof2'];
%% load and displacement vector
P0=-1/length(loadid);  % load
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,P0,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);




















