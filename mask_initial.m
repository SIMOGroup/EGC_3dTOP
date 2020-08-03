function [mask,nc,nc2,Phi,np]=mask_initial(Lx,Ly,Lz)
%% define centre axis
db=3.5; % distance from boundary of design domain
Nx=3; 
Ny=3; 
nex=2;  Dx=(Lx-2*db)/Nx; 
ney=2;  Dy=(Ly-2*db)/Ny;
nc=2*(nex*Nx*(Ny+1)+ney*Ny*(Nx+1));     % number of end points of EGCs
XC=zeros(nc,1);
YC=zeros(nc,1);
zzc=linspace(db,Lx-db,nex*Nx+1);
yyc=linspace(db,Ly-db,ney*Ny+1);
for wj=1:Ny+1
XC([1:2:2*nex*Nx]+(wj-1)*2*nex*Nx)=zzc(1:end-1);
XC([2:2:2*nex*Nx]+(wj-1)*2*nex*Nx)=zzc(2:end);
YC([1:2*nex*Nx]+(wj-1)*2*nex*Nx)=(wj-1)*Dy+db;
end
for wi=1:Nx+1
XC([1:2*ney*Ny]+(wi-1)*2*ney*Ny+2*nex*Nx*(Ny+1))=(wi-1)*Dx+db;
YC([1:2:2*ney*Ny]+(wi-1)*2*ney*Ny+2*nex*Nx*(Ny+1))=yyc(1:end-1);
YC([2:2:2*ney*Ny]+(wi-1)*2*ney*Ny+2*nex*Nx*(Ny+1))=yyc(2:end);
end
ZC=0.5*Lz*ones(nc,1);
%%
np=24;  % number of polygonal segments 
Rmat=6.5*ones(np,nc/2);
Phi=linspace(0,2*pi,np+1)';
nc2=np*nc/2;    % number of vertices 
Ra=6*ones(nc/2,1);
mask=[XC;YC;ZC;Rmat(:);Ra]; % geometric parameters



