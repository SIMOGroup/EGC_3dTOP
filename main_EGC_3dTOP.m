%% This code was selcted from [Van-Nam Hoang, H. Nguyen-Xuan, Extruded-geometric-component-based 3D topology optimization
%% ,Computer Methods in Applied Mechanics and Engineering, in press, 2020 (Q1)].

% MMA updates are avalable on http://www.smoptit.se/,
% phi2stl for post-processing is avalable on
% http://me.eng.stonybrook.edu/~chen/index_files/phi2stl_code.html,     

clear all
close all
sc=1;
Lx=150;          Ly=50;          Lz=10;
nelx=Lx*sc;     nely=Ly*sc;     nelz=Lz*sc;
dx=Lx/nelx; dy=Ly/nely; dz=Lz/nelz;
nel = nelx*nely*nelz;                % number of elements
elsize = [nely nelx nelz];
volfrac0=0.3;
%% Large number
alp0=6;    % depend on mesh-size
%% Material properties
E=1;
nu=0.3;     % elastic modulus (only true for E=1)
penal=3;    % penalization factor
%% Initial Mask
[mask,nc,nc2,Phi,np]=mask_initial(Lx,Ly,Lz);
%% Preparation for FEM
[U,F,edofMat,freedofs,Xe,Ye,Ze,iK,jK,posD,iD]=pre_FEM(nelx,nely,nelz,dx,dy,dz,nel);
[ke]=kemefe(E,nu,sc);
%% MMA parameters 
m=1;            % number of general constraints
n=length(mask); % number of variables 
[xmin,xmax]=xmin_max(Lx,Ly,Lz,nc,nc2);
maskb=(mask-xmin)./(xmax-xmin); % normalization of design variables
xminb=zeros(n,1);
xmaxb=ones(n,1);
xold1=maskb;
xold2=maskb;
low=ones(n,1);
upp=ones(n,1);
a0=1;
a=zeros(m,1);
c_mma=1e3*ones(m,1);
d=zeros(m,1);
%-------------
maxloop=120; C=zeros(1,maxloop); Vol=zeros(1,maxloop);
change=1;
loop =0;    
while change >1e-6 & loop< maxloop
    loop=loop+1;    
    volfrac=volfrac0; 
    alp=min(2.5+loop*(alp0-2.5)/80,alp0); 
    alctr=min(0.979+loop*0.021/100,0.999);
    xold=maskb;
    %% Refinement
    if loop==200
    sc=1;
    nelx=Lx*sc;     nely=Ly*sc;     nelz=Lz*sc;
    dx=Lx/nelx; dy=Ly/nely; dz=Lz/nelz;
    nel = nelx*nely*nelz;                % number of elements
    elsize = [nely nelx nelz];
    [U,F,edofMat,freedofs,Xe,Ye,Ze,iK,jK,posD,iD]=pre_FEM(nelx,nely,nelz,dx,dy,dz,nel);
    [ke]=kemefe(E,nu,sc);
    end
    %% FEM and sensitivities
    [x,sens]=density(nel,nelx,nely,nelz,Xe,Ye,Ze,mask,nc,nc2,Phi,np,alp);
    [c,dc]=FEM(ke,iK,jK,x,nel,elsize,edofMat,freedofs,penal,U,F,posD,iD);
    C(loop)=c;   % objective
    Vol(loop)=mean(x(:));   % volume fraction
    DC=dc'*sens;
    DV=sum(sens)/(nel*volfrac);
    %
    if max(abs(DC))>1e6
        DC=10*DC/max(abs(DC));
        DCDC=0
    end
    DC=DC'.*(xmax-xmin); % normalization of design vatiables
    DV=DV'.*(xmax-xmin);
    %% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    xval=maskb;
    f0val=c;
    df0dx=DC;
    df0dx2=0*df0dx;
    fval=sum(x(:))/(nel*volfrac)-1;
    dfdx=DV;
    dfdx2=0*dfdx;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,loop,xval,...
        xminb,xmaxb,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c_mma,d,alctr);
    maskb=xmma;
    mask=maskb.*(xmax-xmin)+xmin;
    
    %%
    xold2=xold1;
    xold1=maskb;
    change = max(max(abs(maskb-xold)));
    figure(1)
    plotFaces(nel,x,elsize,volfrac0,0)
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',full(c)) ...
        ' Vol.: ' sprintf('%6.3f',Vol(loop)) ...
        ' ch.: ' sprintf('%6.3f',change) ])   
end
%% plot
save mask mask
phi2stl(x-0.5, 10, 1, 1, 50)
%
figure(3)
Rmat=reshape(mask(1+3*nc:nc2+3*nc),np,nc/2);
hold on
for wj=1:nc/2
    XP1=Rmat(:,wj).*cos(Phi(1:end-1));
    YP1=Rmat(:,wj).*sin(Phi(1:end-1));
    XP=[XP1;XP1(1)];
    YP=[YP1;YP1(1)];
    plot(XP,YP,'linewidth',1)
end
axis equal
