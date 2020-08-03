function [xmin,xmax]=xmin_max(Lx,Ly,Lz,nc,nc2)
Xmin=0*ones(nc,1);
Xmax=Lx*ones(nc,1);
Ymin=0*ones(nc,1);
Ymax=Ly*ones(nc,1); %%%%%%%%%%%%
Zmin=0*ones(nc,1);
Zmax=Lz*ones(nc,1);
%
Rdmin=4.5*ones(nc2,1);
Rdmax=15*ones(nc2,1);
Ramin=3*ones(nc/2,1);
Ramax=15*ones(nc/2,1);
xmin=[Xmin; Ymin; Zmin; Rdmin;Ramin];
xmax=[Xmax; Ymax; Zmax; Rdmax;Ramax];

