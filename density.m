function [x,sens]=density(nel,nelx,nely,nelz,Xe,Ye,Ze,mask,nc,nc2,Phi,np,alp)
%% Geometric parameters
XC=mask(1:nc);
YC=mask(1+nc:2*nc);
ZC=mask(1+2*nc:3*nc);
Rmat=reshape(mask(1+3*nc:nc2+3*nc),np,nc/2);
Ra=mask(1+3*nc+nc2:end);
%%
delta=3;
n1=zeros(nc/2,3); n1(:,1)=1;    % normal unit vector of global coordinate system
n2=zeros(nc/2,3); n2(:,2)=1;
% n3=zeros(nc/2,3); n3(:,3)=1;
Nz=[XC(2:2:end)-XC(1:2:end), YC(2:2:end)-YC(1:2:end), ZC(2:2:end)-ZC(1:2:end)];
Nz_dis=sqrt(sum(Nz.^2,2));
nz=Nz./Nz_dis;
coin_nxnz=find(n1(:,1)==nz(:,1)&n1(:,2)==nz(:,2)&n1(:,3)==nz(:,3));
Ny=cross(n1,nz);
Ny(coin_nxnz,:)=cross(n2(coin_nxnz,:),nz(coin_nxnz,:));
Nx=cross(nz,Ny);
Ny_dis=sqrt(sum(Ny.^2,2)); ny=Ny./Ny_dis;   % normal unit vector
Nx_dis=sqrt(sum(Nx.^2,2)); nx=Nx./Nx_dis;   % normal unit vector
%%
XPmat1=Rmat.*cos(Phi(1:end-1));
YPmat1=Rmat.*sin(Phi(1:end-1));
XPmat=[XPmat1;XPmat1(1,:)];
YPmat=[YPmat1;YPmat1(1,:)];
%%
T=zeros(np,2);
sens=zeros(nel,3*nc+nc2+nc/2); % sensitiviy
%%
pai=ones(nel,nc/2);   % element density matrix
for wj=1:nc/2
    k1P=[Xe-XC(2*wj-1), Ye-YC(2*wj-1), Ze-ZC(2*wj-1)];
    xe=k1P*nx(wj,:)';
    ye=k1P*ny(wj,:)';
    ze=k1P*nz(wj,:)';
    %% Calculate density
    XP=XPmat(:,wj);
    YP=YPmat(:,wj);
    inb=inpolygon(xe,ye,XP,YP);
    inpoly=find(inb~=0&ze>=-Ra(wj)&ze<=Nz_dis(wj)+Ra(wj)); % index of polygon area
    XPO=ones(nel,1); % element density for each component void
    %% vector xi_xi+1
    T(:,1)=XP(2:end)-XP(1:end-1);
    T(:,2)=YP(2:end)-YP(1:end-1);
    T_dis=sqrt(sum(T.^2,2));
    t=T./T_dis;
    %% sensitivity of local noraml unit vector with respec to xk1, xk2
    nz1=nz(wj,1); nz2=nz(wj,2); nz3=nz(wj,3);
    invertNz_dis=1/Nz_dis(wj);
    dnzdxk1=[-1+nz1^2, nz1*nz2, nz1*nz3]*invertNz_dis;
    dnzdyk1=[nz1*nz2, -1+nz2^2, nz2*nz3]*invertNz_dis;
    dnzdzk1=[nz1*nz3, nz2*nz3, -1+nz3^2]*invertNz_dis;
    if n1(wj,1)==nz(wj,1)&n1(wj,2)==nz(wj,2)&n1(wj,3)==nz(:,3)
        dNydxk1=cross(n2(wj,:),dnzdxk1);
        dNydyk1=cross(n2(wj,:),dnzdyk1);
        dNydzk1=cross(n2(wj,:),dnzdzk1);
        dnydxk1=(dNydxk1-ny(wj,:)*(ny(wj,:)*dNydxk1'))/Ny_dis(wj);
        dnydyk1=(dNydyk1-ny(wj,:)*(ny(wj,:)*dNydyk1'))/Ny_dis(wj);
        dnydzk1=(dNydzk1-ny(wj,:)*(ny(wj,:)*dNydzk1'))/Ny_dis(wj);
        dnydxk2=-dnydxk1;
        dnydyk2=-dnydyk1;
        dnydzk2=-dnydzk1;
        dNxdxk1=cross(dnzdxk1,Ny(wj,:))+cross(nz(wj,:),dNydxk1);
        dNxdyk1=cross(dnzdyk1,Ny(wj,:))+cross(nz(wj,:),dNydyk1);
        dNxdzk1=cross(dnzdzk1,Ny(wj,:))+cross(nz(wj,:),dNydzk1);
        dnxdxk1=(dNxdxk1-nx(wj,:)*(nx(wj,:)*dNxdxk1'))/Nx_dis(wj);
        dnxdyk1=(dNxdyk1-nx(wj,:)*(nx(wj,:)*dNxdyk1'))/Nx_dis(wj);
        dnxdzk1=(dNxdzk1-nx(wj,:)*(nx(wj,:)*dNxdzk1'))/Nx_dis(wj);
        dnxdxk2=-dnxdxk1;
        dnxdyk2=-dnxdyk1;
        dnxdzk2=-dnxdzk1;
    else
        dNydxk1=cross(n1(wj,:),dnzdxk1);
        dNydyk1=cross(n1(wj,:),dnzdyk1);
        dNydzk1=cross(n1(wj,:),dnzdzk1);
        dnydxk1=(dNydxk1-ny(wj,:)*(ny(wj,:)*dNydxk1'))/Ny_dis(wj);
        dnydyk1=(dNydyk1-ny(wj,:)*(ny(wj,:)*dNydyk1'))/Ny_dis(wj);
        dnydzk1=(dNydzk1-ny(wj,:)*(ny(wj,:)*dNydzk1'))/Ny_dis(wj);
        dnydxk2=-dnydxk1;
        dnydyk2=-dnydyk1;
        dnydzk2=-dnydzk1;
        dNxdxk1=cross(dnzdxk1,Ny(wj,:))+cross(nz(wj,:),dNydxk1);
        dNxdyk1=cross(dnzdyk1,Ny(wj,:))+cross(nz(wj,:),dNydyk1);
        dNxdzk1=cross(dnzdzk1,Ny(wj,:))+cross(nz(wj,:),dNydzk1);
        dnxdxk1=(dNxdxk1-nx(wj,:)*(nx(wj,:)*dNxdxk1'))/Nx_dis(wj);
        dnxdyk1=(dNxdyk1-nx(wj,:)*(nx(wj,:)*dNxdyk1'))/Nx_dis(wj);
        dnxdzk1=(dNxdzk1-nx(wj,:)*(nx(wj,:)*dNxdzk1'))/Nx_dis(wj);
        dnxdxk2=-dnxdxk1;
        dnxdyk2=-dnxdyk1;
        dnxdzk2=-dnxdzk1;
    end
    for wi=1:length(inpoly)
        e=inpoly(wi);
        a=[xe(e)-XP(1:end-1),ye(e)-YP(1:end-1)];  % vector xe-Pi
        b=a-sum(a.*t,2).*t;
        at=sum(a.*t,2);
        idmid=find(at>0&at<T_dis);
        idright=find(at<=0);
        b_dis2=sum(b(idmid,:).*b(idmid,:),2);
        a_dis2=sum(a(idright,:).*a(idright,:),2);
        b_dis2min=min(b_dis2);
        a_dis2min=min(a_dis2);
        %%
        if (b_dis2min<=a_dis2min & length(idmid)>0&length(idright)>0) | length(idright)<1
            de=sqrt(b_dis2min);
            index=idmid(find(b_dis2==b_dis2min,1));
            if ze(e)>0 & ze(e)<Nz_dis(wj)
                [XPO,sens]=dens_mid(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2);
            elseif ze(e)<=0
                Re=sqrt((Xe(e)-XC(2*wj-1))^2+(Ye(e)-YC(2*wj-1))^2+(Ze(e)-ZC(2*wj-1))^2);
                dde=Ra(wj)-Re;
                if Re<=Ra(wj)&dde<=de
                    de=dde;
                    XPOe=1/(1+exp(-alp*(-de+delta)));
                    XPO(e)=XPOe;
                    dedxk1=(Xe(e)-XC(2*wj-1))/Re;
                    dedyk1=(Ye(e)-YC(2*wj-1))/Re;
                    dedzk1=(Ze(e)-ZC(2*wj-1))/Re;
                    dphide=-alp*exp(-alp*(-de+delta))*XPOe;
                    sens(e,2*wj-1)=sens(e,2*wj-1)+dphide*dedxk1;
                    sens(e,2*wj-1+nc)=sens(e,2*wj-1+nc)+dphide*dedyk1;
                    sens(e,2*wj-1+2*nc)=sens(e,2*wj-1+2*nc)+dphide*dedzk1;
                    sens(e,wj+3*nc+nc2)=sens(e,wj+3*nc+nc2)+dphide;
                elseif Re<=Ra(wj)&dde>de
                    [XPO,sens]=dens_mid(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2);
                else
                end
            elseif ze(e)>=Nz_dis(wj)
                Re=sqrt((Xe(e)-XC(2*wj))^2+(Ye(e)-YC(2*wj))^2+(Ze(e)-ZC(2*wj))^2);
                dde=Ra(wj)-Re;
                if Re<=Ra(wj)&dde<=de
                    de=dde;
                    XPOe=1/(1+exp(-alp*(-de+delta)));
                    XPO(e)=XPOe;
                    dedxk2=(Xe(e)-XC(2*wj))/Re;
                    dedyk2=(Ye(e)-YC(2*wj))/Re;
                    dedzk2=(Ze(e)-ZC(2*wj))/Re;
                    dphide=-alp*exp(-alp*(-de+delta))*XPOe;
                    sens(e,2*wj)=sens(e,2*wj)+dphide*dedxk2;
                    sens(e,2*wj+nc)=sens(e,2*wj+nc)+dphide*dedyk2;
                    sens(e,2*wj+2*nc)=sens(e,2*wj+2*nc)+dphide*dedzk2;
                    sens(e,wj+3*nc+nc2)=sens(e,wj+3*nc+nc2)+dphide;
                elseif Re<=Ra(wj)&dde>de
                    [XPO,sens]=dens_mid(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2);
                else                    
                end
            else
            end            
        elseif (b_dis2min>a_dis2min & length(idmid)>0&length(idright)>0) | (length(idmid)<1 & length(idright)>0)
            de=sqrt(a_dis2min);
            index=idright(find(a_dis2==a_dis2min,1));    % index of segment
            if ze(e)>0 & ze(e)<Nz_dis(wj)                
                [XPO,sens]=dens_mid_nonconvex(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2);
            elseif ze(e)<=0
                Re=sqrt((Xe(e)-XC(2*wj-1))^2+(Ye(e)-YC(2*wj-1))^2+(Ze(e)-ZC(2*wj-1))^2);
                dde=Ra(wj)-Re;
                if Re<=Ra(wj)&dde<=de
                    de=dde;
                    XPOe=1/(1+exp(-alp*(-de+delta)));
                    XPO(e)=XPOe;
                    dedxk1=(Xe(e)-XC(2*wj-1))/Re;
                    dedyk1=(Ye(e)-YC(2*wj-1))/Re;
                    dedzk1=(Ze(e)-ZC(2*wj-1))/Re;
                    dphide=-alp*exp(-alp*(-de+delta))*XPOe;
                    sens(e,2*wj-1)=sens(e,2*wj-1)+dphide*dedxk1;
                    sens(e,2*wj-1+nc)=sens(e,2*wj-1+nc)+dphide*dedyk1;
                    sens(e,2*wj-1+2*nc)=sens(e,2*wj-1+2*nc)+dphide*dedzk1;
                    sens(e,wj+3*nc+nc2)=sens(e,wj+3*nc+nc2)+dphide;
                elseif Re<=Ra(wj)&dde>de
                    [XPO,sens]=dens_mid_nonconvex(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2);
                else
                end
            elseif ze(e)>=Nz_dis(wj)
                Re=sqrt((Xe(e)-XC(2*wj))^2+(Ye(e)-YC(2*wj))^2+(Ze(e)-ZC(2*wj))^2);
                dde=Ra(wj)-Re;
                if Re<=Ra(wj)&dde<=de
                    de=dde;
                    XPOe=1/(1+exp(-alp*(-de+delta)));
                    XPO(e)=XPOe;
                    dedxk2=(Xe(e)-XC(2*wj))/Re;
                    dedyk2=(Ye(e)-YC(2*wj))/Re;
                    dedzk2=(Ze(e)-ZC(2*wj))/Re;
                    dphide=-alp*exp(-alp*(-de+delta))*XPOe;
                    sens(e,2*wj)=sens(e,2*wj)+dphide*dedxk2;
                    sens(e,2*wj+nc)=sens(e,2*wj+nc)+dphide*dedyk2;
                    sens(e,2*wj+2*nc)=sens(e,2*wj+2*nc)+dphide*dedzk2;
                    sens(e,wj+3*nc+nc2)=sens(e,wj+3*nc+nc2)+dphide;
                elseif Re<=Ra(wj)&dde>de
                    [XPO,sens]=dens_mid_nonconvex(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2);
                else
                end
            else
            end
        else
        end
    end
    pai(:,wj)=XPO;            
end
X=1-prod(pai,2);
X(find(X<1e-4))=0;
X(find(X>0.9999))=1;
x=reshape(X,nely,nelx,nelz);
sens=-(1-X).*sens;


