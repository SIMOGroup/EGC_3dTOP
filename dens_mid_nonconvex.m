function [XPO,sens]=dens_mid_nonconvex(nx,ny,de,t,T_dis,k1P,b,a,index,XPO,sens,nc,Phi,np,alp,delta,e,wj,dnxdxk1,dnxdyk1,dnxdzk1,dnxdxk2,dnxdyk2,dnxdzk2,dnydxk1,dnydyk1,dnydzk1,dnydxk2,dnydyk2,dnydzk2)
XPOe=1/(1+exp(-alp*(-de+delta)));
XPO(e)=XPOe;
dxedxk1=dnxdxk1*k1P(e,:)'-nx(wj,1);
dxedyk1=dnxdyk1*k1P(e,:)'-nx(wj,2);
dxedzk1=dnxdzk1*k1P(e,:)'-nx(wj,3);
dxedxk2=dnxdxk2*k1P(e,:)';
dxedyk2=dnxdyk2*k1P(e,:)';
dxedzk2=dnxdzk2*k1P(e,:)';
dyedxk1=dnydxk1*k1P(e,:)'-ny(wj,1);
dyedyk1=dnydyk1*k1P(e,:)'-ny(wj,2);
dyedzk1=dnydzk1*k1P(e,:)'-ny(wj,3);
dyedxk2=dnydxk2*k1P(e,:)';
dyedyk2=dnydyk2*k1P(e,:)';
dyedzk2=dnydzk2*k1P(e,:)';
% sensitivity of de with respect to coordinates Pi on plane x-y
invertde=1/de;
dedXi=-a(index,1)*invertde;
dedYi=-a(index,2)*invertde;
dedRi=dedXi*cos(Phi(index))+dedYi*sin(Phi(index));
% sensitivity of de with respect to coordinates xk1, xk2
dedxk1=(a(index,1)*dxedxk1+a(index,2)*dyedxk1)*invertde;
dedyk1=(a(index,1)*dxedyk1+a(index,2)*dyedyk1)*invertde;
dedzk1=(a(index,1)*dxedzk1+a(index,2)*dyedzk1)*invertde;
dedxk2=(a(index,1)*dxedxk2+a(index,2)*dyedxk2)*invertde;
dedyk2=(a(index,1)*dxedyk2+a(index,2)*dyedyk2)*invertde;
dedzk2=(a(index,1)*dxedzk2+a(index,2)*dyedzk2)*invertde;
%
dphide=-alp*exp(-alp*(-de+delta))*XPOe;
sens(e,index+(wj-1)*np+3*nc)=sens(e,index+(wj-1)*np+3*nc)+dphide*dedRi;
sens(e,2*wj-1)=sens(e,2*wj-1)+dphide*dedxk1;
sens(e,2*wj-1+nc)=sens(e,2*wj-1+nc)+dphide*dedyk1;
sens(e,2*wj-1+2*nc)=sens(e,2*wj-1+2*nc)+dphide*dedzk1;
sens(e,2*wj)=sens(e,2*wj)+dphide*dedxk2;
sens(e,2*wj+nc)=sens(e,2*wj+nc)+dphide*dedyk2;
sens(e,2*wj+2*nc)=sens(e,2*wj+2*nc)+dphide*dedzk2;