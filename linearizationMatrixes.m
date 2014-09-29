%linearization 2
clear all;clc
syms Fx_hat Fy_hat phiB_hat % inputvector
syms Bm_hat % input parameters
syms Ihx Ihy Imx Imy theta thetaD % state vector
syms Rx Ry Lx Ly M h kappa alpha I % constants

Ihx=sym('Ihx','real');
Ihy=sym('Ihy','real');
Imx=sym('Imx','real');
Imy=sym('Imy','real');



X=[Ihx;Ihy;Imx;Imy;theta;thetaD];
U=[Fx_hat; Fy_hat; phiB_hat];

Ihx_hat=Bm_hat/h*cos(phiB_hat);
Ihy_hat=Bm_hat/h*sin(phiB_hat);
Imx_hat=2/(3*M*kappa)*(2*Fx_hat/cos(phiB_hat)+Fy_hat/sin(phiB_hat));
Imy_hat=2/(3*M*kappa)*(Fx_hat/cos(phiB_hat)+2*Fy_hat/sin(phiB_hat));



F_mat=[Rx/Lx*(-Ihx+Ihx_hat);
    Ry/Ly*(-Ihy+Ihy_hat);
    Rx/Lx*(-Imx+Imx_hat);
    Ry/Ly*(-Imy+Imy_hat);
    thetaD;
    M*h/I*sqrt(Ihx^2+Ihy^2)*sin(theta)-alpha/I*thetaD];

phiB=atan2(Ihy,Ihx);
phiM=phiB-theta;

G_mat=[M*kappa*cos(phiM)*(Imx-1/2*Imy);
    M*kappa*cos(phiM)*(-1/2*Imx+Imy);
    phiB];
A=sym(zeros(6));
C=sym(zeros(3,6));
for n=1:numel(X)
    A(:,n)=diff(F_mat,X(n));
    C(:,n)=diff(G_mat,X(n));
end
B=sym(zeros(6,3));
D=sym(zeros(3,3));
for n=1:numel(U)
    B(:,n)=diff(F_mat,U(n));
    D(:,n)=diff(G_mat,U(n));
end

% assumption small angles
A=subs(A,{'theta','thetaD'},{0,0});
B=subs(B,{'theta','thetaD'},{0,0});
C=subs(C,{'theta','thetaD'},{0,0});
D=subs(D,{'theta','thetaD'},{0,0});


save('LinMat','A','B','C','D')
