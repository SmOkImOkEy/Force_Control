function [Bi, dBxi, dByi]=B_Field_Interp(x,y,I)    
% [Bi, dBxi, dByi]=XY_interp(X,Y,x,y,B,dBx,dBy)   
% interpolate Magnetic Field and gradient on x,y point
% I: [Ix_helmholtz; Ix_maxwell; Iy_helmholtz; Iy_maxwell]
% the current flows through each coil
global X Y B dBx dBy
    B_calculated=reshape(B,[],4,2);
    dBx_calculated=reshape(dBx,[],4,2);
    dBy_calculated=reshape(dBy,[],4,2);
%  
    [xI, yI]=findCube2D(x,y);
    ind = sub2ind(size(X),yI,xI) ;
    disp(ind)
    dst=sqrt((x-X(ind)).^2+(y-Y(ind)).^2);
%    
    interpLin=@(Fc) (1./dst(:)')*Fc(:)/sum(1./dst(:));
    if size(ind)==1
        Bi=squeeze(B_calculated(ind,:,:)).'*I;
        dBxi=squeeze(dBx_calculated(ind,:,:)).'*I;
        dByi=squeeze(dBy_calculated(ind,:,:)).'*I;
        disp('on it')
    else
        BInterMat=zeros([size(B,3) size(B,4)]);
        dBxInterMat=BInterMat;dByInterMat=BInterMat;
        disp('interpolating it')
        for n=1:size(B,3)
            for m=1:size(B,4)
                BInterMat(n,m)=squeeze(interpLin(B_calculated(ind,n,m)));
                dBxInterMat(n,m)=squeeze(interpLin(dBx_calculated(ind,n,m)));
                dByInterMat(n,m)=squeeze(interpLin(dBy_calculated(ind,n,m)));
            end
        end
        Bi=BInterMat.'*I;
        dBxi=dBxInterMat.'*I;
        dByi=dByInterMat.'*I;
    end
    
    
end