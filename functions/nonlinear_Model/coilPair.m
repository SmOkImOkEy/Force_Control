function [B, gradBx, gradBy]=coilPair(l,rho,R,Ihelm,Imaxw,axis,L)
% coilPair(l,rho,R,Ihelm,Imaxw,axis,L)   
% l: on axis samples to calculate field
% rho: rho samples
% Ihelm: Helmholtz current (1 or 0)
% Imaxw: Maxwell Current (1 or 0)
% axis: 'x' or 'y', axis connects two coils
% L: distance between coils
        [B_rho1, BrR1, BrZ1]=CalculateBr(rho,l+L/2,R,Ihelm+Imaxw);
        [B_rho2, BrR2, BrZ2]=CalculateBr(rho,l-L/2,R,Ihelm-Imaxw);
        Br=B_rho1+B_rho2;
        BrR=BrR1+BrR2;
        BrZ=BrZ1+BrZ2;
        
        [B_z1, BzR1, BzZ1]=CalculateBz(rho,l+L/2,R,Ihelm+Imaxw);
        [B_z2, BzR2, BzZ2]=CalculateBz(rho,l-L/2,R,Ihelm-Imaxw);
        Bz=B_z1+B_z2;
        BzR=BzR1+BzR2;
        BzZ=BzZ1+BzZ2;
               
    % create total fields
    dupZ=@(Az) [Az(end:-1:2,:,:); Az(1:end,:,:)];
    dupR=@(Ar) [-Ar(end:-1:2,:,:); Ar(1:end,:,:)];
    % Converting to cartesian coordinates
    if axis=='x'
        B.x=dupZ(Bz);
        B.y=dupR(Br);

       % gradients:
        gradBx.x=dupZ(BzZ);
        gradBx.y=dupR(BzR);

        gradBy.x=dupR(BrZ);
        gradBy.y=dupZ(BrR);
    
    elseif axis=='y' % if coils on y axis, transpose
        B.y=permute(dupZ(Bz),[2 1 3]);
        B.x=permute(dupR(Br),[2 1 3]);

       % gradients:
        gradBy.y=permute(dupZ(BzZ),[2 1 3]);
        gradBy.x=permute(dupR(BzR),[2 1 3]);

        gradBx.y=permute(dupR(BrZ),[2 1 3]);
        gradBx.x=permute(dupZ(BrR),[2 1 3]);  
    end
end

