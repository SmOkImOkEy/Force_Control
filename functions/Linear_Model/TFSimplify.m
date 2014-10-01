function [Ngal_s, Dgal_s]=TFSimplify(Ngal,Dgal,tol)
% [Ngal_s, Dgal_s]=TFSimplify(Ngal,Dgal,tol)
%
% Zero coefficients smaller than tol
% shorten simplified array 

    Ngal(abs(Ngal)<tol)=0; % zeroing small elements
    Dgal(abs(Dgal)<tol)=0; % zeroing small elements
    
    degD=max(max(sum(Dgal~=0,3))); % degree of D
    Ngal_s=Ngal(:,:,1:degD); % shorening arrays
    Dgal_s=Dgal(:,:,1:degD);
end