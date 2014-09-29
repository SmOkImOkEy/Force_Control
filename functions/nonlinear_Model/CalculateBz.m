function [Bl, BlR, BlL]=CalculateBz(rho,l,R,I,phiAcc)
global mu0
if ~exist('phiAcc','var')
    phiAcc=50;
end
if ~exist('I','var')
    I=1;
end

phi=linspace(0,2*pi,phiAcc);
[L, Rho, Phi]=meshgrid(l,rho,phi);

Z=(R^2+Rho.^2+L.^2-2*R.*Rho.*cos(Phi));
Im=repmat(I,[size(L(:,:,1)) 1]);

% magnetic Field on l Direction
        trapzBl=trapz(phi,(R^2-R*Rho.*cos(Phi))./Z.^(3/2),3);
        Bl=mu0/(2*pi)*Im.*repmat(trapzBl,[1 1 numel(I)]);

% derivative by Rho
        trapzBlR=trapz(phi,(R*cos(Phi).*Z+3*(R^2-R*Rho.*cos(Phi)).*(Rho-R*cos(Phi)))./Z.^2.5,3);
        BlR=-mu0/(2*pi)*Im.*repmat(trapzBlR,[1 1 numel(I)]);
        
% derivative by l
        trapzBlL=trapz(phi,(-3*L.*(R^2-R*Rho.*cos(Phi)))./Z.^2.5,3);
        BlL=mu0/(2*pi)*Im.*repmat(trapzBlL,[1 1 numel(I)]);
end