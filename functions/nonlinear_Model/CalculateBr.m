function [Br BrR BlL]=CalculateBr(rho,l,R,I,phiAcc)
global mu0
if ~exist('phiAcc','var')
    phiAcc=50;
end
if ~exist('I','var')
    I=1;
end

phi=linspace(0,2*pi,phiAcc);
[L Rho Phi]=meshgrid(l,rho,phi);

Z=(R^2+Rho.^2+L.^2-2*R*Rho.*cos(Phi));
Im=repmat(I,[size(L(:,:,1)) 1]);
% magnetic Field on Rho Direction
    trapB=trapz(phi,R*L.*cos(Phi)./Z.^(3/2),3);
    Br=mu0/(2*pi)*Im.*repmat(trapB,[1 1 numel(I)]);

% derivative by Rho
    trapBr=trapz(phi,-3*R*L.*cos(Phi).*(Rho-R*cos(Phi))./Z.^2.5,3);
    BrR=mu0/(2*pi)*Im.*repmat(trapBr,[1 1 numel(I)]);
% derivative by l
    trapBl=trapz(phi,(R*cos(Phi).*Z-3*R*L.^2.*cos(Phi))./Z.^2.5,3);
    BlL=mu0/(2*pi)*Im.*repmat(trapBl,[1 1 numel(I)]);
end