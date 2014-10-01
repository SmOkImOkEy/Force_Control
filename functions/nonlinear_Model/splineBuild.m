function [an, bn]=splineBuild(Xt,Yt,tt)
% input: Xt Yt route points
%    *time dependency is on column 
% output: cubic spline coefficients
tt=tt';
% reset spline coefficients
    an=zeros(3,numel(tt));
    bn=zeros(3,numel(tt)); % spline coefficients: (cubic spline)
    
    for n=2:numel(tt)-1
        T0=[tt(n-1)^2   tt(n-1) 1;
            tt(n)^2     tt(n)   1;
            tt(n+1)^2   tt(n+1) 1];
        
        an(:,n)=T0\Xt(n-1:n+1)';
        bn(:,n)=T0\Yt(n-1:n+1)';
    end
    an(:,1)=an(:,2);
    bn(:,1)=bn(:,2);
    an(:,numel(tt))=an(:,numel(tt)-1);
    bn(:,numel(tt))=bn(:,numel(tt)-1);
end