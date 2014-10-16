% function indepInd=indepCols(Sm,tol)
tol=1e-10
% Sm is square
    dR=diag(qr(Sm));
%     if rank(Sm(:,1:9))<9
%         error('first 9 columns dependent')
%     end
    rem=[];
    n=1;
    indexes=10:numel(dR);
    while n<=numel(indexes) && ~isempty(indexes)
        if abs(dR(indexes(n))/dR(indexes(n)-6))<tol
            rem=cat(2,rem,-9+(indexes(n):6:numel(dR)));
        else
            n=n+1;
        end
        indexes=10:numel(dR);
        indexes(rem)=[];
    end
    indepInd=cat(2,1:9,indexes);
% end