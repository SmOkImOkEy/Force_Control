function indepInd=indepCols(Sm,tol)
    dR=diag(qr(Sm));
%     if rank(Sm(:,1:9))<9
%         error('first 9 columns dependent')
%     end
    rem=[];
    n=1;
   Nindexes=10:numel(dR);
    while n<=numel(Nindexes) && ~isempty(Nindexes)
        if abs(dR(Nindexes(n))/dR(Nindexes(n)-6))<tol
            rem=cat(2,rem,-9+(Nindexes(n):6:numel(dR)));
        else
            n=n+1;
        end
        Nindexes=10:numel(dR);
        Nindexes(rem)=[];
    end
    indepInd=cat(2,1:9,Nindexes);
end