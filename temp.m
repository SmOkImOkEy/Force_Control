 dep=[];
    if ~exist('tol','var')
        tol=1e-10;
    end
    rB=abs(diag(qr(Sm)));
    for n=1:numel(rB)
        r=rB(1:n);
        rn=r(n)
        r(dep)=[];
        rat=rn/max(r);
        if rat<tol
            dep=cat(2,dep,n);
        end
    end