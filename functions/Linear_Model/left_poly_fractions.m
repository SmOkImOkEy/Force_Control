function [N, D]=left_poly_fractions(Gtf)
% [N, D]=leftFractions(Gtf)
% calculating N and D so that:  Gtf=D^-1*N
% N, D are not coprime, just left factors

    subMat_len=size(Gtf,1);
    [~,p,~]=zpkdata(Gtf); 
    pLCD=cell(subMat_len);
    for n=1:subMat_len
        prow=p{1,n};
        for m=1:subMat_len
            ila=ismember(p{n,m},prow);
            prow=[prow; p{n,m}(~ila)];
            pLCD{n,m}=0;
        end
        pLCD{n,n}=poly(prow);
    end
    D=minreal(tf(pLCD,1));
    N=minreal(D*Gtf);
end