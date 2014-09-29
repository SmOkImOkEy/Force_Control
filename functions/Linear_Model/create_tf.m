function [Ntf, Dtf]=create_tf(Ngal,Dgal)
% [Ntf, Dtf]=create_tf(Ngal,Dgal)
% turn Ngal,Dgal to transfer functions:
%            N(s)=Ngal(:,:,1)+Ngal(:,:,2)s+...+Ngal(:,:,p)s^p
%            D(s)=Dgal(:,:,1)+Dgal(:,:,2)s+...+Dgal(:,:,p)s^p
%
    subMat_len=size(Dgal,1);

    N_cp=cell(subMat_len);
    D_cp=cell(subMat_len);
    for n=1:subMat_len
        for m=1:subMat_len
            N_cp{n,m}=permute(Ngal(n,m,end:-1:1),[1 3 2]);
            D_cp{n,m}=permute(Dgal(n,m,end:-1:1),[1 3 2]);
        end
    end
     Ntf=tf(N_cp,1); 
     Dtf=tf(D_cp,1); 
end