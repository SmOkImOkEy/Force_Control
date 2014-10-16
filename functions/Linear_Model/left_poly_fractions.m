% function [N, D]=left_poly_fractions(Gtf)
% [N, D]=leftFractions(Gtf)
% calculating N and D so that:  Gtf=D^-1*N
% N, D are not coprime, just left factors

    subMat_len=size(Gtf,1);
    [~,p,~]=zpkdata(Gtf); 
    cellLen=cellfun('length',p);
    [~, sI]=sort(cellLen,2,'descend');
    pLCD=cell(subMat_len);
   for row=1:subMat_len
        prow=p{row,sI(row,1)};
        for col=1:subMat_len
           ila=ismember(p{row,sI(row,col)},prow);
           prow=[prow; p{row,sI(row,col)}(~ila)];
           pLCD{row,sI(row,col)}=0;
        end
        pLCD{row,row}=poly(prow);
    end
    D=(tf(pLCD,1));
    N=minreal(D*Gtf);
% end