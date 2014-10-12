% for n=1:numel(dR)
%     disp(dR(n));
% end
% %%
% GG(1,1)=zpk(-2,[1,-1],1);
% GG(1,2)=zpk(0,[1 -1 -1],1);
% GG(1,3)=zpk([],[2,1],1);
% 
% GG(2,2)=zpk([],3,3);
% GG(2,3)=zpk([],0,1);
% 
% GG(3,1)=zpk([],2,1);
% GG(3,3)=zpk([],-4,1);
% 
% Gtf=(GG);
% [N, D]=left_poly_fractions(Gtf);

%%
    [~,p,~]=zpkdata(Gtf); 
    
    n=1
    cellLen=cellfun('length',p);
    [~, sI]=sort(cellLen,2,'descend');
    pLCD=cell(subMat_len);
    for row=1:subMat_len
        prow=p{row,sI(row,1)};
        for col=2:subMat_len
           ila=ismember(p{row,sI(row,col)},prow);
           prow=[prow; p{row,sI(row,col)}(~ila)];
           pLCD{row,sI(row,col)}=0;
        end
        pLCD{row,row}=poly(prow);
    end
    D=minreal(tf(pLCD,1));
    N=minreal(D*Gtf);
    
zpk(D)