function [Ngal, Dgal]=Coefficient_Matrixes(LN,LD)  
% [Ngal, Dgal]=Coefficient_Matrixes(LN,LD)  
% 
% Calculate Coefficient matrix
% N=Ngal(:,:,1)+Ngal(:,:,2)s+...+Ngal(:,:,p)s^p
% D=Dgal(:,:,1)+Dgal(:,:,2)s+...+Dgal(:,:,p)s^p
%
% (Dgal(:,:,p)~=0 (stricly proper))

    [num,~]=tfdata(LN);
    [den,~]=tfdata(LD);
    subMat_len=size(LD,1);
   
   cellsz = cellfun(@numel,den,'uni',false);
   degree=max(cell2mat(cellsz(:))); % size of 3rd dimention of ngal,dgal
   Ngal=zeros([subMat_len subMat_len degree]);
   Dgal=zeros([subMat_len subMat_len degree]);
   
    for n=1:subMat_len
        for m=1:subMat_len
            Ngal(n,m,1:numel(num{n,m}))=permute(num{n,m}(end:-1:1),[1,3,2]); 
            Dgal(n,m,1:numel(den{n,m}))=permute(den{n,m}(end:-1:1),[1,3,2]);    
        end
    end
end