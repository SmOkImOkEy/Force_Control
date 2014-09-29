function [Ngal_CP,Dgal_CP]=coprime_Factorization(Ngal,Dgal)

    Sm=double(create_Sm(Ngal,Dgal)); % create coefficient matrix
    subMat_len=size(Dgal,1);
    colLen=2*subMat_len;

    indepInd=indepCols(Sm,subMat_len);
    zEigLoc=1:size(Sm,2);
    zEigLoc(indepInd)=0;
    zEigLoc=reshape(zEigLoc,subMat_len,[]);
    
    % find first nonzero element of each row:
      [row, col]=find(zEigLoc);
      FzeroLoc = zEigLoc((1-subMat_len:0)'+accumarray(row,col,[],@min)*subMat_len);
      z=zeros(size(Sm,1),subMat_len)
    for n=1:subMat_len
        Smt=Sm(:,1:FzeroLoc(n));
        temp=null(Smt);
        weight=sum((Smt*temp).^2,1);
        [v,wI]=min(weight);
        disp(v)
        z(1:size(temp,1),n)=temp(:,wI)/temp(end,wI);
    end
    z(z<size(z,1)*eps(norm(z)))=0;
        
        
    