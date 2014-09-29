function [Ngal_CP,Dgal_CP]=coprime_Factorization(Ngal,Dgal)

    Sm=double(create_Sm(Ngal,Dgal)); % create coefficient matrix
    subMat_len=size(Dgal,1);
    colLen=2*subMat_len;

    modu=@(a,b) mod(a,b)+b*(mod(a,b)==0);

    indepInd=indepCols(Sm,subMat_len);
    zEigLoc=1:size(Sm,2)
    zEigLoc(indepInd)=[];
    zEigLoc=reshape(zEigLoc,subMat_len,[])

    group=modu(zEigLoc,subMat_len)
    
    
    
    