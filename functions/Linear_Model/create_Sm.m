function Sm=create_Sm(nn,dn,matSize)
% Sm=create_Sm(nn,dn)
% input: nn,dn coefficients matrixes of s. [size(G) ;; rank]
%        of nominator and denominator respectively
% output: Sm matrix - [D0 N0  0  0 0 ... 0 ...           ;
%                      D1 N1 D0 N0 0 ... 0 ...           ;
%                       ...
%                      Dn Nn D(n-1) N(n-1) ... D0 N0 ....;
%                      0   0 Dn     Nn     ...             ]

    nnT=permute(nn,[1 3 2]); % reordering to fit into Sm matrix
    nnS=(reshape(nnT,[],size(nn,1)));
    dnT=permute(dn,[1 3 2]);
    dnS=(reshape(dnT,[],size(nn,1)));

    powerCol_len=size(nnS,1); % length of each block
    subMat_len=size(nnS,2); % length of each sub-matrix
    colLen=2*subMat_len; % block column length
    subMat_amount=size(nn,3); % number of submatrixes
    if ~exist('matSize','var')
        matSize=colLen*(subMat_amount-1); % default matSize value
    end
    Sm=sym(zeros(matSize));
    Sm(1:powerCol_len,1:colLen)=[dnS,nnS];

    for n=(colLen+1):colLen:size(Sm,2)
        Sm(:,n:n+colLen-1)=circshift(Sm(:,n-colLen:n-1),[subMat_len 0]);
    end
end
