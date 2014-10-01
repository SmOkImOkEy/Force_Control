function Sm=create_Sm_rows(Ngal,Dgal,colNum)
% Sm=create_Sm(Ngal,Dgal)
% input: Ngal,Dgal coefficients matrixes of s. [size(G) ;; rank]
%        of nominator and denominator respectively
% output: Sm matrix - [D0 N0  0  0 0 ... 0 ...           ;
%                      D1 N1 D0 N0 0 ... 0 ...           ;
%                       ...
%                      Dn Nn D(n-1) N(n-1) ... D0 N0 ....;
%                      0   0 Dn     Nn     ...             ]

    NgalS=reshape(Ngal,size(Dgal,1),[]);
    DgalS=reshape(Dgal,size(Dgal,1),[]);

    powerRow_len=size(DgalS,2); % length of each block
    subMat_len=size(DgalS,1); % length of each sub-matrix
    blockWidth=2*subMat_len; % block row length
    subMat_amount=size(Ngal,3); % number of submatrixes
    
    Sm=zeros(2*(colNum-powerRow_len)+blockWidth,colNum);
    Sm(1:blockWidth,1:powerRow_len)=[DgalS;NgalS];

    for n=(blockWidth+1):blockWidth:size(Sm,1)
        Sm(n:n+blockWidth-1,:)=circshift(Sm(n-blockWidth:n-1,:),[0 subMat_len]);
    end
end
