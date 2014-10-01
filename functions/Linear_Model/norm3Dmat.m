function norma=norm3Dmat(A)
    noVec=zeros(size(A,3),1);
    for n=1:size(A,3)
        noVec(n)=norm(A(:,:,n));
    end
    norma=norm(noVec);
end