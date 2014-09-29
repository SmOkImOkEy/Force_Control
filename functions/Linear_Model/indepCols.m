function indepInd=indepCols(Sm,subMat_len)
    dR=diag(abs(qr(Sm)));
    indexes=reshape(1:numel(dR),subMat_len,[]);
    inD=indexes(:,1:2:end);
    inN=indexes(:,2:2:end);

    [~,sortIndex] = sort(dR(inN(:)),'descend');  %# Sort the values in
    topNInd = inN(sortIndex(1:rank(Sm)-numel(inD)));  %# Get a linear index into A of the 5 largest values
    indepInd=sort([inD(:); topNInd(:)]);
end