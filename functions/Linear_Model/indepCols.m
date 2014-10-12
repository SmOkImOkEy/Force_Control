function indepInd=indepCols(Sm,subMat_len)
    dR=diag(abs(qr(Sm)));
%     dR=[dR ;zeros(size(Sm,2)-numel(dR),1)];
    indexes=reshape(1:size(Sm,2),subMat_len,[]);
    inD=indexes(:,5:2:end);
    inN=indexes(:,4:2:end);

    [~,sortIndex] = sort(dR(inN(:)),'descend');  %# Sort the values in
    topNInd = inN(sortIndex(1:rank(Sm)-9-numel(inD)));  %# Get a linear index into A of the 5 largest values
    indepInd=sort([(1:9)';inD(:); topNInd(:)]);
end