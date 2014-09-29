function dep=rowDependency(Sm)
    % dep=rowDependency(Sm)
    % outputs indexes of dependent rows

    dep=[];
    for n=6:size(Sm,1)
        curr=1:n;
        curr(dep)=[];
        Smt=Sm(curr,:);
        if rank(Smt,sqrt(size(Smt,2)*eps(norm(Smt))))<size(Smt,1)
            dep=[dep n];
        end
    end
end