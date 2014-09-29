function [Ngal_CP,Dgal_CP]=coprime_Factorization(LNgal,LDgal,tol)

    Sm=double(create_Sm(LNgal,LDgal)); % create coefficient matrix
    subMat_len=size(LDgal,1);
    colLen=2*subMat_len;

    modu=@(a,b) mod(a,b)+b*(mod(a,b)==0);
    dR=diag(qr(Sm)); % diagonal of qr factorization
    if ~exist('tol','var')
        tol=1e-10;%sqrt(max(size(Sm))*eps(norm(Sm)));
    end
    diagZerosLoc=find(abs(dR)<tol); % locate zeros on diagonal
    
%     diagZerosLoc=rowDependency(Sm.');
    % unique unsorted: find all main dependences and only them.
      diagZeroLocGroups=modu(diagZerosLoc,subMat_len);
      [~, i1, ~]=unique(diagZeroLocGroups); % groupIndex - which group belong value
      FZeroloc=diagZerosLoc(sort(i1));
    % --
    % reset null matrix
        z=zeros(size(Sm,1),subMat_len);
    % --
    remove=[]; % columns to remove from Sm
    indexes=1:size(Sm,2); %  default index - no remove
    for zeroIndex=1:numel(FZeroloc) % calculate for each zero eigen value.
        currIndex=indexes(indexes<=FZeroloc(zeroIndex)); % pick items from left
                  % to the 'zeroIndex's' main dependence, exclude removed items
        Smt=Sm(:,currIndex); % pick submatrix
        temp=null(Smt); % solve homogenus equation.
      % error handler:
        if size(temp,2)>1 
            disp('nullity larger than 1 - check tolerance')
        elseif isempty(temp)
            disp('no nullity')
            break;
        else
            temp=temp/temp(end); % normalizing last item - monic null vector.
            group=modu(FZeroloc(zeroIndex),subMat_len);
            if group>0
                if all(z(:,group)==0)
                    z(currIndex(1:size(temp,1)),group)=temp;
                else
                    disp('dependency on denominator - check properness')
                    return 
                end % if all(z(:,group)==0)
            end % if group>0
            remove=cat(2,remove,FZeroloc(zeroIndex):colLen:size(Sm,2));
            indexes=1:size(Sm,2);
            indexes(remove)=[];
        end % if size(temp,2)>1
      % ---
    end % for
    % creating co-prime matrixes using null results
    mat=permute(reshape(z.',subMat_len,subMat_len,[]),[2 1 3]);
    Ngal_CP=-mat(:,:,1:2:end);
    Dgal_CP=mat(:,:,2:2:end);
   
     degD=max(max(sum(Dgal_CP~=0,3))); % degree of D
     Ngal_CP=Ngal_CP(:,:,1:degD);
     Dgal_CP=Dgal_CP(:,:,1:degD);
end