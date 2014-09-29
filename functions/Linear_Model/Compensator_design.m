function [Agal, Bgal, Ftf]=Compensator_design(Ngal,Dgal,poles)
% [Agal, Bgal]=Compensator_design(Ngal,Dgal,poles)
% creates compensator to set poles of the system G=ND^-1
% arbitrarily, 
% number of poles must be larger than number of original poles
    subMat_len=size(Dgal,1);
    m=size(Dgal,3); % highest power in compensator denominator
    mu_d=diag(sum(Dgal~=0,3)-1); % degree of each diagonal element
    Fdeg=mu_d+m; % minimum degree of F (number of poles on each diagonal element)
% create F polynom matrix
    Fpol=cell(subMat_len);
    Fpol(:) = {0};
    Fgal=zeros([subMat_len,subMat_len,max(Fdeg)]);
    for n=1:subMat_len
       pol=real(poly(poles(1:Fdeg(n))));
       Fpol{n,n}=pol(end:-1:1); % desired denominator shaped for TF
       Fgal(n,n,1:numel(pol))=permute(pol,[1 3 2]); 
       % /\ desired denominator shaped as matrix (Fgal)
    end
    Ftf=tf(Fpol,1);
    Fnr=reshape(Fgal,subMat_len,[]); % re-ordering to solve

% creating new Sm matrix of co-prime fractions
      Sm=double(create_Sm_rows(Ngal,Dgal,size(Fnr,2)));  % in order to find compensator we need Sm'  
      dep=rowDependency(Sm); % find linearly dependeices on N block rows
      indepInd=1:size(Sm,1); %  Removing dependent indexes 
      indepInd(dep)=[];      %

      FnonZero=any(Fnr,1)~=0; % nonzero column on F (desired equations)

      Smt=Sm(indepInd,FnonZero); % create a full rank square:
                        % lineary independent rows
                        % non homogenic equations

      W=zeros(subMat_len,size(Sm,1)); % reseting to insert answers correctly
      W(:,indepInd)=Fnr(:,FnonZero)*Smt^-1; % solving system

  
  % solving for A,B
      Compensator_resh=reshape(W,subMat_len,subMat_len,[]);
      Agal=Compensator_resh(:,:,1:2:end);
      Bgal=Compensator_resh(:,:,2:2:end); 
end
  