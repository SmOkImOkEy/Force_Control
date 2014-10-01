function [N,flag,bound] = nulls(A,tol)
% calculate the ker(A), some notes:
% 1. primary dependent col must be the last one
% 2. the function returns only one solution

wei=norm(A)*(size(A,2):-1:1)';

Beq=-wei.*A(:,end);
Aeq=A(:,1:end-1);

N=Aeq\Beq;
N=[N; 1; zeros(size(A,2)-size(N,1)-1,1)]./wei;

end