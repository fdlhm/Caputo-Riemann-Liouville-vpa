% sylvesterND.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
function X=sylvesterND(AA,X)
N=length(AA);
UU=cell(1,N);
TT=cell(1,N);
NN=zeros(1,N);
for j=1:N
    NN(j)=size(AA{j},1);
    [UU{j},TT{j}]=schur(AA{j},'complex');
    X=multND(UU{j}',X,j);
end
cumprodN=cumprod([1 NN(1:end-1)]');
N=length(NN);
ii=NN;
for index=prod(NN):-1:1
    num=X(index);
    den=0;
    for j=1:N
        Tj=TT{j};
        den=den+Tj(ii(j),ii(j));
        for k=ii(j)+1:NN(j)
            num=num-X(index+((k-ii(j))*cumprodN(j)))*Tj(ii(j),k);
        end
    end
    X(index)=num/den;
    if ii(1)>1
        ii(1)=ii(1)-1;
    else
        k=2;
        while k<=N&&ii(k)==1
            k=k+1;
        end
        if k<=N
            ii(k)=ii(k)-1;
            ii(1:k-1)=NN(1:k-1);
        end
    end
end
for j=1:N
    X=multND(UU{j},X,j);
end
