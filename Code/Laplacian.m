function [L]=Laplacian(simMatrix)
    [m,~]=size(simMatrix);
    D=zeros(m,m);
    for i=1:m
        D(i,i)=sum(simMatrix(i,:));
    end    
    L=D-simMatrix;
end