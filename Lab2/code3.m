clc;
clear;
close all;

v=[1,0,0;
    1,1,0;
    1,1,1;
    1,1,0];

[m,n]=size(v);
U=zeros(m,n);
Q=zeros(m,n);

for i=1:n
    vector=v(:,i);
    for j=1:(i-1)
        r=Q(:,j)'*vector;
        vector=vector-r*Q(:,j);
    end
    U(:,i)=vector;
    nv=norm(vector);
    Q(:,i)=vector/nv;
end

disp(v);
disp(U);
disp(Q);