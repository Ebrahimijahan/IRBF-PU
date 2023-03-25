function [A,Ax,Axx] = RBF_PU_Mat(data,center,ep)
N = length(data);
M = length(center);
% For each patch, find which node points that lie inside.
Omega = cell(M,1); indice = cell(M,1);
% for j=1:M
%   indice{j} = find(abs((data(:,1) - center(j,1)))<R);
%   Omega{j} = [data(indice{j},1)];
% end
ns=5;
n=3;
[sten,segment] = GenerateStencils(data(:,1),ns,n);
for j=1:M
    indice{j} = sort(segment(j,:));
    Omega{j} = [data(indice{j},1)];
end
% Generate differentiation matrices
[A,Ax,Axx] = deal(sparse(N,N));
for i=1:N
    % Find the active centers for each node
    % Ix = find(abs(data(i,1) - center(:,1))<R);
    Ix=sort(sten(i,:));
    % pn=sort(sten(i,:));
    activecen = center(Ix,:); %activecen=I(x)
    Omega{Ix};
    mat = PU(activecen,data(i,:));
    w=mat.o;
    wx=mat.x;
    wxx=mat.xx;
    for k=1:length(Ix) %Number of the active centers for each node
        [psi,psix,psixx] = IRBF(ep,Omega{Ix(k)},data(i,:));
        %       IRBF(Omega{Ix(k)},data(i,:),ep)
        A(i,indice{Ix(k)}) = A(i,indice{Ix(k)}) + (w(k)*psi);
        Ax(i,indice{Ix(k)}) = Ax(i,indice{Ix(k)}) + (wx(k)*psi + w(k)*psix);
        Axx(i,indice{Ix(k)}) = Axx(i,indice{Ix(k)}) + (wxx(k)*psi + 2*wx(k)*psix + w(k)*psixx);
    end
end