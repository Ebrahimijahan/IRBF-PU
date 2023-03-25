function [w,wx,wxx] = IRBF(ep,x,xc)
dc=0.05;
minK=1e1;
maxK=1e14;
phi= IRBF_Basis('MQ',1);
N=length(x);
Nx=floor(sqrt(N));
%   minK  min condition number of the irbf matrix, e.g. 1e+2
%   maxK  max condition number of the irbf matrix, e.g. 1e+12
%     dc  shape parameter increment, e.g. 0.1
K = 1;
%============================
% [Columnx,Columny] = ExtraMat(x+y,x+y);
%============================
rx=x-x';
%============================
while (K<minK || K>maxK)
    A = [phi.IIr(rx,ep),x,x.^0];
    K=cond(A);
    if K<minK,   ep = ep+dc;
    elseif K>maxK, ep =ep-dc;
    end
end
%============================
I=[phi.IIr(rx,ep),x,x.^0];
I_inv=pinv(I);
ip=find(xc==x);
dx = xc(1,1) - x';
Ae = [phi.IIr(dx,ep),xc,1];
w = Ae*I_inv;
Aex =[phi.Ir(dx,ep),1,0];
wx = Aex*I_inv;
Aexx =[phi.o(dx,ep),0,0];
wxx = Aexx*I_inv;
%============================
end