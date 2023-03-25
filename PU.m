function w = PU(ptrial,ptest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global delta   % size of supports of MLS weight function
global h % fill distance
[~, dim]=size(ptrial);
R=delta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = ptrial(:,1) - ptest(:,1);
rr= (dx./R);
for k=1:length(rr)
    r = (dx(k)./R);
    phi(k,1) = (4*r + 1).*(1 - r).^4;
    phix(k,1)=20*dx(k).*(-1+r).^3./R;
    if r==0
        phixx(k,1)=-20/R.^2;
    else
      phixx(k,1)=20*(1./R.^2).*(4./R.*r-1).*(1-1./R.*r).^2;
    end
end
s = sum(phi);
w.o = phi/s;
sx = sum(phix);
w.x = phix/s - phi*sx/s^2;
sxx = sum(phixx);
w.xx = -2*phix*sx/s^2 + phixx/s + phi*(2*sx^2/s^3 - sxx/s^2);
end % end switch


