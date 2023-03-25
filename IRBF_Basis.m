function Phi = IRBF_Basis(kind,Dim)
switch Dim
    case 1
        switch kind
            case 'Gaussian'
                Phi.o=@(r,c) exp(-c.^2.*r.^2);
                Phi.Ir=@(r,c) sqrt(pi).*erf(c.*r)./(2.*c);
                Phi.IIr=@(r,c)(erf(c.*r).*c.*r.*sqrt(pi) ...
                    +exp(-c.^2.*r.^2))./(2.*c.^2);
                Phi.IIIr=@(r,c) (2.*exp(-c.^2.*r.^2).*r.*c ...
                    +2.*erf(c.*r).*(c.^2.*r.^2+1./2).*sqrt(pi))./(8.*c.^3);
                Phi.IIIIr=@(r,c) ((2.*c.^2.*r.^2+2).*exp(-c.^2.*r.^2) ...
                    +2.*erf(c.*r).*c.*sqrt(pi).*(c.^2.*r.^2+3./2).*r)./(24.*c.^4);
            case 'IMq'
                Phi.o=@(r,c) 1./sqrt(c.^2.*r.^2+1);
                Phi.Ir=@(r,c) log(c.*r+sqrt(c.^2.*r.^2+1))./c;
                Phi.IIr=@(r,c)(log((c.^2.*r+sqrt(c.^2.*r.^2+1).*sqrt(c.^2)) ...
                    ./sqrt(c.^2)).*r-sqrt(c.^2).*sqrt(c.^2.*r.^2+1)./c.^2)./sqrt(c.^2);
            case 'GIMq'
                Phi.o=@(r,c)  1./(c.^2.*r.^2+1).^2;
                Phi.Ir=@(r,c)  r./(2.*(c.^2.*r.^2+1))+atan(c.*r)./(2.*c);
                Phi.IIr=@(r,c)  r.*atan(c.*r)./(2.*c);
                Phi.IIIr=@(r,c) (atan(c.*r).*c.^2.*r.^2-c.*r+atan(c.*r))./(4.*c.^3);
                Phi.IIIIr=@(r,c) (r.^3.*atan(c.*r).*c.^3-2.*c.^2.*r.^2 ...
                    +3.*r.*atan(c.*r).*c-log(c.^2.*r.^2+1))./(12.*c.^4);
            case 'Iq'
                Phi.o=@(r,c)  1./(c.^2.*r.^2+1);
                Phi.Ir=@(r,c)  atan(c.*r)./c;
                Phi.IIr=@(r,c)  (2.*r.*atan(c.*r).*c-log(c.^2.*r.^2+1))./(2.*c.^2);
                Phi.IIIr=@(r,c) (atan(c.*r).*c.^2.*r.^2 ...
                    -r.*log(c.^2.*r.^2+1).*c+c.*r-atan(c.*r))./(2.*c.^3);
                Phi.IIIIr=@(r,c) (1./(12.*c.^4)).*(2.*r.^3.*atan(c.*r).*c.^3 ...
                    -3.*log(c.^2.*r.^2+1).*r.^2.*c.^2+5.*c.^2.*r.^2 ...
                    -6.*r.*atan(c.*r).*c+log(c.^2.*r.^2+1)+3);
            case 'Matern'
                Phi.o=@(r,c)  exp(-c.*sqrt(r.^2));
                Phi.Ir=@(r,c)  (1-exp(-c.*r))./c;
                Phi.IIr=@(r,c)  (c.*r+exp(-c.*r))./c.^2;
            case 'LMatern'
                Phi.o=@(r,c)  exp(-c.*sqrt(r.^2)).*(c.*sqrt(r.^2)+1);
                Phi.Ir=@(r,c)  (-exp(-c.*r).*c.*r-2.*exp(-c.*r)+2)./c;
                Phi.IIr=@(r,c)  (exp(-c.*r).*c.*r+2.*c.*r+3.*exp(-c.*r))./c.^2;
            case 'QMatern'
                Phi.o=@(r,c)  exp(-c.*sqrt(r.^2)).*(3+3.*c.*sqrt(r.^2)+c.^2.*r.^2);
                Phi.Ir=@(r,c)  (8+(-c.^2.*r.^2-5.*c.*r-8).*exp(-c.*r))./c;
                Phi.IIr=@(r,c)  ((c.^2.*r.^2+7.*c.*r+15).*exp(-c.*r)+8.*c.*r)./c.^2;
            case 'CMatern'
                Phi.o=@(r,c)  exp(-c.*sqrt(r.^2)).*(15+15.*c.*sqrt(r.^2) ...
                    +6.*c.^2.*r.^2+c.^3.*(r.^2).^(3./2));
                Phi.Ir=@(r,c)  (48+(-c.^3.*r.^3-9.*c.^2.*r.^2 ...
                    -33.*c.*r-48).*exp(-c.*r))./c;
                Phi.IIr=@(r,c)  ((c.^3.*r.^3+12.*c.^2.*r.^2+57.*c.*r ...
                    +105).*exp(-c.*r)+48.*c.*r)./c.^2;
            case 'Mq'
                Phi.o=@(r,c)  sqrt(c.^2.*r.^2+1);
                Phi.Ir=@(r,c)  (r.*sqrt(c.^2.*r.^2+1).*c ...
                    +log(c.*r+sqrt(c.^2.*r.^2+1)))./(2.*c);
                Phi.IIr=@(r,c)  (sqrt(c.^2.*r.^2+1).*c.^2.*r.^2+3.*log(c.*r ...
                    +sqrt(c.^2.*r.^2+1)).*r.*c-2.*sqrt(c.^2.*r.^2+1))./(6.*c.^2);
            case 'MQ'
                Phi.o=@(r,c)  sqrt(c.^2+r.^2);
                Phi.Ir=@(r,c)  (1./2).*r.*sqrt(c.^2+r.^2) ...
                    +(1./2).*c.^2.*log(r+sqrt(c.^2+r.^2));
                Phi.IIr=@(r,c)  (1./2).*log(r+sqrt(c.^2+r.^2)).*r.*c.^2 ...
                    +(1./6).*(-2.*c.^2+r.^2).*sqrt(c.^2+r.^2);
            case 'Linear'
                Phi.o=@(r,c)  r;
                Phi.Ir=@(r,c)  (1/2).*r.^2;
                Phi.IIr=@(r,c)  (1/6).*r.^3;
                Phi.IIIr=@(r,c) (1/24).*r.^4;
                Phi.IIIIr=@(r,c) (1/120).*r.^5;
            case 'TPS'
                Phi.o=@(r,c)  c.^2.*r.^2;
                Phi.Ir=@(r,c)  (1/3)*c.^2.*r.^3;
                Phi.IIr=@(r,c) (1/12)*c.^2.*r.^4;
                Phi.IIIr=@(r,c) (1/60)*c.^2.*r.^5;
        end
    case 2
        switch kind
            case 'Gaussian'
                Phi.o=@(rx,ry,c) exp(-c.^2.*(rx.^2+ry.^2));
                Phi.Ix=@(rx,ry,c) exp(-c.^2.*ry.^2).*sqrt(pi).*erf(c.*rx)./(2.*c);
                Phi.Iy=@(rx,ry,c) exp(-c.^2.*rx.^2).*sqrt(pi).*erf(c.*ry)./(2.*c);
                Phi.IIx=@(rx,ry,c) (erf(c.*rx).*c.*rx.*sqrt(pi) ...
                    +exp(-c.^2.*rx.^2)).*exp(-c.^2.*ry.^2)./(2.*c.^2);
                Phi.IIy=@(rx,ry,c) (erf(c.*ry).*c.*ry.*sqrt(pi) ...
                    +exp(-c.^2.*ry.^2)).*exp(-c.^2.*rx.^2)./(2.*c.^2);
            case 'IMq'
                Phi.o=@(rx,ry,c) 1./sqrt(1+c.^2.*(rx.^2+ry.^2));
                Phi.Ix=@(rx,ry,c) log(c.*rx+sqrt(1+c.^2.*(rx.^2+ry.^2)))./c;
                Phi.Iy=@(rx,ry,c) log(c.*ry+sqrt(1+c.^2.*(rx.^2+ry.^2)))./c;
                Phi.IIx=@(rx,ry,c) (log(c.*rx+sqrt(1+c.^2.*(rx.^2+ry.^2))).*rx.*c ...
                    -sqrt(1+c.^2.*(rx.^2+ry.^2)))./c.^2;
                Phi.IIy=@(rx,ry,c) (log(c.*ry+sqrt(1+c.^2.*(rx.^2+ry.^2))).*ry.*c ...
                    -sqrt(1+c.^2.*(rx.^2+ry.^2)))./c.^2;
            case 'GIMq'
                Phi.o=@(rx,ry,c) 1./(1+c.^2.*(rx.^2+ry.^2)).^2;
                Phi.Ix=@(rx,ry,c) rx./((2.*(c.^2.*ry.^2+1)).*(c.^2.*rx.^2+c.^2.*ry.^2+1)) ...
                    +atan(c.*rx./sqrt(c.^2.*ry.^2+1))./(2.*c.*(c.^2.*ry.^2+1).^(3./2));
                Phi.Iy=@(rx,ry,c) ry./((2.*(c.^2.*rx.^2+1)).*(c.^2.*rx.^2+c.^2.*ry.^2+1)) ...
                    +atan(c.*ry./sqrt(c.^2.*rx.^2+1))./(2.*c.*(c.^2.*rx.^2+1).^(3./2));
                Phi.IIx=@(rx,ry,c) (2.*rx.*atan(c.*rx./sqrt(c.^2.*ry.^2+1)).*c ...
                    +sqrt(c.^2.*ry.^2+1).*log(c.^2.*ry.^2+1))./(4.*c.^2.*(c.^2.*ry.^2+1).^(3./2));
                Phi.IIy=@(rx,ry,c) (2.*ry.*atan(c.*ry./sqrt(c.^2.*rx.^2+1)).*c ...
                    +sqrt(c.^2.*rx.^2+1).*log(c.^2.*rx.^2+1))./(4.*c.^2.*(c.^2.*rx.^2+1).^(3./2));
            case 'Iq'
                Phi.o=@(rx,ry,c) 1./(1+c.^2.*(rx.^2+ry.^2));
                Phi.Ix=@(rx,ry,c) atan(c.*rx./sqrt(c.^2.*ry.^2+1))./(c.*sqrt(c.^2.*ry.^2+1));
                Phi.Iy=@(rx,ry,c) atan(c.*ry./sqrt(c.^2.*rx.^2+1))./(c.*sqrt(c.^2.*rx.^2+1));
                Phi.IIx=@(rx,ry,c) (rx.*atan(c.*rx./sqrt(c.^2.*ry.^2+1)) ...
                    -sqrt(c.^2.*ry.^2+1).*log(c.^2.*rx.^2./(c.^2.*ry.^2+1)+1)./(2.*c))./(c.*sqrt(c.^2.*ry.^2+1));
                Phi.IIy=@(rx,ry,c) (ry.*atan(c.*ry./sqrt(c.^2.*rx.^2+1)) ...
                    -sqrt(c.^2.*rx.^2+1).*log(c.^2.*ry.^2./(c.^2.*rx.^2+1)+1)./(2.*c))./(c.*sqrt(c.^2.*rx.^2+1));
            case 'Matern'
                Phi.o=@(rx,ry,c) exp(-c.^2.*(rx.^2+ry.^2));
                Phi.Ix=@(rx,ry,c) exp(-c.^2.*ry.^2).*sqrt(pi).*erf(c.*rx)./(2.*c);
                Phi.Iy=@(rx,ry,c) exp(-c.^2.*rx.^2).*sqrt(pi).*erf(c.*ry)./(2.*c);
                Phi.IIx=@(rx,ry,c) (erf(c.*rx).*c.*rx.*sqrt(pi) ...
                    +exp(-c.^2.*rx.^2)).*exp(-c.^2.*ry.^2)./(2.*c.^2);
                Phi.IIy=@(rx,ry,c) (erf(c.*ry).*c.*ry.*sqrt(pi) ...
                    +exp(-c.^2.*ry.^2)).*exp(-c.^2.*rx.^2)./(2.*c.^2);
            case 'Mq'
                Phi.o=@(rx,ry,c) sqrt(1+c.^2.*(rx.^2+ry.^2));
                Phi.Ix=@(rx,ry,c) (1./(2.*c)).*(log(c.*rx+sqrt(1+c.^2.*(rx.^2+ry.^2))).*c.^2.*ry.^2 ...
                    +rx.*sqrt(1+c.^2.*(rx.^2+ry.^2)).*c+log(c.*rx+sqrt(1+c.^2.*(rx.^2+ry.^2))));
                Phi.Iy=@(rx,ry,c) (1./(2.*c)).*(log(c.*ry+sqrt(1+c.^2.*(rx.^2+ry.^2))).*c.^2.*rx.^2 ...
                    +c.*ry.*sqrt(1+c.^2.*(rx.^2+ry.^2))+log(c.*ry+sqrt(1+c.^2.*(rx.^2+ry.^2))));
                Phi.IIx=@(rx,ry,c) (1./(6.*c.^2)).*(3.*c.*rx.*(c.^2.*ry.^2+1).*log(c.*rx ...
                    +sqrt(1+c.^2.*(rx.^2+ry.^2)))+(-2+(rx.^2-2.*ry.^2).*c.^2).*sqrt(1+c.^2.*(rx.^2+ry.^2)));
                Phi.IIy=@(rx,ry,c) (1./(6.*c.^2)).*(3.*c.*ry.*(c.^2.*rx.^2+1).*log(c.*ry ...
                    +sqrt(1+c.^2.*(rx.^2+ry.^2)))+(-2+(-2.*rx.^2+ry.^2).*c.^2).*sqrt(1+c.^2.*(rx.^2+ry.^2)));
            case 'MQ'
                Phi.o=@(rx,ry,c) sqrt(c.^2+rx.^2+ry.^2);
                Phi.Ix=@(rx,ry,c) (1./2).*(c.^2+ry.^2).*log(rx ...
                    +sqrt(c.^2+rx.^2+ry.^2))+(1./2).*rx.*sqrt(c.^2+rx.^2+ry.^2);
                Phi.Iy=@(rx,ry,c) (1./2).*(c.^2+rx.^2).*log(ry ...
                    +sqrt(c.^2+rx.^2+ry.^2))+(1./2).*ry.*sqrt(c.^2+rx.^2+ry.^2);
                Phi.IIx=@(rx,ry,c) (1./2).*rx.*(c.^2+ry.^2).*log(rx+sqrt(c.^2 ...
                    +rx.^2+ry.^2))-(1./3).*(c.^2-(1./2).*rx.^2+ry.^2).*sqrt(c.^2+rx.^2+ry.^2);
                Phi.IIy=@(rx,ry,c) (1./2).*ry.*(c.^2+rx.^2).*log(ry+sqrt(c.^2 ...
                    +rx.^2+ry.^2))-(1./3).*sqrt(c.^2+rx.^2+ry.^2).*(c.^2+rx.^2-(1./2).*ry.^2);
            case 'TPS'
                Phi.o=@(rx,ry,c) c.^2.*(rx.^2+ry.^2);
                Phi.Ix=@(rx,ry,c) ((1./3).*rx.^3+ry.^2.*rx).*c.^2;
                Phi.Iy=@(rx,ry,c) c.^2.*(rx.^2.*ry+(1./3).*ry.^3);
                Phi.IIx=@(rx,ry,c) ((1./12).*rx.^4+(1./2).*ry.^2.*rx.^2).*c.^2;
                Phi.IIy=@(rx,ry,c) c.^2.*((1./2).*ry.^2.*rx.^2+(1./12).*ry.^4);
            case 'Linear'
                Phi.o=@(rx,ry,c) c.^2+rx.^2+ry.^2;
                Phi.Ix=@(rx,ry,c) c.^2.*rx+(1./3).*rx.^3+ry.^2.*rx;
                Phi.Iy=@(rx,ry,c) c.^2.*ry+rx.^2.*ry+(1./3).*ry.^3;
                Phi.IIx=@(rx,ry,c) (1./2).*c.^2.*rx.^2+(1./12).*rx.^4+(1./2).*ry.^2.*rx.^2;
                Phi.IIy=@(rx,ry,c) (1./2).*c.^2.*ry.^2+(1./2).*ry.^2.*rx.^2+(1./12).*ry.^4;
        end
end
end
