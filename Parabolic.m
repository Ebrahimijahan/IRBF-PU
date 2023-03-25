%%% This code writen by Dr. Ali Ebrahimi Jahan
clc; 
clear all;
close all;
global h delta % h: fill distance of trial points, delta: support size of MLS weight
global ep % shape parameter in Gaussian weights
%%

N = 100; 
xe = linspace(0,1,N)'; 
data = [xe(:)];
xm=data(:,1);
xc = linspace(0,1,10)';
H = min(diff(xc));
h=abs(xc(2)-xc(1));
delta=2*h; % the size of PU supports
ep = 0.05; % shape parameter in Gaussian weight (experimentally)

xcen = xc(1:end-1) + 0.5*H; 
center = [xe(:)]; 
M = length(center);
delta =1.5*(0.5*H);
R=delta;
[D.o,D.x,D.xx] = IRBF_PU_Mat(data,center,ep);
border = find(data(:,1)==0 | data(:,1)==1);
%% time discrete
T_final=1;
dt=0.01;
T=(0:dt:T_final)';
Nt=length(T);

%% u_t = u_xx + f(t,x)
u_ex=@(t,x) exp(t).*sin(pi*x);
f=@(t,x) exp(t).*(sin(pi*x)+pi^2*sin(pi*x));

I=eye(N);

U=u_ex(0,data);
L=I-0.5*dt*D.xx;
R=I+0.5*dt*D.xx;
L(border,:)=I(border,:);
for n=1:Nt-1
    rhs=R*U+dt*f(T(n)+0.5*dt,data);
    rhs(border)=u_ex(T(n+1),data(border));
    U=L\rhs;
    U(border)=u_ex(T(n+1),data(border));
    Uex=u_ex(T(n+1),data);
    plot(data,U,data,Uex,'k.')
    drawnow
end
norm(U-Uex)