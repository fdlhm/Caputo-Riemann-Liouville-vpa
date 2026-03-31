% pde2D.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via
% shifted Chebyshev polynomials, (2026).
clear
tic
a=0.97; % $\alpha$
m=330; % $m$
Nx=16; % $N_x$
Nt=400; % $N_t$
dig=400; % Number of digits
T=2; % $T$
[~,Dta,t]=CaputoMatrix(Nt,a,T,dig);
b=1.4; % Scale factor $b$
[x,DD]=herdif(Nx,2,b);
Dx=DD(:,:,1); % $\mathbf D_x$
Dx2=DD(:,:,2); % $\mathbf D_x^2$
[tt,xx]=ndgrid(t,x); % Two-dimensional grid
u0=exp(-x.^2).'; % $u_0$
u=exp(1i*m*t-xx.^2); % $u$
E=[eye(Nt);zeros(1,Nt)]; % $\mathbf E$
F=[zeros(Nt,Nx);u0]; % $\mathbf F$
G=Dx2+2*x.*Dx+2*eye(Nx);
H=((1i*m)^a*(1-igamma(ceil(a)-a,1i*m*t)/gamma(ceil(a)-a))) ...
    .*exp(1i*m*tt-xx.^2); % $\mathbf H$
A=E.'*Dta*E; % $\mathbf A$
B=-G.'; % $\mathbf B$
C=-E.'*Dta*F+E.'*H; % $\mathbf C$
Uinner=lyap(A,B,-C); % $\mathbf U_{inner}$
U=E*Uinner+F; % $\mathbf U$
norm(u(:)-U(:),inf) % Error in max norm
toc % Elapsed time