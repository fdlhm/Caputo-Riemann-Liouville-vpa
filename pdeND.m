% pdeND.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via
% shifted Chebyshev polynomials, (2026).
clear
tic
d=5; % Number of spatial dimensions
a=0.97; % $\alpha$
m=29; % $m$
Nx=16; % $N_x$
Nt=60; % $N_t$
dig=60; % Number of digits
T=2; % $T$
[~,Dta,t]=CaputoMatrix(Nt,a,T,dig);
b=1.4; % Scale factor $b$
[x,DD]=herdif(Nx,2,b);
Dx=DD(:,:,1); % $\mathbf D_x$
Dx2=DD(:,:,2); % $\mathbf D_x^2$
aux=[{t};repmat({x},d,1)];
ttxx=cell(d+1,1);
[ttxx{:}]=ndgrid(aux{:}); % (d+1)-dimensional grid
xx=cell(d,1);
aux=aux(2:end);
[xx{:}]=ndgrid(aux{:});
u0=1; % Initial data $u_0$
u=exp(1i*m*t); % Exact solution $u$
for j=1:d
    u0=u0.*exp(-xx{j}.^2);
    u=u.*exp(-ttxx{j+1}.^2);
end
E=[eye(Nt);zeros(1,Nt)]; % $\mathbf E$
F=zeros([Nt+1,Nx*ones(1,d)]); % $\mathbf F$
index=[{Nt+1};repmat({':'},d,1)];
F(index{:})=u0;
H=((1i*m)^a*(1-igamma(ceil(a)-a,1i*m*t)/gamma(ceil(a)-a)))...
    .*exp(1i*m*t); % $\mathbf H$
for j=1:d
    H=H.*exp(-ttxx{j+1}.^2);
end
B=-multND1(E.'*Dta,F)+multND1(E.',H); % $\mathbf B$
G=Dx2+2*x.*Dx+2*eye(Nx); % $\mathbf G$
AA=cell(d+1,1); % Store the matrices $\mathbf A_j$
AA{1}=E.'*Dta*E;
for j=1:d
    AA{j+1}=-G;
end
Uinner=sylvesterND(AA,B); % $\mathbf U_{inner}$
U=multND1(E,Uinner)+F; % $\mathbf U$
norm(U(:)-u(:),inf) % Error in max norm
toc % Elapsed time