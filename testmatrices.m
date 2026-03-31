% testmatrices.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via
% shifted Chebyshev polynomials, (2026).
clear
tic
format short
a=0.37; % $\alpha$
N=100; % $N$
dig=100; % Number of digits
T=1.2; % $T$
% $\hat{\mathbf D}^\alpha_t$, $\mathbf D^\alpha_t$ and $t$
[hatDta,Dta,t]=CaputoMatrix(N,a,T,dig);
% $\hat{\mathbf E}^\alpha_t$ and $\mathbf E^\alpha_t$
[hatEta,Eta]=RiemannLiouvilleMatrix(N,a,T,dig);
m=2;
f=exp(1i*m*t); % $f(t)=e^{imt}$
% Exact values of $D^\alpha_t f(t_j)$
if ceil(a)==a % if $\alpha$ is an integer
    Dtaf=(1i*m)^a*exp(1i*m*t);
else
    Dtaf=(1i*m)^a*exp(1i*m*t).*(1-igamma(ceil(a)-a,1i*m*t)/gamma(ceil(a)-a));
end
% Exact values of $I^\alpha_t f(t_j)$
if a == 0 % if $\alpha=0$
    Itaf=exp(1i*m*t);
else
    Itaf=(1i*m)^(-a)*exp(1i*m*t).*(1-igamma(a,1i*m*t)/gamma(a));
end
g=[f;f(N:-1:2)]; % $g_j$
hatg=fft(g)/N; % $\hat g_k$
hatf=hatg(1:N+1); % $\hat f_k$
hatf([1 N+1])=hatf([1 N+1])/2;
M=2*cos(pi*(0:N)'*(0:N)/N); % Generate $M$ using 64-bit precision
M(:,[1 N+1])=M(:,[1 N+1])/2;
M([1 N+1],:)=M([1 N+1],:)/2;
M=M/N;
% Errors in $L^\infty$ norm
norm(hatf-M*f,inf) % error in the approximation of $\hat f$ by $M$
norm(hatDta*hatf-Dtaf,inf) % error without Krasny's filter
norm(hatEta*hatf-Itaf,inf) % error without Krasny's filter
hatf(abs(hatf)<eps)=0; % Krasny's filter
norm(hatDta*hatf-Dtaf,inf) % error with Krasny's filter
norm(hatEta*hatf-Itaf,inf) % error with Krasny's filter
norm(Dta*f-Dtaf,inf)
norm(hatDta*M*f-Dtaf,inf)
norm(Eta*f-Itaf,inf)
norm(hatEta*M*f-Itaf,inf)
toc % Elapsed time