% RiemannLiouvilleMatrix.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via
% shifted Chebyshev polynomials, (2026).
%
% Create the Riemann-Liouville integration matrices
% $\hat{\mathbf E}_t^\alpha$ and $\mathbf E_t^\alpha$.
% N, a, T and dig denote respectively
% $N$, $\alpha$, $\T$ and the number of digits for vpa.
function [hatEta,Eta,t]=RiemannLiouvilleMatrix(N,a,T,dig)
digits(dig); % Set the number of digits for vpa
a_vpa=vpa(a); % vpa version of $\alpha$
T_vpa=vpa(T); % vpa version of $\T$
t=T_vpa*(1+cos(vpa(pi)*(0:N)/N)')/2; % $t\in[0,T]$
C=sym(zeros(N+1)); % $\mathbf C$
C(1,1)=1;
C(1:2,2)=[-1;2];
for j=3:N+1
    C(1:j,j)=[(-1)^(j-1);4*C(1:j-1,j-1)-2*C(2:j,j-1)-C(2:j,j-2)];
end
E1=t.^((0:N)+a_vpa); % $\mathbf E_2$
e2=cumprod([1/gamma(1+a_vpa),...
    ((1:N)./((1:N)+a_vpa))/T_vpa]).'; % $\mathbf e_2$
hatEta=E1*(e2.*C); % $\hat{\mathbf E}_t^\alpha$
M=2*cos(vpa(pi)*(0:N)'*(0:N)/N);
M(:,[1 N+1])=M(:,[1 N+1])/2;
M([1 N+1],:)=M([1 N+1],:)/2;
M=M/N;
Eta=double(hatEta*M); % $\mathbf E_t^\alpha$
hatEta=double(hatEta);
t=double(t);