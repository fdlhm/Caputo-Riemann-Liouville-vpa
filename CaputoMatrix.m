% CaputoMatrix.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via
% shifted Chebyshev polynomials, (2026).
%
% Create the Caputo differentiation matrices
% $\hat{\mathbf D}_t^\alpha$ and $\mathbf D_t^\alpha$.
% N, a, T and dig denote respectively
% $N$, $\alpha$, $\T$ and the number of digits for vpa.
function [hatDta,Dta,t]=CaputoMatrix(N,a,T,dig)
digits(dig); % Set the number of digits for vpa
a_vpa=vpa(a); % vpa version of $\alpha$
ceila=double(ceil(a_vpa)); % $\lceil\alpha\rceil$
T_vpa=vpa(T); % vpa version of $\T$
t=T_vpa*(1+cos(vpa(pi)*(0:N)/N)')/2; % $t\in[0,T]$
C=sym(zeros(N+1)); % $\mathbf C$
C(1,1)=1;
C(1:2,2)=[-1;2];
for j=3:N+1
    C(1:j,j)=[(-1)^(j-1);4*C(1:j-1,j-1)-2*C(2:j,j-1)-C(2:j,j-2)];
end
D1=t.^((ceila:N)-a_vpa); % $\mathbf D_1$
% $\mathbf d_2$
if double(a_vpa==ceil(a_vpa)) % Test whether $\alpha\in\mathbb Z$
    d2=cumprod([prod(1:a_vpa)/T_vpa^a_vpa,...
        ((a_vpa+1):N)./(1:(N-a_vpa))/T_vpa]).';
else
    d2aux=cumprod([1/gamma(1-a_vpa),((1:N)./((1:N)-a_vpa))/T_vpa]).'; 
    d2=d2aux(ceila+1:N+1);
end
hatDta=(D1*(d2.*C(ceila+1:N+1,:))); % $\hat{\mathbf D}_t^\alpha$
M=2*cos(vpa(pi)*(0:N)'*(0:N)/N);
M(:,[1 N+1])=M(:,[1 N+1])/2;
M([1 N+1],:)=M([1 N+1],:)/2;
M=M/N;
Dta=double(hatDta*M); % $\mathbf D_t^\alpha$
hatDta=double(hatDta);
t=double(t);