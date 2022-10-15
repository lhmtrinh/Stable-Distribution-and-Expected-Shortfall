function [ES, VaR] = asymstableES(xi, a, b, mu, scale, method)
if nargin < 3 , b=0; end , if nargin < 4 , mu=0; end
if nargin < 5, scale = 1; end , if nargin < 6 , method = 1; end

% Get q , the quantile from the S(0,1) distribution
opt=optimset('Display', 'off', 'TolX' , 1e-6) ;
q=fzero(@stabcdfroot, -6, opt, xi, a, b); VaR=mu+scale*q;


if (q == 0)
t0 = (1/a)*atan(b*tan(pi*a/2));
ES = ((2*gamma((a-1)/a))/(pi-2*t0) )*(cos(t0)/cos(a*t0)^(1/a));
return ;
end
if (method==1), ES=(scale*Stoy(q,a,b)/xi)+mu;
else ES=( scale*stabletailcomp(q,a,b)/xi)+mu;
end


function diff = stabcdfroot (x, xi, a, b)
if exist('stableqkcdf.m','file') , F = stableqkcdf(x, [a, b],1) ; % Nolan routine
else [~ , F] = asymstab (x, a, b);
end
diff = F - xi;


function tailcomp = stabletailcomp(q, a, b )
% direct integration of x∗ f ( x ) and use of asymptotic tail behavior
K = (a/pi)*sin( pi*a/2)*gamma(a)*(1-b); % formula f o r K_
ell = -120; M = ell ; display = 0; term1=K*(-M)^(1-a)/(1-a);
term3 = quadl(@stableCVARint, ell, q, 1e-5, display, a, b) ;
tailcomp = term1 + term3;


function [g] = stableCVARint(x, a, b)
if exist('stableqkpdf.m' ,'file'), den = stableqkpdf(x, [a, b],1) ;
else den = asymstab(x, a, b);
end
g = x.* den;

function S = Stoy(cut, alpha, beta)
if nargin<3 , beta=0; end
cut = -cut ; % we use a different sign convention
bbar = -sign(cut)*beta ;
t0bar = (1/alpha)*atan(bbar*tan( pi*alpha/2));
% ' beta ==0 ' => ' bbar ==0 ' => ' t0bar==0'
small = 1e-8; tol = 1e-8; abscut = abs(cut) ; display = 0;
integ = quadl(@stoyint, -t0bar+small, pi/2-small, tol, ...
            display, abscut, alpha, t0bar) ;
S = alpha / (alpha-1) / pi*abscut*integ ;


function I = stoyint(t, cut, a, t0bar)
s = t0bar + t ;
g = sin(a*s-2*t )./sin(a*s)-a*cos(t).^2./sin(a*s).^2;
v = (cos(a*t0bar)).^(1 / (a-1)).*(cos(t)./sin(a*s)).^(a/(a-1)) ...
        .*cos(a*s-t )./cos(t);
term = -(abs(cut)^(a/(a-1)));
I=g.*exp(term.*v);


%alpha=1.7, beta=0, location zero and scale one