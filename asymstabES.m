function [ES, VaR] = asymstabES(xi,a,b,scale,mu)
    % Get q , the quantile from the S(0 ,1) distribution

    opt=optimset('Display','off','TolX',1e-12) ;
    diff = @(x)stabcdfroot(x,xi,a,b);
    q=fzero(diff,-6,opt); 
    VaR=mu+scale*q;

    if ( q == 0)
        t0 = (1 / a) * atan( b * tan( pi * a/2));
        ES = ( ( 2 * gamma((a-1) / a) ) / ( pi - 2*t0) ) * (cos( t0 ) / cos(a*t0 )^(1 / a) );
        return ;
    end
    ES=(scale*Stoy(q,a,b)/xi)+mu;
end


function diff = stabcdfroot (x, xi ,a, b)
    [~ , F] = asymstab(x, a, b);
    diff = F - xi ;
end

function S = Stoy(cut,alpha , beta)
    cut = -cut ; % we use a different sign convention
    bbar = -sign(cut) * beta ;
    t0bar = (1 / alpha) * atan( bbar * tan( pi * alpha / 2) ) ;
    % ' beta ==0 ' => ' bbar ==0 ' => ' t0bar==0'
    small = 1e-6; 
    abscut = abs(cut) ; 
    fun = @(t)stoyint(t,abscut,alpha,t0bar);
    integ = integral(fun , -t0bar+small , pi /2-small,'RelTol',0,'AbsTol',1e-12) ;
    S = alpha / (alpha -1) / pi * abscut * integ ;
end

function I = stoyint ( t , cut , a, t0bar)
    s = t0bar + t ;
    g = sin(a*s - 2* t ) ./ sin(a*s ) - a * cos( t ) .^2./ sin (a*s).^2;
    v = (cos(a*t0bar) ) .^(1 / (a - 1)) .* (cos( t )./ sin (a*s)) .^(a / (a - 1)).* cos(a*s - t ) ./ cos( t ) ;
    term = -(abs(cut) ^(a / (a - 1)) ) ;
    I=g.*exp( term.* v ) ;
end

% function [f,F] = asymstab1(xvec,a,b)
%     bordertol=1e-8; lo=bordertol; hi=1-bordertol; tol=1e-7;
%     xl=length(xvec); F=zeros(xl,1); f=F;
%     for loop=1:length(xvec)
%     x=xvec(loop); dopdf=1;
%     f(loop)= quadl(@fff,lo,hi,tol,[],x,a,b,1) / pi;
%     if nargout>1
%         F(loop)=0.5-(1/pi)* quadl(@fff,lo,hi,tol,[],x,a,b,0);
%     end
%     end
% end
% 
% function I=fff(uvec,x,a,b,dopdf)
%     for ii=1:length(uvec)
%     u=uvec(ii); t = (1-u)/u;
%     if a==1
%     cf = exp( -abs(t)*( 1 + i*b*(2/pi)*sign(t) * log(t) ) );
%     else
%     cf = exp( - ((abs(t))^a) *( 1 - i*b*sign(t) * tan(pi*a/2) ) );
%     end
%     z = exp(-i*t*x) .* cf;
%     if dopdf==1, g=real(z); else g=imag(z)./t; end
%     I(ii) = g / u^2;
%     end
% end