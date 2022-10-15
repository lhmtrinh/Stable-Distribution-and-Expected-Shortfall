function [ES, VaR] = asymstabES(xi,a,b,mu,scale)
    % Get q , the quantile from the S(0 ,1) distribution

    opt=optimset('Display','off','TolX',1e-6) ;
    q=fzero(@stabcdfroot,-6,opt,xi,a,b,scale,mu); 
    VaR=mu+scale*q;

    if ( q == 0)
        t0 = (1 / a) * atan( b * tan( pi * a/2));
        ES = ( ( 2 * gamma((a-1) / a) ) / ( pi - 2*t0) ) * (cos( t0 ) / cos(a*t0 )^(1 / a) );
        return ;
    end
    ES=(scale*Stoy(q,a,b)/xi)+mu;
end


function diff = stabcdfroot (x, xi , a, b,scale,mu)
    [~ , F] = asymstab (x, a, b,scale,mu);
    diff = F - xi ;
end

function S = Stoy(cut,alpha , beta)
    cut = -cut ; % we use a different sign convention
    bbar = -sign(cut) * beta ;
    t0bar = (1 / alpha) * atan( bbar * tan( pi * alpha / 2) ) ;
    % ' beta ==0 ' => ' bbar ==0 ' => ' t0bar==0'
    small = 1e-6; 
    tol = 1; 
    abscut = abs(cut) ; 
    display = 0;
    integ = quadl(@stoyint , -t0bar+small , pi /2-small , tol, display , abscut , alpha , t0bar) ;
    S = alpha / (alpha -1) / pi * abscut * integ ;
end

function I = stoyint ( t , cut , a, t0bar)
    s = t0bar + t ;
    g = sin(a*s - 2* t ) ./ sin(a*s ) - a * cos( t ) .^2./ sin (a*s).^2;
    v = (cos(a*t0bar) ) .^(1 / (a - 1)) .* (cos( t )./ sin (a*s)) .^(a / (a - 1)).* cos(a*s - t ) ./ cos( t ) ;
    term = -(abs(cut) ^(a / (a - 1)) ) ;
    I=g.*exp( term.* v ) ;
end


