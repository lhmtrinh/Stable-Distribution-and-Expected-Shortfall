function [f,F] = asymstab (xvec,a,b,c,d)
    % Calculate the pdf and cdf of a asymetric stable distribution
    % alpha = a, beta = b

    bordertol = 1e-8; 
    lo= bordertol; 
    hi= 1-bordertol ; 
    tol =1e-7;
    xl=length(xvec); 
    F=zeros(xl,1) ; 
    f=F;
    for loop=1:length(xvec)
        x=xvec(loop); 
        f(loop) = quadl(@fff ,lo,hi,tol,[],x,a,b,c,d) / pi ;
        F(loop) = 0.5 -(1/pi) * quadl(@fff,lo,hi,tol,[],x,a,b,c,d);
    end
end

function I = fff(uvec,x,a,b,c,d)
    I = zeros(size(uvec));
    for ii =1:length(uvec)
        u=uvec(ii);
        t =(1-u)/u;

        if a==1
            cf = exp( -abs(t)*c*( 1 + 1i*b*(2/pi)*sign( t )*log(t)) +1i*d*t ) ;
        else 
            cf = exp( -((abs( t ) )^a)*c^a * ( 1 - 1i*b*sign(t)*tan(pi*a/2)) +1i*d*t);
        end
        z = exp(-1i*t*x).*cf ; 
        g=real(z) ; 
        I(ii)=g*u^(-2);
    end
end