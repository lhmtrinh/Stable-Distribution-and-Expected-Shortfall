function [f,F] = asymstab (xvec,a,b,c,d)
    % Calculate the pdf and cdf of a asymetric stable distribution
    % alpha = a, beta = b
    
    if nargin<4
        c=1;
        d=0;
    end
    
    bordertol = 1e-8; 
    lo= bordertol; 
    hi= 1-bordertol ; 
    xl=length(xvec); 
    F=zeros(xl,1) ; 
    f=F;
    for loop=1:length(xvec)
        x=xvec(loop); 
        fun_pdf = @(u)fff(u,x,a,b,c,d,1);
        fun_cdf = @(u)fff(u,x,a,b,c,d,0);
        f(loop) = integral(fun_pdf,lo,hi,'RelTol',0,'AbsTol',1e-12)/pi ;
        
        F(loop) = 0.5 -(1/pi) * integral(fun_cdf,lo,hi,'RelTol',0,'AbsTol',1e-12);
    end
end

function I = fff(uvec,x,a,b,c,d,dopdf)
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
        if dopdf==1
            g=real(z) ; 
        else
            g=imag(z)./t;
        end
        I(ii)=g*u^(-2);
    end
end