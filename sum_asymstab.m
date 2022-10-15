function [f,F] = sum_asymstab(xvec,a1,b1,c1,d1,a2,b2,c2,d2)
    xl=length(xvec);
    f = zeros(size(xvec));
    F=zeros(xl,1) ; 
    bordertol = 1e-8; 
    lo= bordertol; 
    hi= 1-bordertol ; 
    for loop = 1 :xl
        x=xvec(loop);
        fun_pdf = @(u)fff(u,x,a1,b1,c1,d1,a2,b2,c2,d2,1);
        f(loop) = integral(fun_pdf,lo,hi)/pi;
        fun_cdf = @(u)fff(u,x,a1,b1,c1,d1,a2,b2,c2,d2,0);
        F(loop) = 0.5 -(1/pi) * integral(fun_cdf,lo,hi);
    end
end

function I = fff(uvec,x,a1,b1,c1,d1,a2,b2,c2,d2,dopdf)
    I = zeros(size(uvec));
    for ii =1:length(uvec)
        u=uvec(ii);
        t =(1-u)/u; 
        if a1==1
            cf1 = exp( -abs(t)*c1*( 1 + 1i*b1*(2/pi)*sign( t )*log(t)) +1i*d2*t );
        else 
            cf1 = exp( -((abs( t ) )^a1)*c1^a1 * ( 1 - 1i*b1*sign(t)*tan(pi*a1/2)) +1i*d1*t);
        end
        if a2 ==1
            cf2 = exp( -abs(t)*c2*( 1 + 1i*b2*(2/pi)*sign( t )*log(t)) +1i*d2*t );
        else        
            cf2 = exp( -((abs( t ) )^a2)*c2^a2 * ( 1 - 1i*b2*sign(t)*tan(pi*a2/2)) +1i*d2*t);
        end
        z = exp(-1i*t*x)*cf1*cf2 ; 
        if dopdf==1
            g = real(z); 
        else 
            g=imag(z)./t;
        end
        I(ii)=g*u^(-2); 
    end
end

