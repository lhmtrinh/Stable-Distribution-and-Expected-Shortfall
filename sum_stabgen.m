function x=sum_stabgen(nobs,a1,b1,c1,d1,a2,b2,c2,d2,seed1,seed2)
    z=nobs;
    rand('twister', seed1), V1= unifrnd(-pi/2, pi/2, 1, z);
    rand('twister', seed1+42), W1= exprnd(1, 1, z);
    if a1==1
    x1=(2/pi)*(((pi/2)+b1*V1).*tan(V1)-b1*log((W1.*cos(V1))./((pi/2)+b1*V1)));
    x1=c1*x1+d1-(2/pi)*d1*log(d1)*c1*b1;
    else
    Cab=atan(b1*tan(pi*a1/2))/(a1); Sab=(1+b1^2*(tan((pi*a1)/2))^2)^(1/(2*a1));
    A=(sin(a1*(V1+Cab)))./((cos(V1)).^(1/a1));
    B0=(cos(V1-a1*(V1+Cab)))./W1; B=(abs(B0)).^((1-a1)/a1);
    x1=Sab*A.*(B.*sign(B0)); x1=c1*x1+d1;
    end
    rand('twister', seed2), V2= unifrnd(-pi/2, pi/2, 1, z);
    rand('twister', seed2+42), W2= exprnd(1, 1, z);
    if a2==1
    x2=(2/pi)*(((pi/2)+b2*V2).*tan(V2)-b2*log((W2.*cos(V2))./((pi/2)+b2*V2)));
    x2=c2*x2+d2-(2/pi)*d2*log(d2)*c2*b2;
    else
    Cab=atan(b2*tan(pi*a2/2))/(a2); Sab=(1+b2^2*(tan((pi*a2)/2))^2)^(1/(2*a2));
    A=(sin(a2*(V2+Cab)))./((cos(V2)).^(1/a2));
    B0=(cos(V2-a2*(V2+Cab)))./W2; B=(abs(B0)).^((1-a2)/a2);
    x2=Sab*A.*(B.*sign(B0)); x2=c2*x2+d2;
    end
    x=x1+x2;
end

