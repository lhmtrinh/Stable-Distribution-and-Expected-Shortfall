
%global parameters
alpha1 = 1.6; %set to 1.5
alpha2 = 1.8; %set to 1.9
beta1 = 0;
beta2 = 0;

%create grid
x_initial = -15;
x_final = 15;
nPoints = 1e3;
dx = (x_final-x_initial)/(nPoints-1);
x = x_initial : dx : x_final;

%method 1: Convolution formula
s = x;

%find the PDFs for the two stable rvs
stable1PDF = @(t) asymstabpdf(t, alpha1, beta1);
stable2PDF = @(t) asymstabpdf(t, alpha1, beta1);
res1 = zeros(size(s));

%use quadrature to calculate the convolution integrand 
for i=1:1:size(s')
    f = @(t) stable1PDF(t).*stable2PDF(s(i)-t); 
    res1(i) = integral(f, x_initial, x_final, 'ArrayValued', true);
end
plot(s, res1, "LineWidth", 2.5)
hold on

%Method 2: Inversion formula
res2 = sum_asymstab(x, alpha1, beta1, alpha2, beta2);
plot(x,res2, "LineWidth", 2.2)
hold on

%Method 3: Simulation
nSimuls = 10e6;
%simulate the first stable rvs
simulStable1 = stabgen(nSimuls, alpha1, beta1);
%simulate the second stable rvs
simulStable2 = stabgen(nSimuls, alpha2, beta2);
%since the two above rvs are independent, we can just sum them
sumSimulStable = simulStable1 + simulStable2;
[res3, xi]= ksdensity(sumSimulStable, 'NumPoints', nSimuls);
plot(xi, res3, "LineWidth", 1.75)
hold on;

%Method 4: conv(.) function
stable1 = asymstabpdf(x, alpha1, beta1);
stable2 = asymstabpdf(x, alpha2, beta2);
%multiply by dx to normalize
res4 = conv(stable1, stable2, 'same') * dx;
plot(x, res4, "LineWidth", 1.5)
hold on;

title('Density of the convolution between S_{(\alpha=1.6,\beta=0)}(c=1,\mu=0) and S_{(\alpha=1.8,\beta=0)}(c=1,\mu=0)')
legend({'Density with Convolution Formula', 'Density with Inversion Formula', 'Density with Simulation', 'Density with MATLAB conv(.)'},'Location', 'northwest', 'NumColumns', 1)
ax = gca;
ax.FontSize = 7;
xlabel('x') 
ylabel('pdf') 
xlim([-15 15])
hold off;
