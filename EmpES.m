%alpha=1.7, beta=0, location zero and scale one
alpha=1.7;
beta=0;
mu=0;
scale=1;
q=-5.1519;
N=1000000; X=[];
X = stabgen(N,alpha,beta,scale,mu,1)'; X = mu+scale*X;
use=X(X<q);
empiricalES = mean(use)