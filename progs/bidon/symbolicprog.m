clear all
syms mu1 mu2 mu3 gamma1 gamma2 gamma3 lambda
P=(mu1+lambda)^2*(mu2+lambda)^2*(mu3+lambda)^2 ...
    - gamma1*(mu2+lambda)^2*(mu3+lambda)^2 ...    
    - gamma2*(mu3+lambda)^2*(mu1+lambda)^2 ...
     - gamma3*(mu1+lambda)^2*(mu2+lambda)^2;
     


d1 = diff(P,lambda,1);
d2 = diff(P,lambda,2);
d3 = diff(P,lambda,3);
d4 = diff(P,lambda,4);
d5 = diff(P,lambda,5);
d6 = diff(P,lambda,6);
 
lambda=0;
eval(P)
eval(d1)
eval(d2)/factorial(2)
eval(d3)/factorial(3)
eval(d4)/factorial(4)
eval(d5)/factorial(5)
eval(d6)/factorial(6)

