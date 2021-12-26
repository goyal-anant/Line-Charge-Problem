tic
clc; clear all; close all;

L = 1;  %length of the wire
a = 0.001;  %radius of the wire
N = 5;  %discretization level
delta = L/N;    %discretization step
y_m = delta/2 : delta : L-delta/2;  %matching points
epsilon_o = 8.87*10^-12;    %permitivity
b = 4*pi*epsilon_o*ones(N,1);   
syms y_prime    %variable of integration

%constructing A and basis function 'g'
for i = 1 : N
    for j = 1 : N
        integrand = @(y_prime) 1./sqrt((y_m(i)-y_prime).^2 + a^2);
        A(i,j) = integral(integrand,(j-1)*delta,j*delta);
        g(i,:) = rectangularPulse((i-1)*delta,i*delta,y_prime);
    end
end

%getting the coefficients by gauss elimination
a_n = (A\b)*10^12; %in picoColoumbs/m

%calculating charge density 'rho' 
for y_prime_v = 1:N
    rho(y_prime_v) = a_n' * g;
end

% plotting rho vs length for a given 'N'
fplot(rho)
axis([0 L 0 10])
toc
