A = [1, 2, 3, 4, 5, 4, 6, 8, 10, 9, 12, 15, 16, 20, 25 ];

x = [0.5,0.4,-1.3,2.8,-0.75]';

n = 5;
y =zeros(n,1);

for j = 1:n
    for i = 1:j-1
        y(i) = y(i) + A((i-1)*n - i*(i-1)/2 + j)*x(j);
    end
    for i = j:n
        y(i) = y(i) + A((j-1)*n - j*(j-1)/2 +i)*x(j);
    end
end

y'


%%%%%%%%%%%%%%%% BFGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xnext = [0.5, 0.15, 2.25]';
x = [0.2, 0.3, 0.4]';
gnext = [-0.75, 1.45, -0.35]';
g = [0.25, -0.35, 0.12]';

B = eye(3);

s = xnext - x
y = gnext - g
%BFGS update
B = B - B*s*s'*B/(s'*B*s) + y*y'/(y'*s)
%Bs = B*s
%dy1 = Bs*Bs'/(s'*Bs)
%dy2 = y*y'/(y'*s)
%B = B -  dy1+dy2 
% 


%%%%%%%% damped BFGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = eye(3); % reset of B
theta = 0;
if s'*y >= 0.2*s'*B*s
    theta = 1;
else
    theta = 0.8*s'*B*s/(s'*B*s - s'*y)
end

r = theta*y + (1-theta)*B*s;

B = B - B*s*s'*B/(s'*B*s) + r*r'/(s'*r)
    

%%%%%%%%%%%%%%%%%%% inverse Hessian approximation %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('===============================================================')
disp('Approximation of INVERSE of Hessian (should be symmetric as well):')
I = eye(3);
H = y'*s/(y'*y)*I
rho = 1/(y'*s);

proleft = I - rho*s*y'   %non-symmetric
proright = I - rho*y*s'    %proleft'
H = proleft*H*proright  + rho*s*s'


