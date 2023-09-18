clc; clear; close all;
% params
m = 8; 
n = 6;
kx = 1.02;
ky = 0.98;
lx = 1.0;
ly = 1.0;

% x,y,dx,dy
x = zeros(1,m+1)';
y = zeros(1,n+1)';
dx = zeros(1,m)';
dy = zeros(1,n)';
% X, DX gem
dx1 = lx  * (1-kx) /(1-kx^m);
if (abs(kx - 1) < 1e-9) 
    dx1= lx/m;
end
x(1) = 0.0;
x(m+1) = lx;
dx(1) = dx1;
for i = 2:m
    dx(i) = kx*dx(i-1);
    x(i) = x(i-1) + dx(i-1);d
end

% Y, DY gen
dy1 = ly  * (1-ky) /(1-ky^n);
if (abs(ky - 1) < 1e-9) 
    dy1= ly/n;
end
y(1) = 0.0;
y(n+1) = ly;
dy(1) = dy1;
for i = 2:n
    dy(i) = ky*dy(i-1);
    y(i) = y(i-1) + dy(i-1);
end

% R GEN
dyq = zeros(m*n, 1);
for i = 1:n
    dyq((i-1)*(m)+1:i*(m)) = repmat(dy(i), m, 1);
end
dxq = zeros((n-1)*m, 1); % Preallocate dxq with zeros
index = 1; % Initialize an index variable
for i = 1:n-1
    dxq(index:index+length(dx)-1) = dx; % Fill dxq
    index = index + length(dx); % Update the index
end
R = [diag([dyq; dxq])];


% M GEN
M = zeros((m-1)*n+m*(n-1),1);
ctr = 0;
c = 0;
for j = 1:n
    for i = 2:m
        M((i-1)+c) = 0.5*(dx(i)+dx(i-1));
        ctr = ctr+1;
    end
    c = ctr;
end
c = (m-1)*n;
ctr = 0;
for i = 2:n
    for j = 1:m
        M(c+(i-2)*(m)+j) = 0.5*(dy(i)+dy(i-1));
    end
end
M = diag(M);


% A = (M*L*inv(R));
% spy(round((A-A')*100)/100);