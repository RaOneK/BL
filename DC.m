clc; clear; close all;

m = 6; 
n = 6;
D = zeros(m*n,m*n+m*(n-1));
nodes = size(D,1);
links = size(D,2);



for i = 1:nodes
    D(i,i) = 1;
end

ctr = 0;
for i = 1:n % rows
    ctr = ctr + 2;
    D(ctr,ctr-1) = -1;
    for j = 1:m-2 % columns
        ctr = ctr + 1;
        D(ctr,ctr-1) = -1;
    end
end

for i = 1:nodes-m
    D(i,m*n+i) = 1;
end

for i = m+1:nodes
    % fprintf("i=%d.mn+1=%d\n",i,m*(n-1)+i);
    D(i,m*(n-1)+i) = -1;
end

x0=50;
y0=550;
width=700;
height=400;
set(gcf,'position',[x0,y0,width,height])
x = [0 28];
y = [0 16];
clim([-1 1])
imagesc((D))
colorbar


C = zeros(links,(m-1)*(n-1));

% bot vec
ctrV = 0;
ctrN = 0;
for i = 1:m-1
    ctrV = ctrV + 1;
    ctrN = ctrN + 1;
    C(ctrV,ctrN) = 1;
    for j = 1:n-1
        ctrV = ctrV + 1;
        ctrN = ctrN + 1;
        C(ctrV,ctrN) = 1;
    end
end

% top vec
ctrV = m;
ctrN = 0;
for i = 1:m-1
    ctrV = ctrV + 1;
    ctrN = ctrN + 1;
    C(ctrV,ctrN) = -1;
    for j = 1:n-1
        ctrV = ctrV + 1;
        ctrN = ctrN + 1;
        C(ctrV,ctrN) = -1;
    end
end
% 
% left vec
ctrV = m*n;
ctrN = 0;
for i = 1:n-1
    ctrN = ctrN + 1;
    ctrV = ctrV + 1;
    C(ctrV,ctrN) = -1;
    for j = 1:n-1
        ctrN = ctrN + 1;
        ctrV = ctrV + 1;
        C(ctrV,ctrN) = -1;
    end
end

% right vec
ctrV = m*n;
ctrN = -1;
for i = 1:n-1
    ctrN = ctrN + 2;
    ctrV = ctrV + 2;
    C(ctrV,ctrN) = 1;
    for j = 1:n-2
        ctrN = ctrN + 1;
        ctrV = ctrV + 1;
        C(ctrV,ctrN) = 1;
    end
end



%bdry nodes

figure();
x0=50;
y0=550;
width=500;
height=700;
set(gcf,'position',[x0,y0,width,height])
x = [0 28];
y = [0 16];
clim([-1 1])
imagesc((C))
colorbar
