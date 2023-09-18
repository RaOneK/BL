clear; close all; clc;
format short;
figure('Renderer', 'painters', 'Position', [1000 1200 1650 800]);
figure(1);
filename = 'out.gif';
ctr = 1;
hx = 0.1;
hy = 0.01;
L = 50;
x = 2.5:hx:7.5;
y = -5:hy:0;
[X,Y] = meshgrid(x,y);
dt=0.1;
t=1e-5;


% plot3(X,-Y,real(z(X,Y)));
xm = [5.3 -.55];
syms x y;
deltaglob=1;
ctr = 0;
xarr = [];
m = 1;
drawnow; 
while deltaglob>1e-5
    z = @(x,y) -(-(x).^0.5 - y)- (1/t)*(log(-(sin(x)-y-3)) + log(-((x-3*pi/2).^2+y)));
    % z = f(x)  - log(-f1) - log(-f2) -> find min in this surface
    deltaloc = 1;
    % Start Newton
    while deltaloc > 1e-5
        hx = hessian(z,[x,y]); % SYMBOLIC HESSIAN
        grad = gradient(z,[x,y]); % GRADIENT SYMBOLIC
        gradnum = vpa(subs(grad,{x,y},{xm(1),xm(2)})); % NUMERICAL HESSIAN
        hnum = vpa(subs(hx,{x,y},{xm(1),xm(2)})); % NUMERICAL GRADIENT
        xmp1 = xm - (hnum\gradnum)'; % x_(m + 1) = x_m  - inverse(Hessian(x_m)) * gradient(f(x_m))
        % inv(A) * B = A\B 
        deltaloc = max(abs(xmp1-xm)); % CAUCHY CONVERGENCE
        xarr = [xarr xm']; % ARRAY OF MINIMUM POINTS TO PLOT LATER
        fprintf("%3d. x(m) = (%6f,%6f) - minimum; t = %4f\n",ctr,xm(1),xm(2),t); 
        xm = xmp1;
    end
    % end Newton

    tiledlayout(1,2); % KARTINOK 1x2
    %3d plot
    nexttile;
    plot3(X,Y,real(z(X,Y)),"b"); hold on;
    plot3(xarr(1,:),xarr(2,:),z(xarr(1,:),xarr(2,:)),'or'); hold off;
    % 2d plot
    nexttile;
    plot(xarr(1,:),xarr(2,:),'or'); hold on;
    xgrid = 2.5:0.1:7.5;
    plot(xgrid,-(xgrid).^0.5+(xarr(1,length(xarr))^0.5)+xarr(2,length(xarr))); hold on; % y = -sqrt(x) + sqrt(x_i)-y_i
    plot(xgrid,sin(xgrid)-3,"b"); hold on;
    plot(xgrid,-(xgrid-3*pi/2).^2,"b");hold off;
    xlim([2.5 7.5]); % granicy plota
    ylim([-5 0]);
    drawnow; 
    % START SAVING TO GIF FILE (TOOK FROM INTERNET)
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ctr == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    deltaglob = m/t;
    t = t + dt; 
    dt = dt * 2;
    ctr = ctr + 1;
end