clc; clear; close all;
figure('Renderer', 'painters', 'Position', [1000 200 800 800]);%plotcoord
t = readtable('/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/StreamFunction.txt', 'ReadVariableNames', false);
n = max(size(t(:,1)));
nt = max(size(t(1,:)));
filename = ['streamfunction_dim' num2str(sqrt(n)) '_it' num2str(nt) '_colored.gif'];
s = t(:,1);
e = t(:,nt);
oneD = sqrt(n);
step = floor(nt/100);
ctr = 1;
i = 1;
iterationFactor = 1.5;
while i < nt
    ctr = ctr + 1;
    if floor(i*iterationFactor)~=i
        i = floor(i*iterationFactor); 
    else
        i = i+1;% added 1 to break out through first iteration
    end
    if i >=nt % last iterand
        i = nt;
    end
    table = t(:,i);
    arr = table2array(table);
    fprintf("%d. %e\n", i, max(arr));
    mat = reshape(arr,oneD,oneD);
    mat = mat';
    contour(mat, 'LineColor','k', 'LevelStep',0.005);
    drawnow;
        % gif creator
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 2
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
end
fprintf("total frames = %d\n",ctr);


