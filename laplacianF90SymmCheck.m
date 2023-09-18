clc; clear; close all;
% figure('Renderer', 'painters', 'Position', [1000 200 800 800]);%plotcoord

Ltable = readtable('/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/L.txt', 'ReadVariableNames', false);
L = table2array(Ltable);

% dxtable = readtable('LaplacianCheck/Xgrid.txt', 'ReadVariableNames', false);
% dx = table2array(dxtable);
% 
% dytable = readtable('LaplacianCheck/Ygrid.txt', 'ReadVariableNames', false);
% dy = table2array(dytable);

filename = '/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/Xgrid.txt'; % specify the name of your text file
fileID = fopen(filename, 'r'); % open the file for reading
x = fscanf(fileID, '%f'); % read the data from the file and store it in a vector
fclose(fileID); % close the file

filename = '/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/Ygrid.txt'; % specify the name of your text file
fileID = fopen(filename, 'r'); % open the file for reading
y = fscanf(fileID, '%f'); % read the data from the file and store it in a vector
fclose(fileID); % close the file

dx = diff(x);
dy = diff(y);
n = length(dy);
m = length(dx);

Dtable = readtable('/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/D.txt', 'ReadVariableNames', false);
D = table2array(Dtable);

Ctable = readtable('/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/C.txt', 'ReadVariableNames', false);
C = table2array(Ctable);

for i = 1:n
    D(i*m,(n-1)*i) = -1;
end
null(D,'rational')
sum(sum(abs(D*C)))

BCinhom = zeros(m*n,1);
% CHECK Du=bc2
for i = 0:n-1
    BCinhom(i*(m)+1) = 1;
end
% 
% for i = 1:n
%     D(i*m,(n-1)*i) = 0;
% end



X = lsqminnorm(D,BCinhom);
max(abs(D*X-BCinhom))



% MAR^-1
dx = diff(x);
dy = diff(y);

% R GEN
n = length(dy);
dyq = zeros(n*(n-1), 1);
for i = 1:n
    dyq((i-1)*(n-1)+1:i*(n-1)) = repmat(dy(i), n-1, 1);
end
m = length(dx);
dxq = [];
for i = 1:m-1
    dxq = [dxq; dx];
end
R = [diag([dyq; dxq])];

% M GEN
m = length(dx);
n = length(dy);
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


A = (M*L*inv(R));
spy(round((A-A')*100)/100);




% CTMACtable = readtable('/Users/sarabek/Library/Mobile Documents/com~apple~CloudDocs/UofA/RAUAN/proj2/code/HALL_METHOD/CTAC.txt', 'ReadVariableNames', false);
% CTMAC = table2array(Ltable);
% spy(CTMAC'-CTMAC)