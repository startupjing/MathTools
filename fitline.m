% ESE403 Operation Research
% Final exam: Problem 6(b)
% Jing Lu: jinglu@wustl.edu

% load data points
dat = load('final.csv');
x = dat(:,1);
y = dat(:,2);
% number of data points
n = length(y);

% convert the problem into Ax<=b form
% constraints on RHS
b = [-y; y];
% coefficients in objective
c = [ones(n,1); 0; 0; 0; 0];
temp = [-x x -ones(n,1) ones(n,1)]
% constraints coefficient matrix
A = [-eye(n) temp; -eye(n) -temp];
% use linprog to get solution
soln = linprog(c,A,b);
% compute slope and interp
len = length(soln);
slope = soln(len-3) - soln(len-2);
interp = soln(len-1) - soln(len);
% compute predicted values
pred = slope.*x + interp;

% build string for the line
plotTitle = sprintf('y = %.4fx + %.4f', slope, interp);

% make plot
figure
% plot data points
scatter(x,y,'filled');
hold on;
% plot the line
plot(x,pred,'--','LineWidth',3);
title(plotTitle);
xlabel('x');
ylabel('y');



