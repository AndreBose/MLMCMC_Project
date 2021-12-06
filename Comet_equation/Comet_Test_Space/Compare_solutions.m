%% Compute Matlab solution

clear
clc

mu    = 2;
theta = 3.14;

nref  = 5;

[errors,solutions,femregion,Dati] = C_main2D('Test1', nref, mu, theta);

grid = femregion.coord;
uh   = full(solutions.uh);

matlabtable = [grid, uh];

%% Load FEniCS solution and compare the two

pythontable = table2array(readtable('comet_solution.csv'));

figure(2)
scatter3(pythontable(:,1), pythontable(:,2),pythontable(:,3), 'bo')
hold on
scatter3(matlabtable(:,1), matlabtable(:,2),matlabtable(:,3), 'rx')
legend('fenics solver', 'PoliFEM SD solver')

%% 