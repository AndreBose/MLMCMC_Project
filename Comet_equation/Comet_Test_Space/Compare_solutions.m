%% Compute Matlab solution

clear
clc

mu    = 1;
theta = 1;

nref  = 6;

[errors,solutions,femregion,Dati] = C_main2D('Test1', nref, mu, theta);

grid = femregion.coord;
uh   = full(solutions.uh);

matlabtable = [grid, uh];

%% Load FEniCS solution and compare the two

pythontable = table2array(readtable('comet_solution.csv'));

figure(1)
hold on

scatter3(pythontable(:,1), pythontable(:,2),pythontable(:,3), 'bo')
scatter3(matlabtable(:,1), matlabtable(:,2),matlabtable(:,3), 'rx')
legend('fenics solver', 'PoliFEM SD solver')

%% 