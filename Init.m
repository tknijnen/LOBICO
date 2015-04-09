%% paths
restoredefaultpath;
if isunix
    addpath(genpath('/users/tknijnen/cplex/ILOG/CPLEX_Studio124/cplex/matlab'));
    addpath('/users/tknijnen/cplex/ILOG/CPLEX_Studio124/cplex/examples/src/matlab');
else
    addpath(genpath(('C:\Program Files\IBM\ILOG\CPLEX_Studio125\cplex\matlab\x64_win64')));
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio125\cplex\examples\src\matlab');
end
addpath(genpath(pwd));

%% Variables and Figs
clc; clear all; close all;