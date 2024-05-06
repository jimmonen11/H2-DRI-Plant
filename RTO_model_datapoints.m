%Script to run Simulink model to create datapoints for RTO
%model. 

clc; clear all; close all;

% n_samp = 500; %number of samples to get points from the simulink model at
% n_dim = 4; %number of dimensions for the model, 2 decision vars
% n_match = 10000; %number of 'matches' for to create low discrepency grid
% 
% winbox = minimax(n_samp, n_dim, n_match);
% 
% winbox(:,1) = scalefun(563,663,winbox(:,1));
% winbox(:,2) = scalefun(0,750, winbox(:,2));

%%

n_loops = 10;

obj_function_inputs = zeros(n_loops,12);
ndot_H2O_feed = linspace(30, 145, n_loops);

for i = 1:n_loops
    obj_function_inputs(i, :) = model_runner(ndot_H2O_feed(i));
    i
end
