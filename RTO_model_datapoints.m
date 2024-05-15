%Script to run Simulink model to create datapoints for RTO
%model. 

clc; clear all; close all;

%%

n_loops = 12;

obj_function_inputs = zeros(n_loops,12);
ndot_H2O_feed = linspace(35, 155, n_loops);

for i = 1:n_loops
    obj_function_inputs(i, :) = model_runner(ndot_H2O_feed(i));
    i
end
