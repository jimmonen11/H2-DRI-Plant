Project: H2-DRI Plant

Run the H2_DRI_plant_loader.m file to load relevant variables in Matlab workspace. Then can run SOE_H2_DRI_plant.slx simulations. 
The loader file pulls locational marginal price and carbon emissions data from CAISO. ndot_path is the optimal path found from external optimization using GEKKO.
Can change ndot_path to be more economic or environmental. Can also change simulation to do step changes on ndot_path. To utilize on feedback control, simply update the model and remove the 'feedforward' contribution.

There is also a H2_DRI_plant.slx Simulink file but it is not fully completed. This represents a plant using PEM/AEL. 