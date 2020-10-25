function [all_steps] = RecallOdeFcn(fun, T, Y, varargin)
%{

RECALLODEFCN - This function allows to compute some parameters used
inside the ODE integrations

INPUTS:     
            - fun, ode function used;
            - T, integration time vector;
            - Y, state matrix.

OUTPUTS:
            - all_steps, structure which contains all the parameters needed. 

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept 
email: adriano.filippo.inno@skywarder.eu
Release date: 16/11/2018

%}

NT = length(T);
fun_info = functions(fun);

for i = 1:NT
    
    [~,single_step] = fun(T(i),Y(i,:),varargin{:});
    
    
    all_steps.velocities(i, 1:3) = single_step.velocities;
    
    all_steps.accelerations.body_acc(i, 1:3) = single_step.accelerations.body_acc;
    
    all_steps.forces.AeroDyn_Forces(i, 1:3) = single_step.forces.AeroDyn_Forces;
    
    all_steps.coeff.CA(i, 1) = single_step.coeff.CA;
    
end