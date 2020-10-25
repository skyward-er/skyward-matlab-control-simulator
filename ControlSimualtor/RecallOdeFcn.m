function [all_steps] = RecallOdeFcn(fun, T, Y, settings, C, varargin)
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
% fun_info = functions(fun);

for i = 1:NT
    
    [~,single_step] = fun(T(i),Y(i,:), settings, C(i), varargin{:});
    
    all_steps.integration.t(i) = single_step.integration.t;
    all_steps.interp.alt(i) = single_step.interp.alt;
    all_steps.wind.body_wind(1:3,i) = single_step.wind.body_wind;
    all_steps.wind.NED_wind(1:3,i) = single_step.wind.NED_wind;
    all_steps.velocities(i, 1:3) = single_step.velocities;
    
    all_steps.air.rho(i) = single_step.air.rho;
    all_steps.air.P(i) = single_step.air.P;
    
    all_steps.accelerations.body_acc(i, 1:3) = single_step.accelerations.body_acc;
    
    all_steps.interp.M(i) = single_step.interp.M;
    all_steps.interp.alpha(i) = single_step.interp.alpha;
    all_steps.interp.beta(i) = single_step.interp.beta;
    
    all_steps.forces.AeroDyn_Forces(i, 1:3) = single_step.forces.AeroDyn_Forces;
    all_steps.forces.T(i) = single_step.forces.T;
    
    all_steps.accelerations.ang_acc(i, 1:3) = single_step.accelerations.ang_acc;
    
    all_steps.coeff.CA(i) = single_step.coeff.CA;
%     all_steps.coeff.CYB(i) = single_step.coeff.CYB;
%     all_steps.coeff.CNA(i) = single_step.coeff.CNA;
%     all_steps.coeff.Cl(i) = single_step.coeff.Cl;
%     all_steps.coeff.Clp(i) = single_step.coeff.Clp;
%     all_steps.coeff.Cma(i) = single_step.coeff.Cma;
%     all_steps.coeff.Cmad(i) = single_step.coeff.Cmad;
%     all_steps.coeff.Cmq(i) = single_step.coeff.Cmq;
%     all_steps.coeff.Cnb(i) = single_step.coeff.Cnb;
%     all_steps.coeff.Cnr(i) = single_step.coeff.Cnr;
%     all_steps.coeff.Cnp(i) = single_step.coeff.Cnp;
    
    if isfield(single_step.coeff, 'XCP')
        all_steps.coeff.XCP(i) = single_step.coeff.XCP;
    end
end