function y = simBasedObjectiveFcn(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Close figures

save('interationData.mat','x')
run( '../start_simulation.m');
load('interationData.mat')

load( '../altitude_velocity');
% figure;
% plot(altitude_velocity.Z_ref,altitude_velocity.V_ref,altitude_velocity.Z_real,altitude_velocity.V_real);

% sort Z_ref
[altitude_velocity.Z_ref,iZ_ref] = sort(altitude_velocity.Z_ref);
altitude_velocity.V_ref = altitude_velocity.V_ref(iZ_ref);
[altitude_velocity.Z_real,iZ_real] = sort(altitude_velocity.Z_real);
altitude_velocity.V_real = altitude_velocity.V_real(iZ_real);

% Check uniqueness
[~,ref_idx] = unique(altitude_velocity.Z_ref); 
Z_ref = altitude_velocity.Z_ref(ref_idx);
V_ref = altitude_velocity.V_ref(ref_idx);

[~,real_idx] = unique(altitude_velocity.Z_real); 
Z_real = altitude_velocity.Z_real(real_idx);
V_real = altitude_velocity.V_real(real_idx);

Zmin=max([min(Z_ref);min(Z_real)]);
Zmax=min([max(Z_ref);max(Z_real)]);
z = linspace(Zmin,Zmax,200);

vref = interp1(Z_ref,V_ref,z);
vreal = interp1(Z_real,V_real,z);

figure(5);
plot(z,vref,z,vreal);

errorPenalty=trapz(z,abs(vref-vreal))

load( '../control_inputs');

figure(6)
plot(control_inputs.time,control_inputs.Angle);

dt=control_inputs.time(2)-control_inputs.time(1);
controlPenalty=trapz(control_inputs.time,abs([0;diff(control_inputs.Angle)]/dt))
% figure
% plot(control_inputs.time,[0;diff(control_inputs.Angle)]/dt):

fprintf('Actual Kp: %g\n', x(1));
fprintf('Actual Ki: %g\n', x(2));
fprintf('Actual errorPenalty: %g\n', errorPenalty);
fprintf('Actual controlPenalty: %g\n', 100*controlPenalty);
fprintf('Actual hole penalty: %g\n', errorPenalty+100*controlPenalty);

y=(errorPenalty+100*controlPenalty);
end

