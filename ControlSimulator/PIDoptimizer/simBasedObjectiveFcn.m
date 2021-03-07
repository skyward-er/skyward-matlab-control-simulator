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
figure;
plot(Z_ref,V_ref);
vreal = interp1(Z_real,V_real,z);
% figure;
% plot(z,vref,z,vreal);
% figure;
% plot(z,abs(vref-vreal));
dif=trapz(z,abs(vref-vreal))

y=(dif);
end

