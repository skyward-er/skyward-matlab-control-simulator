% Salva le traiettorie

M = load('Trajectories.mat')

for i=1:11
    name=sprintf('Traiettoria_%d.txt',i);
    t = [M.trajectories_saving(i).Z_ref , M.trajectories_saving(i).V_ref]
    csvwrite(name,t)
end
