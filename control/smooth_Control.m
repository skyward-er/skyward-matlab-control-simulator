function u = smooth_Control(u, uprev, filter)

% Smooth the control variable with a filter

u = filter*u + (1 - filter)*uprev;

u = round(u);
end