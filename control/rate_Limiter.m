function U = rate_Limiter(U, Uprev, dU_max, dt)

dU = (U - Uprev) / dt;

if (dU > dU_max)
    U =   dt*dU_max + Uprev;
elseif (dU < - dU_max)
    U = - dt*dU_max + Uprev;
end

end