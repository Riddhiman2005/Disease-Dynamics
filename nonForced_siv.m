% This function creates a non-forced SIV system of ODE
function dx = nonForced_siv(t, x)
    dx = [0; 0; 0];
    alpha = 0.0035;
    beta = 0.56;
    gamma = 0.0056;
    kappa = 0.0028;
    
    dx(1) = -alpha * x(1) * x(2) + beta * x(2) - gamma * x(1) + kappa * x(3);
    dx(2) = alpha * x(1) * x(2) - beta * x(2) - gamma * x(2);
    dx(3) = gamma * x(1) + gamma * x(2) - kappa * x(3);
end
