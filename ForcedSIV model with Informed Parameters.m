clear all;
close all;
clc;
options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3]);
figure
[t, x] = ode45('siv', [0 180], [8000 1000 1000], options);
plot(t, x(:,1), 'b.', t, x(:,2), 'r.', t, x(:,3), 'g.')
legend('S_1', 'I_1', 'V_1')
hold on
[t, x] = ode45('siv', [0 180], [1000 2000 7000], options);
plot(t, x(:,1), 'b:', t, x(:,2), 'r:', t, x(:,3), 'g:')
hold on
[t, x] = ode45('siv', [0 180], [5500 500 4000], options);
plot(t, x(:,1), 'b-', t, x(:,2), 'r-', t, x(:,3), 'g-')
xlabel('Months')
ylabel('Individuals Infected')
legend('S_1', 'I_1', 'V_1', 'S_2', 'I_2', 'V_2', 'S_3', 'I_3', 'V_3', 'Location','best')
title('Forced SIV Model with Informed Parameters')
