>> clear all;
close all;
clc;
options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3]);
figure
[t, x] = ode45('nonForced_siv', [0 240], [400 100 500], options);
plot(t, x(:,1), 'b.', t, x(:,2), 'r.', t, x(:,3), 'g.')
legend('S_1', 'I_1', 'V_1')
hold on
[t, x] = ode45('nonForced_siv', [0 240], [500 200 300], options);
plot(t, x(:,1), 'b:', t, x(:,2), 'r:', t, x(:,3), 'g:')
hold on
[t, x] = ode45('nonForced_siv', [0 240], [100 500 400], options);
plot(t, x(:,1), 'b-', t, x(:,2), 'r-', t, x(:,3), 'g-')
xlabel('Months')
ylabel('Individuals Infected')
legend('S_1', 'I_1', 'V_1', 'S_2', 'I_2', 'V_2', 'S_3', 'I_3', 'V_3', 'Location','best')
title('Non-forced SIV model')
