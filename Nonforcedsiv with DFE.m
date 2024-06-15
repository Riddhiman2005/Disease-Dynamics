>> clear all;
close all;
clc;
options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3]);
figure
[t, x] = ode45('nonForced_siv', [0 240], [1000 0 0], options);
plot(t, x(:,1), 'b:', t, x(:, 3), 'g:')
hold on
[t, x] = ode45('nonForced_siv', [0 240], [800 0 200], options);
plot(t, x(:,1), 'b.', t, x(:, 3), 'g.')
hold on
[t, x] = ode45('nonForced_siv', [0 240], [600 0 400], options);
plot(t, x(:,1), 'b-', t, x(:, 3), 'g-')
xlabel('Months')
ylabel('Individuals Infected')
legend('S', 'V', 'S_2', 'V_2', 'S_3', 'V_3', 'Location','best')
title('Non-forced SIV Model with DFE Initial Conditions')
