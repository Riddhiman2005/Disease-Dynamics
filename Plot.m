% Code to plot Influenza infections by year from 2010 to 2024

years = 2010:2024;
infected = [20604, 603, 5044, 5253, 937, 42592, 1786, 38811, 15266, 28798, 2752, 778, 13202, 8125, 4871];

figure;
plot(years, infected, '-o', 'LineWidth', 2);
title('Influenza Infections by Year from 2010 to 2024');
xlabel('Year');
ylabel('Number of Infected Individuals');

xticks(years);
