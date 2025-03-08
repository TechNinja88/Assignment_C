
clc;
clear;

%% Problem setting
x = [10.7522, 5.3221, 10.3951, 5.7331, 6.1495, 18.9560, 8.4689, 9.6489, 11.1000, 7.1000];
y = [59.9139, 60.3913, 63.4305, 58.9690, 62.4722, 69.6496, 61.1145, 59.9494, 60.7945, 58.1599];
names = {'Oslo', 'Bergen', 'Trondheim', 'Stavanger', 'Ålesund', 'Tromsø', 'Kristiansand', 'Drammen', 'Sandnes', 'Haugesund'};
nCities = numel(x);

lb = ones(1, nCities); % Lower bound
ub = nCities * ones(1, nCities); % Upper bound
prob = @(route) TSPCost(route, x, y); % Fitness function

%% Algorithm parameters
Np = 100; % Population size
T = 200; % Number of iterations
limit = 20; % Permissible number of failures

%% Starting of ABC
f = NaN(Np, 1); % Vector to store the objective function value of the population
fit = NaN(Np, 1); % Vector to store the fitness function value of the population
trial = zeros(Np, 1); % Trial vector
D = nCities; % Number of decision variables in the problem
P = zeros(Np, D); % Population

% Generate the initial population
for i = 1:Np
    P(i, :) = randperm(D);
end

BestCost = zeros(T, 1);

for p = 1:Np
    f(p) = prob(P(p, :)); % Objective function value
    fit(p) = CalFit(f(p)); % Fitness function value
end

[bestobj, ind] = min(f); % Determine and memorize the best objective value
bestsol = P(ind, :); % Determine and memorize the best solution

for t = 1:T
    % Employed Bee Phase
    for i = 1:Np
        [trial, P, fit, f] = GenNewSol(prob, lb, ub, Np, i, P, fit, trial, f, D);
    end

    % Onlooker Bee Phase
    probability = 0.9 * (fit / max(fit)) + 0.1;
    m = 0;
    n = 1;
    while m < Np
        if rand < probability(n)
            [trial, P, fit, f] = GenNewSol(prob, lb, ub, Np, n, P, fit, trial, f, D);
            m = m + 1;
        end
        n = mod(n, Np) + 1;
    end

    [bestobj, ind] = min([f; bestobj]);
    CombinedSol = [P; bestsol];
    bestsol = CombinedSol(ind, :);

    % Scout Bee Phase
    [val, ind] = max(trial);
    if val > limit
        trial(ind) = 0; % Reset trial value
        P(ind, :) = randperm(D); % Generate a random solution
        f(ind) = prob(P(ind, :)); % Determine the objective function value of new solution
        fit(ind) = CalFit(f(ind)); % Determine the fitness value of new solution
    end

    disp(['Iteration ' num2str(t) ': Best solution = ' num2str(bestsol) ' Best cost ' num2str(bestobj)]);
    BestCost(t) = bestobj;
end

%% Results
figure;
plot(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

figure;
hold on;
plot([x(bestsol) x(bestsol(1))], [y(bestsol) y(bestsol(1))], 'o-', 'LineWidth', 2);
text(x, y, names, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('Longitude');
ylabel('Latitude');
title('TSP Route');
grid on;
hold off;

%% Supporting Functions
function cost = TSPCost(route, x, y)
    cost = 0;
    for i = 1:length(route) - 1
        cost = cost + sqrt((x(route(i)) - x(route(i+1)))^2 + (y(route(i)) - y(route(i+1)))^2);
    end
    cost = cost + sqrt((x(route(end)) - x(route(1)))^2 + (y(route(end)) - y(route(1)))^2);
end

function fitness = CalFit(f)
    if f >= 0
        fitness = 1 / (1 + f);
    else
        fitness = 1 + abs(f);
    end
end

function [trial, P, fit, f] = GenNewSol(prob, lb, ub, Np, i, P, fit, trial, f, D)
    phi = -1 + 2 * rand(1, D);
    k = i;
    while k == i
        k = randi([1 Np]);
    end
    newSol = P(i, :) + phi .* (P(i, :) - P(k, :));
    newSol = max(lb, min(ub, newSol)); % Ensure new solution is within bounds

    newSol = round(newSol); % Round to nearest integer to keep it as a permutation
    newSol = unique(newSol, 'stable'); % Ensure it remains a permutation

    if length(newSol) < D
        missing = setdiff(1:D, newSol);
        newSol = [newSol missing];
    end

    newf = prob(newSol);
    newfit = CalFit(newf);

    if newfit > fit(i)
        P(i, :) = newSol;
        f(i) = newf;
        fit(i) = newfit;
        trial(i) = 0;
    else
        trial(i) = trial(i) + 1;
    end
end