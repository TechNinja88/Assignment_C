%{
 Let us understand the implementation of PSO using a simple
 function (Sphere) which calculates the sum of squares. Our objective is to
 minimize the sum of squares and find out the value of decision variables.
 However, we know that the minimum (global minima) of this function will be
 when all decision variables' values will be 0, but here
 the objective is to understand the PSO. Later on, you can use PSO
 for complex problems.
 References:
1.Eberhart, R., & Kennedy, J. (1995, November). 
Particle swarm optimization. In Proceedings of the
IEEE international conference on neural networks (Vol. 4, pp. 1942-1948).
2. Yarpiz (2023). Video Tutorial of Particle Swarm Optimization (PSO) 
%}
%% Clear screen, workspace and close figure if opened
clc;
clear;
close all;

taxi_locations = [
    10.7522, 59.9139; % Oslo
    5.3221, 60.3913;  % Bergen
    10.3951, 63.4305; % Trondheim
    5.7331, 58.9690;  % Stavanger
    6.1495, 62.4722;  % Ålesund
    18.9560, 69.6496; % Tromsø
    8.4689, 61.1145;  % Kristiansand
    9.6489, 59.9494;  % Drammen
    11.1000, 60.7945; % Sandnes
    7.1000, 58.1599   % Haugesund
];

% Customer locations
customer_locations = [
    10.9000, 59.9900; % Customer 1
    5.5000, 60.5000; % Customer 2
    10.4900, 63.4900; % Customer 3
    5.7900, 58.9900;  % Customer 4
    6.1900, 62.4900;  % Customer 5
    18.9900, 69.6900; % Customer 6
    8.5900, 61.1190;  % Customer 7
    9.6900, 59.9900;  % Customer 8
    11.0000, 60.8900; % Customer 9
    7.0000, 58.1900   % Customer 10
];

nTaxis = size(taxi_locations, 1); % Number of taxis
nCustomers = size(customer_locations, 1); % Number of customers

% Calculate Euclidean distance matrix
distance_matrix = zeros(nTaxis, nCustomers);
for i = 1:nTaxis
    for j = 1:nCustomers
        distance_matrix(i, j) = sqrt((taxi_locations(i, 1) - customer_locations(j, 1))^2 + ...
                                (taxi_locations(i, 2) - customer_locations(j, 2))^2);
    end
end

% Define parameters for PSO
num_particles = 30;
max_iter = 100;
w = 0.7; % Inertia weight
c1 = 1.5; % Cognitive (personal) weight
c2 = 1.5; % Social (global) weight

% Initialize particles (binary matrices)
particles = zeros(num_particles, nTaxis, nCustomers); % Each particle is a binary matrix
velocities = rand(num_particles, nTaxis, nCustomers) - 0.5; % Initialize velocities

% Initialize personal best positions and global best
pbest = particles;
pbest_scores = inf(num_particles, 1);
gbest = zeros(nTaxis, nCustomers); % Global best position
gbest_score = inf; % Global best score

% Fitness function
fitness_function = @(x) sum(sum(x .* distance_matrix));

% Main PSO loop
convergence_curve = zeros(1, max_iter);
for iter = 1:max_iter
    for i = 1:num_particles
        % Update velocity
        r1 = rand(nTaxis, nCustomers);
        r2 = rand(nTaxis, nCustomers);
        velocities(i, :, :) = w * squeeze(velocities(i, :, :)) + ...
                              c1 * r1 .* (squeeze(pbest(i, :, :)) - squeeze(particles(i, :, :))) + ...
                              c2 * r2 .* (gbest - squeeze(particles(i, :, :)));

        % Update position (convert velocity to probabilities)
        probabilities = 1 ./ (1 + exp(-squeeze(velocities(i, :, :))));
        particles(i, :, :) = rand(nTaxis, nCustomers) < probabilities;

        % Ensure one-to-one assignment
        particles(i, :, :) = EnsureOneToOne(squeeze(particles(i, :, :)));

        % Evaluate fitness
        fitness = fitness_function(squeeze(particles(i, :, :)));
        if fitness < pbest_scores(i)
            pbest_scores(i) = fitness;
            pbest(i, :, :) = particles(i, :, :);
        end

        % Update global best
        if fitness < gbest_score
            gbest_score = fitness;
            gbest = squeeze(particles(i, :, :));
        end
    end

    % Store convergence data
    convergence_curve(iter) = gbest_score;

    % Display iteration info
    fprintf('Iteration %d, Best Fitness: %f\n', iter, gbest_score);
end

% Plot convergence curve with improvements
figure;
plot(1:max_iter, convergence_curve, 'LineWidth', 2, 'Color', 'b');
xlabel('Iteration');
ylabel('Best Fitness (Total Distance)');
title('Convergence of Binary PSO');
subtitle(sprintf('Minimum Distance: %.2f km', gbest_score));
grid on;

% Visualize results with improved plot
figure;
hold on;

% Plot taxis
plot(taxi_locations(:, 1), taxi_locations(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Taxis');

% Plot customers
plot(customer_locations(:, 1), customer_locations(:, 2), 'bx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Customers');

% Add legend, title, and labels
legend('Location', 'best');
title('Taxi-Customer Assignment using Binary PSO');
subtitle(sprintf('Total Distance: %.2f km', gbest_score));
xlabel('Longitude');
ylabel('Latitude');
grid on;
axis equal; % Ensure equal scaling for x and y axes
hold off;

% Function to ensure one-to-one assignment
function particle = EnsureOneToOne(particle)
    [nTaxis, nCustomers] = size(particle);
    % Ensure each taxi is assigned to exactly one customer
    for i = 1:nTaxis
        if sum(particle(i, :)) ~= 1
            [~, idx] = max(particle(i, :)); % Assign to the most probable customer
            particle(i, :) = 0;
            particle(i, idx) = 1;
        end
    end
    % Ensure each customer is assigned to exactly one taxi
    for j = 1:nCustomers
        if sum(particle(:, j)) ~= 1
            [~, idx] = max(particle(:, j)); % Assign to the most probable taxi
            particle(:, j) = 0;
            particle(idx, j) = 1;
        end
    end
end

initial_cost = convergence_curve(1);
minimum_cost = gbest_score;

fprintf('Initial Cost: %f\n', initial_cost);
fprintf('Minimum Cost: %f\n', minimum_cost);
disp('Best Allocation (Global Best Position):');
disp(gbest);