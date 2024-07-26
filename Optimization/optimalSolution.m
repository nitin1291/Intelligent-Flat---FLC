% Prompt user to select the directory containing the scripts
folderPath = uigetdir('D:\Study Materials\Modelling Uncertainty & Optimization\CW-Modelling\', 'Select Folder');

if folderPath == 0
    error('No folder selected. Please select the correct folder.');
else
    cd(folderPath);
    addpath(folderPath);
end

% Load the shift vectors for the benchmark functions
load('sphere_func_data.mat'); % This file should contain the vector 'o'
load('schwefel_102_data.mat'); % This file should contain the vector 'o'

% Assuming the shift vectors are named 'o' in the files
shift_vector_sphere = o;
shift_vector_schwefel = o;
fixed_bias_sphere = -450; % Fixed bias for the Shifted Sphere function
fixed_bias_schwefel = -450; % Fixed bias for the Shifted Schwefel function

% Define the Shifted Sphere function with fixed bias
shifted_sphere = @(x) sum((x - shift_vector_sphere(1:length(x))).^2) + fixed_bias_sphere;

% Define the Shifted Schwefel’s Problem 1.2 function with fixed bias
shifted_schwefel = @(x) sum(arrayfun(@(i) sum(x(1:i) - shift_vector_schwefel(1:i)).^2, 1:length(x))) + fixed_bias_schwefel;

% Define the Shifted Schwefel’s Problem 1.2 with Noise in Fitness function with fixed bias
shifted_schwefel_noise = @(x) (sum(arrayfun(@(i) sum(x(1:i) - shift_vector_schwefel(1:i)).^2, 1:length(x))) * (1 + 0.4 * abs(normrnd(0, 1)))) + fixed_bias_schwefel;

% Define the optimization techniques
optimizeGA = @(fun, numDimensions) ga(fun, numDimensions, [], [], [], [], [], [], [], ...
    optimoptions('ga', 'MaxGenerations', 100, 'PopulationSize', 50));

optimizePSO = @(fun, numDimensions) particleswarm(fun, numDimensions, [], [], ...
    optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100));

optimizeSA = @(fun, numDimensions) simulannealbnd(fun, rand(1, numDimensions), [], [], ...
    optimoptions('simulannealbnd', 'MaxIterations', 100, 'InitialTemperature', 100));

% Helper function to create the fitness functions
createFitnessFunction = @(fun, dim) @(x) fun(x(1:dim));

% Select functions
functions = {shifted_sphere, shifted_schwefel, shifted_schwefel_noise};
functionNames = {'ShiftedSphere', 'ShiftedSchwefel', 'ShiftedSchwefelNoise'};

% Define dimensions
dimensions = [2, 10];

% Number of runs for statistical analysis
numRuns = 15;

% Store results
results = struct();

for i = 1:length(functions)
    fun = functions{i};
    for j = 1:length(dimensions)
        dim = dimensions(j);

        results.(functionNames{i}).(['D' num2str(dim)]) = struct();

        % Initialize result storage
        gaResults = zeros(numRuns, 1);
        psoResults = zeros(numRuns, 1);
        saResults = zeros(numRuns, 1);

        for k = 1:numRuns
            [~, gaResults(k)] = optimizeGA(createFitnessFunction(fun, dim), dim);
            [~, psoResults(k)] = optimizePSO(createFitnessFunction(fun, dim), dim);
            [~, saResults(k)] = optimizeSA(createFitnessFunction(fun, dim), dim);
        end

        % Calculate statistics
        results.(functionNames{i}).(['D' num2str(dim)]).GA.mean = mean(gaResults);
        results.(functionNames{i}).(['D' num2str(dim)]).GA.std = std(gaResults);
        results.(functionNames{i}).(['D' num2str(dim)]).GA.best = min(gaResults);
        results.(functionNames{i}).(['D' num2str(dim)]).GA.worst = max(gaResults);

        results.(functionNames{i}).(['D' num2str(dim)]).PSO.mean = mean(psoResults);
        results.(functionNames{i}).(['D' num2str(dim)]).PSO.std = std(psoResults);
        results.(functionNames{i}).(['D' num2str(dim)]).PSO.best = min(psoResults);
        results.(functionNames{i}).(['D' num2str(dim)]).PSO.worst = max(psoResults);

        results.(functionNames{i}).(['D' num2str(dim)]).SA.mean = mean(saResults);
        results.(functionNames{i}).(['D' num2str(dim)]).SA.std = std(saResults);
        results.(functionNames{i}).(['D' num2str(dim)]).SA.best = min(saResults);
        results.(functionNames{i}).(['D' num2str(dim)]).SA.worst = max(saResults);
    end
end

% Save results to a file
save('optimization_results789.mat', 'results');

% Load the saved results
load('optimization_results789.mat');

% Function to display results
function displayResults(results, functionNames, dimensions)
    for i = 1:length(functionNames)
        funcName = functionNames{i};
        for j = 1:length(dimensions)
            dim = dimensions(j);
            fprintf('Results for %s Function with D=%d\n', funcName, dim);
            fprintf('Genetic Algorithm - Mean: %f, Std: %f, Best: %f, Worst: %f\n', ...
                results.(funcName).(['D' num2str(dim)]).GA.mean, ...
                results.(funcName).(['D' num2str(dim)]).GA.std, ...
                results.(funcName).(['D' num2str(dim)]).GA.best, ...
                results.(funcName).(['D' num2str(dim)]).GA.worst);
            fprintf('Particle Swarm Optimization - Mean: %f, Std: %f, Best: %f, Worst: %f\n', ...
                results.(funcName).(['D' num2str(dim)]).PSO.mean, ...
                results.(funcName).(['D' num2str(dim)]).PSO.std, ...
                results.(funcName).(['D' num2str(dim)]).PSO.best, ...
                results.(funcName).(['D' num2str(dim)]).PSO.worst);
            fprintf('Simulated Annealing - Mean: %f, Std: %f, Best: %f, Worst: %f\n', ...
                results.(funcName).(['D' num2str(dim)]).SA.mean, ...
                results.(funcName).(['D' num2str(dim)]).SA.std, ...
                results.(funcName).(['D' num2str(dim)]).SA.best, ...
                results.(funcName).(['D' num2str(dim)]).SA.worst);
            fprintf('\n');
        end
    end
end

% Define function names and dimensions
functionNames = {'ShiftedSphere', 'ShiftedSchwefel', 'ShiftedSchwefelNoise'};
dimensions = [2, 10];

% Display the results
displayResults(results, functionNames, dimensions);
