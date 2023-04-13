% Define parameters and observations
gamma = 1;
mu = 5;
sigma = 0.016;
rho = 0.95;
g_x = 1.0029;
i_y = 0.2133;
n = 1/3;
k_y = 10;
a = 0.36;
delta = 0.06;
z = 1;

% Set a range for rho and sigma_z
rho_vals = linspace(0.9, 1, 20);
sigma_z_vals = linspace(0.01, 0.05, 20);

% Initialize matrices for results
fit_results = zeros(length(rho_vals), length(sigma_z_vals));

% Loop over possible values of rho and sigma_z
for i = 1:length(rho_vals)
    for j = 1:length(sigma_z_vals)
        rho_z = rho_vals(i);
        sigma_z = sigma_z_vals(j);
        
        % Find steady-state values using fsolve function
        steady_state = fsolve(@(x) model_steady_state(x, gamma, mu, i_y, n, k_y, sigma, rho, g_x, rho_z, sigma_z, a, delta, z), [1, 1, 1]);
        
        % Calculate fit using an appropriate metric (e.g., mean squared error)
        fit = fit_metric(steady_state, i_y, n, k_y, a);
        
        % Store fit result in the matrix
        fit_results(i, j) = fit;
    end
end

% Find the best fit and corresponding rho and sigma_z values
[min_fit, idx] = min(fit_results(:));
[rho_idx, sigma_z_idx] = ind2sub(size(fit_results), idx);
best_rho = rho_vals(rho_idx);
best_sigma_z = sigma_z_vals(sigma_z_idx);

% Display results
fprintf('Best fit: %f\n', min_fit);
fprintf('Best rho: %f\n', best_rho);
fprintf('Best sigma_z: %f\n', best_sigma_z);

% Function to solve for steady-state values
function res = model_steady_state(x, gamma, mu, i_y, n, k_y, sigma, rho, g_x, rho_z, sigma_z, a, delta, z)
    c_hat = x(1);
    k_hat = x(2);
    theta = x(3);
    
    % Calculate i_hat from the given parameters
    i_hat = i_y * (theta * k_hat^a * n^(1-a));
    
    % FOC equations
    res(1) = (1-n)^(-mu) - c_hat^(-gamma) * theta * k_hat^a * (1-a) * n^(-a);
    res(2) = z * g_x - theta * a * k_hat^(a-1) * n^(1-a) + (1-delta);
    res(3) = c_hat + i_hat - theta * k_hat^a * n^(1-a);
end

% Function to calculate the fit using mean squared error
function fit = fit_metric(steady_state, i_y, n, k_y, a)
    c_hat = steady_state(1);
    k_hat = steady_state(2);
    theta = steady_state(3);
    
    % Calculate the steady-state values of interest
    i_hat = i_y * (theta * k_hat^a * n^(1-a));
    y_hat = theta * k_hat^a * n^(1-a);
    
    % Calculate the observed ratios
    obs_i_y = i_y;
    obs_n = n;
    obs_k_y = k_y;
    
    % Calculate the model predictions
    pred_i_y = i_hat / y_hat;
    pred_n = n;
    pred_k_y = k_hat / y_hat;
    
    % Calculate the mean squared error
    mse = (obs_i_y - pred_i_y)^2 + (obs_n - pred_n)^2 + (obs_k_y - pred_k_y)^2;
    
    % Set the fit to the mean squared error
    fit = mse;
end