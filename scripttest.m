clear
clc
% Example script to run the two_site_dmrg function

% Parameters for the Hubbard model
N = 6;             % Number of sites
bd = 2;             % Local space dimension (spin up and spin down)
max_dims=8;           % Bond dimension
D_periodic = 4;
D = 10;              %Max Bond dimension
U = 0;             % Onsite interaction strength
t = 1;             % Hopping parameter
max_sweeps = 10;   % Maximum number of sweeps
%  tol = 1e-6;        % Convergence tolerance
 NrEl=2;            %Number of elements
mu =0;


% Run the two_site_dmrg function
% [lowest_energy, energy_values, M, E_exact] = two_site_dmrg(N, bd, U, mu, t, NrEl,max_sweeps);
 % [lowest_energy, full_sweep_energy, M, E_exact,occupation_number_old,occupation_number_new] = two_site_dmrg(N, bd, max_dims, U, mu, t, NrEl, max_sweeps);
 [E_iter,new_energy_values, M, E_exact,occupation_number_old,occupation_number_new] = two_site_dmrg(N, bd, U,max_dims,D,D_periodic, mu, t, NrEl, max_sweeps);

% % Initialize a matrix to hold the full sweep energies for all D values
% energy_per_site = zeros(length(max_dims), max_sweeps);
% groundstate_errors = zeros(length(max_dims), max_sweeps);
% 
% % Iterate over the maximum bond dimensions
% for i = 1:length(max_dims)
%     % Run the two_site_dmrg function with the current D and fixed D_trunc
%    [energy_values_temp, M, E_exact,occupation_number_old,occupation_number_new] = two_site_dmrg(N, bd, U,max_dims,D,D_periodic, mu, t, NrEl, max_sweeps);
% 
% 
%       % Calculate and store the energy per site values
%     energy_per_site(i,:) = energy_values_temp;  % N is the number of sites
% 
%  % Assuming energy_values_temp correctly represents energy values for each sweep
%     for sweep = 1:max_sweeps
%         groundstate_errors(i, sweep) = (energy_per_site(i,sweep)-E_exact)/abs(E_exact); 
%     end
% 
% end
% 
% % Plot for relative error
% figure; % Create a new figure for relative error
% for i = 1:length(max_dims)
%     loglog(1:numel(energy_per_site(i,:)), groundstate_errors(i,:), 'LineWidth', 1, 'DisplayName', ['D = ' num2str(max_dims(i))]);
%     hold on; % Allows overlaying of multiple lines
% end
% hold off; % End adding lines to this figure
% set(gca, 'FontSize', 13, 'LineWidth', 1); % Set the axis properties
% xlim([0 max_sweeps]);
% xlabel('# of sweeps', 'FontSize', 13); % X-axis label
% ylabel('Relative error', 'FontSize', 13); % Y-axis label
% title('Relative Error vs. Number of Sweeps', 'FontSize', 14); % Title
% legend('show', 'Location', 'northeastoutside'); % Show legend
% grid on; % Add grid
% % set(gca, 'YDir', 'reverse'); % Reverse the direction of the y-axis
%    % axis([0 5 -0.19 -0.18])
% 
% 
% 
% 
% 
% 
% 
% 
% 






























% clear
% clc
% 
% % Common parameters
% N = 6;
% bd = 2;
% t = 1;
% NrEl = 2;
% max_sweeps = 20;
% mu = 0;

% % Varying U
% U_values = [0.5, 1, 2, 4];
% lowest_energies = zeros(1, length(U_values)); % Array to store lowest energies
% E_exacts = zeros(1, length(U_values)); % Array to store exact energies
% 
% % Run simulations
% for i = 1:length(U_values)
%     U = U_values(i);
%     [lowest_energy, full_sweep_energy, M, E_exact] = two_site_dmrg(N, bd, U, mu, t, NrEl, max_sweeps);
%     lowest_energies(i) = lowest_energy; % Store lowest energy
%     E_exacts(i) = E_exact; % Store exact energy
% end
% % Plotting
% figure;
% plot(U_values, lowest_energies, '-o', 'DisplayName', 'DMRG Lowest Energy');
% hold on;
% plot(U_values, E_exacts, '-x', 'DisplayName', 'Exact Energy');
% xlabel('Onsite Interaction Strength (U)');
% ylabel('Energy');
% title('Energy vs Onsite Interaction Strength');
% legend;
% grid on;


% % Varying t
% t_values = [0.5, 1, 1.5, 2];
% lowest_energies = zeros(1, length(t_values)); % Array to store lowest energies
% E_exacts = zeros(1, length(t_values)); % Array to store exact energies
% 
% % Run simulations
% for i = 1:length(t_values)
%     U = t_values(i);
%     [lowest_energy, full_sweep_energy, M, E_exact] = two_site_dmrg(N, bd, U, mu, t, NrEl, max_sweeps);
%     lowest_energies(i) = lowest_energy; % Store lowest energy
%     E_exacts(i) = E_exact; % Store exact energy
% end
% 
% 
% 
% % Plotting
% figure;
% plot(t_values, lowest_energies, '-o', 'DisplayName', 'DMRG Lowest Energy');
% hold on;
% plot(t_values, E_exacts, '-x', 'DisplayName', 'Exact Energy');
% xlabel('hopping parameter (t)');
% ylabel('Energy');
% title('Energy vs Hopping paramter');
% legend;
% grid on;
% 





























