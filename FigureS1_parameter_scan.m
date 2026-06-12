%% ================================================================
% Gray-Scott Reaction-Diffusion Model
% 2D periodic domain
% Fourth-order Runge-Kutta (RK4)
%
% Suitable for large parameter scans
%
% Author: Xinyu Wang
% Date: 2026
%% ================================================================

clear;
clc;
close all;

%% ================================================================
% For Reproducibility
% rng(1,'twister');

%% ================================================================
% Output folders
%% ================================================================

result_folder = 'results';
figure_folder = fullfile(result_folder,'figures');
data_folder   = fullfile(result_folder,'data');

if ~exist(result_folder,'dir')
    mkdir(result_folder);
end

if ~exist(figure_folder,'dir')
    mkdir(figure_folder);
end

if ~exist(data_folder,'dir')
    mkdir(data_folder);
end

%% ================================================================
% Numerical parameters
%% ================================================================

Nx = 200;
Ny = 200;

Lx = Nx;
Ly = Ny;

dx = Lx/Nx;
dy = Ly/Ny;

dt = 0.04;
steps = 400000;
plot_interval = 1000;

%% ================================================================
% Convergence criteria
%% ================================================================

convergence_tol = 1e-4;
uniform_tol = 1e-10;

%% ================================================================
% Parameter scan
%% ================================================================

Du_list = 0.1:0.1:1.4;
Dv_list = 0.05:0.05:1.0;

Feed_list = 0.045;
kill_list = 0.06;

%% ================================================================
% Storage preallocation
%% ================================================================

max_cases = length(Du_list) * ...
            length(Dv_list) * ...
            length(Feed_list) * ...
            length(kill_list);

record_u = zeros(Nx,Ny,max_cases,'single');
record_v = zeros(Nx,Ny,max_cases,'single');

record_parameters = zeros(max_cases,4);

record_metadata = strings(max_cases,1);

count = 0;

%% ================================================================
% Spatial grid
%% ================================================================

[x,y] = meshgrid(1:Nx,1:Ny);

%% ================================================================
% Laplacian operator (periodic boundary)
%% ================================================================

lap = @(M) ( ...
    circshift(M,[0 1])  + ...
    circshift(M,[0 -1]) + ...
    circshift(M,[1 0])  + ...
    circshift(M,[-1 0]) - ...
    4*M ) / dx^2;

%% ================================================================
% Main parameter scan
%% ================================================================

for Du_gs = Du_list

    for Dv_gs = Dv_list

        for Feed = Feed_list

            for kill = kill_list

                fprintf('\n=======================================\n');
                fprintf('Running simulation:\n');
                fprintf('Du = %.3f\n',Du_gs);
                fprintf('Dv = %.3f\n',Dv_gs);
                fprintf('Feed = %.3f\n',Feed);
                fprintf('Kill = %.3f\n',kill);
                fprintf('=======================================\n');

                %% ==========================================
                % Initial condition
                %% ==========================================

                u_gs = ones(Ny,Nx);
                v_gs = zeros(Ny,Nx);

                % Small random perturbation
                u_gs = u_gs + 0.01*randn(Ny,Nx);
                v_gs = v_gs + 0.01*randn(Ny,Nx);

                % Gaussian seeds
                sigma = 3;
                amplitude = 0.30;
                number_of_seeds = 20;

                for seed_id = 1:number_of_seeds

                    cx = randi([1 Nx]);
                    cy = randi([1 Ny]);

                    gaussian_seed = exp( ...
                        -((x-cx).^2 + (y-cy).^2) ...
                        /(2*sigma^2));

                    v_gs = v_gs + gaussian_seed;

                end

                v_gs = amplitude * ...
                    v_gs / max(v_gs(:));

                %% ==========================================
                % Define RHS operators
                %% ==========================================

                rhs_u = @(u,v) ...
                    Du_gs*lap(u) ...
                    - u.*v.^2 ...
                    + Feed*(1-u);

                rhs_v = @(u,v) ...
                    Dv_gs*lap(v) ...
                    + u.*v.^2 ...
                    - (Feed+kill)*v;

                %% ==========================================
                % Visualization setup
                %% ==========================================

                fig = figure('Position',[100 100 1200 500]);

                subplot(1,2,1);
                h1 = imagesc(u_gs);
                axis equal off;
                title('u field');
                colorbar;

                subplot(1,2,2);
                h2 = imagesc(v_gs);
                axis equal off;
                title('v field');
                colorbar;

                colormap(turbo);

                %% ==========================================
                % Convergence monitoring
                %% ==========================================

                Ulast = u_gs;

                simulation_status = "unfinished";

                %% ==========================================
                % Main time evolution
                %% ==========================================

                for n = 1:steps

                    %% ======================================
                    % RK4 integration
                    %% ======================================

                    % ---------- k1 ----------

                    k1u = rhs_u(u_gs,v_gs);
                    k1v = rhs_v(u_gs,v_gs);

                    % ---------- k2 ----------

                    u2 = u_gs + 0.5*dt*k1u;
                    v2 = v_gs + 0.5*dt*k1v;

                    k2u = rhs_u(u2,v2);
                    k2v = rhs_v(u2,v2);

                    % ---------- k3 ----------

                    u3 = u_gs + 0.5*dt*k2u;
                    v3 = v_gs + 0.5*dt*k2v;

                    k3u = rhs_u(u3,v3);
                    k3v = rhs_v(u3,v3);

                    % ---------- k4 ----------

                    u4 = u_gs + dt*k3u;
                    v4 = v_gs + dt*k3v;

                    k4u = rhs_u(u4,v4);
                    k4v = rhs_v(u4,v4);

                    % ---------- Final update ----------

                    u_gs = u_gs + ...
                        (dt/6)*( ...
                        k1u + 2*k2u + 2*k3u + k4u);

                    v_gs = v_gs + ...
                        (dt/6)*( ...
                        k1v + 2*k2v + 2*k3v + k4v);

                    %% ======================================
                    % Positivity constraint
                    %% ======================================

                    u_gs = max(u_gs,0);
                    v_gs = max(v_gs,0);

                    %% ======================================
                    % Visualization and diagnostics
                    %% ======================================

                    if mod(n,plot_interval)==0

                        set(h1,'CData',u_gs);
                        set(h2,'CData',v_gs);

                        sgtitle(sprintf( ...
                            'step=%d | F=%.3f | k=%.3f | Du=%.2f | Dv=%.2f',...
                            n,Feed,kill,Du_gs,Dv_gs));

                        drawnow limitrate;

                        %% ==================================
                        % Relative convergence error
                        %% ==================================

                        relative_error = ...
                            norm(u_gs(:)-Ulast(:)) ...
                            / max(norm(u_gs(:)),1e-12);

                        fprintf('step = %d | error = %.3e\n',...
                            n,relative_error);

                        %% ==================================
                        % Convergence detection
                        %% ==================================

                        if relative_error < convergence_tol

                            fprintf('Pattern converged.\n');

                            simulation_status = "converged";

                            break;

                        end

                        %% ==================================
                        % Uniform state detection
                        %% ==================================

                        if var(v_gs(:)) < uniform_tol

                            fprintf('Uniform steady state detected.\n');

                            simulation_status = "uniform_state";

                            break;

                        end

                        %% ==================================
                        % No-pattern detection
                        %% ==================================

                        % if sum(v_gs(:) > 0.2) < 0.01*numel(v_gs)
                        % 
                        %     fprintf('No pattern detected.\n');
                        % 
                        %     simulation_status = "no_pattern";
                        % 
                        %     break;
                        % 
                        % end

                        Ulast = u_gs;

                    end

                end

                %% ==========================================
                % Save simulation result
                %% ==========================================

                count = count + 1;

                record_u(:,:,count) = single(u_gs);
                record_v(:,:,count) = single(v_gs);

                record_parameters(count,:) = ...
                    [Feed kill Du_gs Dv_gs];

                metadata_text = sprintf(...
                    'status=%s | steps=%d | dt=%.4f',...
                    simulation_status,n,dt);

                record_metadata(count) = metadata_text;

                %% ==========================================
                % Save figure
                %% ==========================================

                figure_name = sprintf(...
                    'pattern_%03d_F%.3f_k%.3f_Du%.2f_Dv%.2f.png',...
                    count,Feed,kill,Du_gs,Dv_gs);

                saveas(fig,...
                    fullfile(figure_folder,figure_name));

                close(fig);

                %% ==========================================
                % Intermediate save
                %% ==========================================

                save(fullfile(data_folder,'gray_scott_results.mat'),...
                    'record_u',...
                    'record_v',...
                    'record_parameters',...
                    'record_metadata',...
                    '-v7.3');

            end

        end

    end

end

%% ================================================================
% Final summary
%% ================================================================

fprintf('\n=======================================\n');
fprintf('All simulations completed.\n');
fprintf('Total saved cases = %d\n',count);
fprintf('=======================================\n');
