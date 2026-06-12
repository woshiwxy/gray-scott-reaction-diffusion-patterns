%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure1_simulation.m
%
% Generate Figure 1c–e in the manuscript.
%
% Gray–Scott reaction–diffusion model.
%
% Fixed parameters:
%   F  = 0.047
%   k  = 0.060
%   Du = 1.20
%
% Variable parameter:
%   Dv = [0.60, 0.36, 0.22]
%
% Output:
%   Three simulated patterns corresponding to Figure 1c–e.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc; close all;

% -------------------------------
% 参数与网格
% -------------------------------
Nx = 200; Ny = 200;
dx = 1.0; dy = 1.0;
dt = 0.02;
steps = 100000;
plot_interval = 1000;

[x, y] = meshgrid(1:Nx, 1:Ny);
record_parameters=zeros(1,4);
count=0;
record_u=zeros(Nx,Ny,1);
record_v=zeros(Nx,Ny,1);

% test=1:10000

for Feed = 0.047
    for k = 0.06:0.01:0.06
        for Du_gs = 1.2
            for Dv_gs = [0.60, 0.36, 0.22]

                u_gs = ones(Ny, Nx);
                v_gs = zeros(Ny, Nx);
                u_gs = u_gs - 0.02*(rand(Ny, Nx));
                v_gs = v_gs + 0.02*(rand(Ny, Nx));


                rng('shuffle');
                sigma=5;
                amp=0.5;
                for i = 1:50
                    cx = randi([1 Nx]);
                    cy = randi([1 Ny]);
                    g = exp(-((x-cx).^2 + (y-cy).^2)/(2*sigma^2));
                    v_gs = v_gs + g;
                end
                v_gs = amp * v_gs / max(abs(v_gs(:)));

                v_gs = v_gs + 0.15*cos(2*pi*x/35+2*pi*y/35);

                Ulast=u_gs;
                Flast=0;

                % -------------------------------
                % 4-point Laplace
                % -------------------------------
                % lap = @(M) (circshift(M,[0 1]) + circshift(M,[0 -1]) + ...
                %     circshift(M,[1 0]) + circshift(M,[-1 0]) - 4*M)/(dx^2);

                % -------------------------------
                % 9-point Laplacian (periodic BC)
                % -------------------------------
                lap = @(M) ( 4*( circshift(M,[0 1])  + circshift(M,[0 -1]) + ...
                    circshift(M,[1 0])  + circshift(M,[-1 0]) ) + ( ...
                    circshift(M,[1 1])   + circshift(M,[1 -1]) + ...
                    circshift(M,[-1 1])  + circshift(M,[-1 -1]) ...
                    ) - 20*M ) / (6*dx^2);


                figure('Position',[100,100,1200,500]);

                % -------------------------------
                % main
                % -------------------------------

                for n = 1:steps
                    % ---------- Gray–Scott ----------
                    Lu = lap(u_gs);
                    Lv = lap(v_gs);
                    du = Du_gs*Lu - u_gs.*v_gs.^2 + Feed*(1 - u_gs); % +0.001*randn(size(u_gs));
                    dv = Dv_gs*Lv + u_gs.*v_gs.^2 - (Feed + k)*v_gs; % +0.001*randn(size(u_gs));
                    u_gs = u_gs + dt*du;
                    v_gs = v_gs + dt*dv;

                    if mod(n, plot_interval) == 0
                        subplot(1,2,1);
                        imagesc(u_gs); axis equal off;
                        title(['Gray–Scott (u), step = ', num2str(n)]);
                        colormap(gray); colorbar;

                        subplot(1,2,2);
                        imagesc(v_gs); axis equal off;
                        title(['Gray–Scott (v), step = ', num2str(n)]);
                        colormap(gray); colorbar;
                        sgtitle(['parameters: Feed=',num2str(Feed),' k=',num2str(k),' Du=',num2str(Du_gs),' Dv=',num2str(Dv_gs),''])
                        
                        drawnow;


                        if mod(n,plot_interval)==0
                            U=u_gs;
                            V=v_gs;
                            EU = norm(U(:)-Ulast(:)) / norm(U(:));
                            F = abs(fft2(U));
                            EF = norm(F(:)-Flast(:))/norm(F(:));

                            EpsU = 1e-3;
                            EpsF = 2e-3;

                            if EU < EpsU && EF < EpsF
                                disp(['Converged at n = ' num2str(n)])
                                break
                            end

                            if sum(U(:) > 0.2) < numel(U)*0.01 || sum(V(:) > 0.01) < numel(V)*0.01
                                disp(['No pattern detected → stop at n=' num2str(n)])
                                close all
                                break
                            end
                            if (var(u_gs(:)) < 1e-8 && var(v_gs(:)) < 1e-8) % && mean(v_gs(:)) > 0.01
                                break
                            end
                            if n>=steps
                                break
                            end

                            Ulast = U;
                            Flast = F;
                        end
                    end
                end
                filename = sprintf( ...
                    'pattern_Du%.2f_Dv%.2f_F%.3f_k%.3f.png', ...
                    Du_gs,Dv_gs,Feed,k);
                print(gcf,filename,'-dpng','-r300');
            end
        end
    end
end

