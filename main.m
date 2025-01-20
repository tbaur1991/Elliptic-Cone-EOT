% *************************************************************************
%    I   SSSS   DDDD
%    I   S      D   D 
%    I   SSSS   D   D    Institut f�r Systemdynamik
%    I      S   D   D
%    I   SSSS   DDDD   
% *************************************************************************
% Author: Tim Baur
% University: University of Applied Sciences "HTWG Konstanz"
% E-Mail: tbaur@htwg-konstanz.de
% Creation date: 20.01.2025
% *************************************************************************
% Content: In this simulation environment, the elliptic cylinder 3D EOT
% filter is implemented in a static scenario. A Monte Carlo simulation can
% be started by setting 'mc = true' and 'nMC' as the number of runs to be
% conducted. By setting 'artificial_noise' as true or false, an artificial
% measurement noise can be estimated to process measurements from the
% interior of the object. The variable 'filter' can be set as 'ERHM' for an
% extrusion random hypersurface model or as 'GAM' for a greedy association
% measurement source association model. As measurement 'source'
% approximation a radial or projected association can be chosen. Plots can
% be visualized by setting 'do_plot' as true.
%
% !! The package for the sampling procedure for the S2KF must be downloaded
% from: https://nonlinearestimation.bitbucket.io/ !!
%
% *************************************************************************

close all
clear
clc

% simulation environment
nSamples = 100;
nMeas = 50;
mc = false;
do_plot = true;

% filter settings. filter either 'ERHM' or 'GAM'. source either 'radial' or
% 'projected'.
artificial_noise = true;
filter = 'ERHM';
source = 'radial';
n_upd = 50;

% triangular distributions
pd1 = makedist("Triangular","a",0,"b",0,"c",1);
pd2 = makedist("Triangular","a",0,"b",1,"c",1);

% reference
pos_ref = [2;1;0];
or_ref = 0;
R_ref = [cos(or_ref) -sin(or_ref) 0; sin(or_ref) cos(or_ref) 0; 0 0 1];
a_ref = 3;
b_ref = 1.5;
h_ref = 6;

% memory allocation
nx = 7;
if mc
    nMC = 100;
else
    nMC = 1;
end
X = zeros(nx,nSamples,nMC);
P = zeros(nx,nx,nSamples,nMC);
times = zeros(nSamples,nMC);

% measurement and process noise
sig_r = 0.1;
Q = 1e-4*eye(nx);

% get samples of S2KF
if strcmp(filter,'ERHM')
    nxu = nx + 2*n_upd;
    nxu2 = nx + 2*n_upd/2;
elseif strcmp(filter,'GAM')
    nxu = nx; nxu2 = nx;
end
S.name = 'Symmetric LCD';
S.sampling = GaussianSamplingLCD();
S.sampling.setNumSamplesByFactor(10);
[samples,~,numSamples] = S.sampling.getStdNormalSamples(nxu);
[samples2,~,numSamples2] = S.sampling.getStdNormalSamples(nxu2);

% start plot
figure(1)
ax = gca;

% start filtering
for j = 1:nMC
    j
    % artificial measurement noise 
    meanInX = 0; meanInY = 0; varInX = 0; varInY = 0;
    tao = 200;
    for k = 1:nSamples
        % generate measurement
        thetas = 2*pi*rand(1,nMeas);
        us = random(pd1,1,nMeas);
        nu = rand(1,nMeas); 
        ss = ones(1,nMeas); 
        ss(nu > 0.7) = random(pd2,1,numel(ss(nu > 0.7)));
        meas = [ss.*(1 - us)*a_ref.*cos(thetas); ss.*(1 - us)*b_ref.*sin(thetas); us*h_ref] + sig_r*randn(3,nMeas);
        % translate and rotate
        meas = pos_ref + R_ref*meas;
    
        % start filter
        if k == 1
            % do position initialization
            pos = zeros(3,1);
            pos(1:2) = mean(meas(1:2,:),2);
            pos(3) = min(meas(3,:));
            % init orientation
            c = cov(meas(1,:),meas(2,:));
            [v,d] = eig(c); [~,ii] = max(diag(d));
            or = atan2(v(2,ii),v(1,ii));
            % init extent
            R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            mm = R_rot'*(meas - pos);
            a = c1(max(abs(mm(1,:))),0,'lower');
            b = c1(abs(max(mm(2,:))),0,'lower');
            h = c1(max(mm(3,:)) - min(mm(3,:)),0,'lower');
            % initialization
            X(:,k,j) = [pos; or; a; b; h];
            P(:,:,k,j) = blkdiag(c, 2*eye(5));   
        else
            % prediction
            X(:,k,j) = X(:,k-1);
            P(:,:,k,j) = P(:,:,k-1) + Q;
        end
    
        % update elliptic cone
        tic
        [X(:,k,j),P(:,:,k,j),meanInX,meanInY,varInX,varInY] = elliptic_cone_update(X(:,k,j),P(:,:,k,j),meas,sig_r,samples,numSamples,samples2,numSamples2,n_upd,...
                                                              meanInX,meanInY,varInX,varInY,tao,artificial_noise,filter,source);
        times(k,j) = toc;

        if k == 50
            aaa=1;
        end
    
        %% plot reference, measurement, and estimate
        if do_plot
            plot_results
        end
    end
end

% plot results
figure(2)
tiledlayout(5,1)
set(gcf,'Color','w')

% constrained values
X(5:7,:,:) = c1(X(5:7,:,:),0,'lower');

% mean values
X_mean = mean(X,3);

% calculate rmses
rmse_pos = sqrt((X_mean(1,:) - pos_ref(1)).^2 + (X_mean(2,:) - pos_ref(2)).^2 + (X_mean(3,:) - pos_ref(3)).^2);
rmse_or = sqrt((X_mean(4,:) - or_ref).^2);
rmse_a = sqrt((X_mean(5,:) - a_ref).^2);
rmse_b = sqrt((X_mean(6,:) - b_ref).^2);
rmse_h = sqrt((X_mean(7,:) - h_ref).^2);

%% plot position rmse
nexttile
plot(1:nSamples,rmse_pos,'LineWidth',2)
ylabel('RMSE','Interpreter','latex')
title('Position RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

%% plot orientation rmse
nexttile
plot(1:nSamples,rmse_or,'LineWidth',2)
ylabel('RMSE','Interpreter','latex')
title('Orientation RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

%% plot major half axes rmse
nexttile
plot(1:nSamples,rmse_a,'LineWidth',2)
ylabel('RMSE','Interpreter','latex')
title('Major half axes a RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

%% plot minor half axes rmse
nexttile
plot(1:nSamples,rmse_b,'LineWidth',2)
ylabel('RMSE','Interpreter','latex')
title('Minor half axes b RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

%% plot height rmse
nexttile
plot(1:nSamples,rmse_h,'LineWidth',2)
xlabel('sample','Interpreter','latex')
ylabel('RMSE','Interpreter','latex')
title('Height RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

% display mean calculation time
disp(['Mean calculation time: ',num2str(mean(times(:)))])

% end of function code
% *************************************************************************
%
%
% *************************************************************************