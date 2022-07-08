% An TM-based alternate projection algorithm for arbitrary glare suppression 
% through scattering media, under phase-only modulation. 

% Author: Shengfu Cheng
% Date: July 08 2022

close all

%% Parameters
N = 64; %size of phase mask 
M = 160; %size of speckle field
L =  30; %side length of target square region
iters = 30; %iteration number
constrNum = 2; %two speckle-domain constraints (ER and HIO)
gamma = 0.8; %feedback coefficient of HIO
gpuFlag = 1; %use GPU or not

if ~exist('TM', 'var'); TM = generate_tm(M^2, N^2, 0); end
if ~gpuFlag && ~exist('TM_inv', 'var'); TM_inv = Tikinv(TM); end
if gpuFlag && ~exist('TM_inv_gpu', 'var')
    TM_gpu = gpuArray(TM); 
    TM_inv_gpu = Tikinv(TM_gpu); 
end

%% Target square region
target = ones(M, M, 'single'); target(M/2-L/2:M/2+L/2-1, M/2-L/2:M/2+L/2-1) = 0; 

Ein_zeroPha = exp(1i*zeros(N^2, 1, 'single')); 
Iout_zeroPha = reshape(abs(TM * Ein_zeroPha).^2, M, M);
I_target = target.*Iout_zeroPha; I_target_gpu = gpuArray(sqrt(I_target(:)));
figure('color', 'w', 'position', [150 250 400 350]), imagesc(I_target), colormap('hot'), title('Target pattern'), pause(1)

idxT = find(target==0); 
idxB = find(target~=0);
TR = ones(M, M, 'single'); TR(idxT) = 0; %target region
BR = ones(M, M, 'single'); BR(idxB) = 0; %complementary to target region
if gpuFlag; TR_gpu = gpuArray(TR(:)); end

%% Alternate projection algorithm ER and HIO constraints
sigmaCurve_ER = zeros(iters, 1); sigmaCurve_HIO = zeros(iters, 1); %integrated intensty in target region
etaCurve_ER = zeros(iters, 1); etaCurve_HIO = zeros(iters, 1); %suppression factor

tic

I_cst = repelem(abs(TM * Ein_zeroPha).*exp(1i*2*pi*rand(M^2, 1)), 1, constrNum); %Initial speckle field
I_previous_HIO = I_cst(:, 2); 

Ih1 = abs(I_cst(:, 1)).^2; sigmaCurve_ER(1) = sum(gather(Ih1(idxT))); etaCurve_ER(1) = mean(Ih1(idxT))/mean(Ih1(idxB)); 
Ih2 = abs(I_cst(:, 2)).^2; sigmaCurve_HIO(1) = sum(gather(Ih2(idxT))); etaCurve_HIO(1) = mean(Ih2(idxT))/mean(Ih2(idxB));  
        
for i=1:iters
    if gpuFlag
        A = TM_inv_gpu * I_cst; 
    else
        A = TM_inv * I_cst; 
    end
    
    A_pha = exp(1i*angle(A)); %phase-only constraint
    
    if gpuFlag
        I = TM_gpu * A_pha; 
    else
        I = TM * A_pha; 
    end
    
    %record
    if i> 1
        Ih1 = abs(I(:, 1)).^2; sigmaCurve_ER(i) = gather(sum(gather(Ih1(idxT)))); etaCurve_ER(i) = gather(mean(Ih1(idxT))/mean(Ih1(idxB))); 
        Ih2 = abs(I(:, 2)).^2; sigmaCurve_HIO(i) = gather(sum(gather(Ih2(idxT)))); etaCurve_HIO(i) = gather(mean(Ih2(idxT))/mean(Ih2(idxB))); 
    end
    
    %speckle-domain constraint
    I_cst(:, 1) = I(:, 1) .* TR_gpu;  % ER constraint 
    I_cst(:, 2) = I(:, 2); I_cst(idxT, 2) = I_previous_HIO(idxT) - gamma * I(idxT, 2); I_previous_HIO = I_cst(:, 2);  % HIO constraint
      
end

toc

%% glare suppression results
sigma_zeroPha = sum(Iout_zeroPha.*BR, 'all'); eta_zeroPha = mean(Iout_zeroPha(idxT))/mean(Iout_zeroPha(idxB));

Ein_ER = exp(1i*angle(A(:, 1))); Ein_HIO = exp(1i*angle(A(:, 2))); 
Iout_ER = abs(reshape(gather(TM_gpu * Ein_ER), M, M)).^2; Iout_HIO = abs(reshape(gather(TM_gpu * Ein_HIO), M, M)).^2;
eta_ER = mean(Iout_ER(idxT))/mean(Iout_ER(idxB)); eta_HIO = mean(Iout_HIO(idxT))/mean(Iout_HIO(idxB));
sigma_ER = sum(Iout_ER(idxT),'all'); sigma_HIO = sum(Iout_HIO(idxT),'all'); 

    
figure('color', 'w', 'position', [150 200 1200 450]),
clim1 = min([min(Iout_zeroPha(:)),min(Iout_ER(:)),min(Iout_HIO(:))]); %minimum image value
clim2 = max([max(Iout_zeroPha(:)),max(Iout_ER(:)),max(Iout_HIO(:))]); %maximum image value
Iout_zeroPha_norm = (Iout_zeroPha-clim1) ./ (clim2-clim1); Iout_ER_norm = (Iout_ER-clim1) ./ (clim2-clim1); Iout_HIO_norm = (Iout_HIO-clim1) ./ (clim2-clim1); 
clim1 = 0; clim2 = 1;  

subplot(131), imshow(Iout_zeroPha_norm, []); colormap(hot); caxis([clim1, clim2]); title(sprintf('All on, eta: %.3f', eta_zeroPha), 'fontsize', 14)
hold on, plot([M/2-L/2, M/2+L/2-1, M/2+L/2-1, M/2-L/2, M/2-L/2], [M/2+L/2-1, M/2+L/2-1, M/2-L/2, M/2-L/2, M/2+L/2-1], 'w--', 'linewidth', 1.8); hold off

subplot(132), imshow(Iout_ER_norm, []), colormap(hot); caxis([clim1, clim2]); title(['ER, ', '$$\eta=$$', sprintf('%.3f', eta_ER)], 'Interpreter','latex', 'fontsize', 18)
hold on, plot([M/2-L/2, M/2+L/2-1, M/2+L/2-1, M/2-L/2, M/2-L/2], [M/2+L/2-1, M/2+L/2-1, M/2-L/2, M/2-L/2, M/2+L/2-1], 'w--', 'linewidth', 1.8); hold off

subplot(133), imshow(Iout_HIO_norm, []), colormap(hot); caxis([clim1, clim2]); title(['HIO, ', '$$\eta=$$', sprintf('%.3f', eta_HIO)], 'Interpreter','latex', 'fontsize', 18)
hold on, plot([M/2-L/2, M/2+L/2-1, M/2+L/2-1, M/2-L/2, M/2-L/2], [M/2+L/2-1, M/2+L/2-1, M/2-L/2, M/2-L/2, M/2+L/2-1], 'w--', 'linewidth', 1.8); hold off

h=colorbar('eastoutside','fontsize',15);
set(h,'Position', [0.93 0.20 0.015 0.65]);  

%% evoluation curves
figure('color', [1 1 1], 'position', [200 200 1000 450]), 

subplot(121), plot(1:iters, sigma_zeroPha*ones(1, iters),  'k--', 1:iters, sigmaCurve_ER,  'b-', 1:iters, sigmaCurve_HIO, 'r-', 'LineWidth', 3); legend('Initial', 'ER', 'HIO', 'Fontname', 'Times New Roman'); 
set(gca,'FontSize',18, 'LineWidth', 2), xlabel('Iterations', 'fontsize', 24, 'Fontname', 'Times New Roman'), ylabel('Integrated intensity', 'fontsize', 24, 'Fontname', 'Times New Roman'); ylim([0 1.05*sigma_zeroPha]); 

subplot(122), semilogy(1:iters, eta_zeroPha*ones(1, iters), 'k--', 1:iters, etaCurve_ER,  'b-', 1:iters, etaCurve_HIO, 'r-', 'LineWidth', 3); legend('Initial', 'ER', 'HIO', 'Fontname', 'Times New Roman');
set(gca,'FontSize',18, 'LineWidth', 2), xlabel('Iterations', 'fontsize', 24, 'Fontname', 'Times New Roman'), ylabel('Suppression factor', 'fontsize', 24, 'Fontname', 'Times New Roman'); ylim([0 1.2*eta_zeroPha]);   
 
