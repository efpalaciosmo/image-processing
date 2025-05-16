clear; close all; clc
disp('Reproducible Kernel Hilbert Space Based Global and Local Image Segmentation')
disp('by Liam Burrows, Weihong Guo, Ke Chen and Francesco Torella')
disp('Inverse Problems and Imaging, 2021, 15(1): 1-25 doi: 10.3934/ipi.2020048')
disp('    (Global version by a single functional)')
%% Image
%%
load im2_40
 
fprintf('\n INPUT image im of size %d x %d\n', size(im))
%%
gamma = 1e-6; %reg of
iota = 1000;
lambda = 1e-6;  lambda = 2e-6;  % alpha in model
mu = 1e-3;
rho2 = 1e-9;

dProxConst = 1e-9;
betaProxConst = 10;
uProxConst = 4.0e-6;
c1ProxConst = 1e-9;
c2ProxConst = 1e-9;

[imNew,imT,imP,u] = rk_seg_global(im,mask,mu,lambda,gamma,...
         rho2,dProxConst,betaProxConst,uProxConst,c1ProxConst,c2ProxConst);

%% figure and save
threshold = 0.5;
FigH = figure('Position', get(0, 'Screensize'));
imagesc(imNew); colormap gray; axis off; axis image; %title("uPc = " + uProxConst + ", lambda = " + lambda);
hold on; contour(u,[threshold,threshold],'r','LineWidth',2);
saveas(gcf,'output.png');
 
fprintf(' Output image u of size %d x %d:   paras alpha=%.2e gamma=%.2e\n',...
         size(u),lambda,gamma)
