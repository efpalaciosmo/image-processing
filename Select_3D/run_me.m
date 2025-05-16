%  Usage: [phi,slice_start,slice_end,phi0]=...
%              Selective_3D_seg(I, x1,y1,z1,x2,y2,z2,x3,y3,z3, alpha,tol,maxit);
%  Notation: Input 3D image I sized m X n X p [ie. p vertical slices of mxn]
%                 (xj,yj,zj) with z_num(j) markers (points) on Slice j=1,2,3
%                 maxit ------ max numer of outer iterations to achieve "tol".
%   Ouput:              phi -- 3D level function of the same size m X n X p
%     slice_start,slice_end -- Deteced object sitting between these 2 slices
%                 phi0 ------- Initial level set function phi0(m,n,p)
%  % International Journal of Computer Mathematics    [3D Segmentation]
%     "A Fast Algorithm for Automatic Segmentation ... by Active Surfaces" 
%           by J P Zhang, K Chen and D Gould, Vol 92 (6), pp.1251-1274, 2015
%  Model:  min_phi J(phi,c1,c2) = Fitting(phi,c1,c2) + alpha*Regularizer(phi)
%(c) http://www.liv.ac.uk/cmit (2015)
clc, clear; close all %%%%%%%%%%%%%%%%%
 load I_CTData150_30(110).mat; % load the initial 3D image (by dicomread etc)
 
 [nrow,ncol,nz]=size(I);    alpha=0.03*nrow*ncol;  alpha=0.01*255^2;  
 
%% Use "DatCursor" on Image Display  or  roiploy
    x1 =[  27.6313   36.9631   41.8018   33.5069];% Slice 1 Markers (x1,y1)
    y1 =[  83.6140   72.6491   76.5965   92.8246];
    z1 =    65; % --------------------Slice 1 (top)
    x2 =[  20.7189 30.0507 44.5668 36.6175 40.4194   28.6682]; % Slice 2
    y2 =[  88.8772 73.9649 75.2807 83.6140 94.1404  101.1579];
    z2 =    49; % --------------------Slice 2 (middle)
    x3 =[  26.5945   39.7281   49.0599   36.2719];% Slice 3 Markers (x3,y3)
    y3 =[  95.4561   83.6140   91.5088  102.9123];
    z3 =    33; % --------------------Slice 3 (bottom)
%% Main Program
    tol = 3.0e-8; maxit= 200;
 [phi,slice_start,slice_end,phi0]=...
     Selective_3D_seg(I, x1,y1,z1,x2,y2,z2,x3,y3,z3,alpha,tol,maxit);
  
%% 4D Image Display
disp('          START of 4D Display ...')
 out_test