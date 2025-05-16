function [imNew,imT,imP,u] = rk_seg_convergencefinal_slower(im,mask,mu,lambda,gamma,rho2,dProxConst,betaProxConst,uProxConst,c1ProxConst,c2ProxConst)
addpath('..');
%%% this function doesn't have the problem of weird lines coresponding to
%%% patches in the results.
[n,m] = size(im);

if nargin < 2 %%%% if mask is not predefined, initialise as zeros
    u = zeros(h,w); wx = zeros(h,w); wy = zeros(h,w);
else
    mask = double(mask);
    u = mask;
    [wx,wy] = gradient(u);    
end

patchsize = 4; overlap = 2;

%%%
sig=4;
K = construct_Gaussian_kernel(patchsize,patchsize,1,sig,0);
[P,~] = consP_2d_arbDim(patchsize,patchsize,patchsize,patchsize,1e-4,12);

 %%% Parameters
 rho1=2; 
 alpha=1; %regularisor on ||beta||_1
 %gamma = 1e-5; %regulariser on c'Kc
 eta = 1000; %fitting parameter
 iota = 1000;
 
param.r1 = rho1; %sparsity beta admm
param.alpha = alpha;
param.gamma = gamma;
param.eta = eta;

[mm1,nn1] = size(K);
 [mm2, nn2] = size(P);

%%%%
Matrix.C1 = (eta)*(K'*K) + 2*gamma*K;
Matrix.C2 = (eta)*K'*P;
Matrix.B1 = (eta)*P'*K;
B1 = Matrix.B1;
Matrix.B2 = (eta)*(P'*P) + rho1*eye(nn2);
Matrix.B2I = inv(Matrix.B2);
B2I = Matrix.B2I;
Matrix.A1 = Matrix.C1 - Matrix.C2*Matrix.B2I*Matrix.B1;
Matrix.A1I = inv(Matrix.A1);
A1I = Matrix.A1I;
Matrix.A2 = Matrix.C2*Matrix.B2I;
A2 = Matrix.A2;
matrixA = eta*(K'*K) + 2*gamma*K + dProxConst*eye(nn1);
invMatrixA = inv(matrixA);
%Matrix.A3 = inv(Matrix.C1 - Matrix.A2);


h=n; w=m;
%%% grid for patch 
gridy = 1:patchsize - overlap: w-(mod(w,patchsize-overlap)+1+patchsize-overlap);
gridy = setdiff(gridy, [(w-patchsize+1):w]);  % delete rest elements
gridx = 1:patchsize - overlap: h-(mod(h,patchsize-overlap)+1+patchsize-overlap);
gridx = setdiff(gridx, [(h-patchsize+1):h]); % delete rest elements
gridx = [gridx h-patchsize+1]; gridy = [gridy w-patchsize+1];

A = zeros(h,w);
B = zeros(h,w);
Pb = zeros(h,w);
Kd = zeros(h,w);
b1Patches = zeros(nn2,size(gridy,2)*size(gridx,2));
thetaPatches = zeros(nn2,size(gridy,2)*size(gridx,2));
cPatches = zeros(nn1,size(gridy,2)*size(gridx,2));
betaPatches = zeros(nn2,size(gridy,2)*size(gridx,2));
oldbetaPatches = zeros(nn2,size(gridy,2)*size(gridx,2));

lambda1 = psf2otf([1,-1],[h,w]); lambda2 = psf2otf([1;-1],[h,w]);
lambda1Conj = conj(lambda1); lambda2Conj = conj(lambda2);
lambdaSum = lambda1.*lambda1Conj + lambda2.*lambda2Conj;

c1 = sum(sum(im.*u))/sum(u(:));
c2 = sum(sum(im.*(1-u)))/(sum(sum(1-u))+eps);
param.c1 = c1;
param.c2 = c2;
b2x = zeros(h,w); b2y = zeros(h,w);


%%
count = 0;
E=[];
tic
for iter = 1:20
    fprintf('ADMM steps: %d\n', iter);
    count = count+1;
    %compute c,beta for each patch.
    patchCounter = 0;
    
    A = zeros(h,w);
    B = zeros(h,w);
    Pb = zeros(h,w);
    Kd = zeros(h,w);
    
    
    %if count < 20 || mod(count,10)==0
for i=1:length(gridx)
for j=1:length(gridy)
    patchCounter = patchCounter+1;
    xx = gridx(i); yy = gridy(j);
    imPatch = im(xx:xx+patchsize-1,yy:yy+patchsize-1);
    uPatch = u(xx:xx+patchsize-1,yy:yy+patchsize-1);
    wxPatch = wx(xx:xx+patchsize-1,yy:yy+patchsize-1);
    wyPatch = wy(xx:xx+patchsize-1,yy:yy+patchsize-1);
    
    imvec = imPatch(:);
    uVec = uPatch(:);
    
    b1 = b1Patches(:,patchCounter);
    theta = thetaPatches(:,patchCounter);
    beta = betaPatches(:,patchCounter);
    c = cPatches(:,patchCounter);
    
    %%% c beta problem
%     CR = eta*K'*imvec;
%     BR = eta*P'*imvec + rho1*(theta+b1);
%    
%     c = A1I*(CR-A2*BR);
%     beta = B2I*(BR-B1*c);
%    

    %c problem
    CR = eta*(K'*imvec) - (eta+2*lambda)*(K'*P)*beta + dProxConst*c + 2*lambda*K'*(c1*uVec + c2*(1-uVec));
    c = invMatrixA*CR;
    
    a1 = K*c;
    
    %%% beta problem
    olderBeta = oldbetaPatches(:,patchCounter);
    betaHat = beta + (beta - olderBeta);
    PBetaHat = P*betaHat;
    gBetaHat = 1./(1+iota.*(PBetaHat).^2); %this is vector form
    
    pHat1 = -P'*(imvec - (K*c + P*betaHat));
    pHat2 = 0; %pHat2 = rho1*(theta - betaHat + b1);
    
    wForP = abs(wxPatch) + abs(wyPatch); wForP = wForP(:);
    pHat3 = -2*mu*iota*P'*((gBetaHat.^2).*(PBetaHat).*wForP);
    pHat4 = 2*lambda*P'*(a1 + PBetaHat - c1*uVec - c2*(1-uVec));
    pHat = pHat1 + pHat2 + pHat3 + pHat4;
    
    beta = (1/(rho1+betaProxConst))*(betaProxConst*betaHat - pHat + rho1*(theta+b1));
    
    
    
    %%% theta problem
   p1 = sign(beta - b1);
   p2 = abs(beta - b1);
   theta = p1.*max(p2 - alpha/rho1, 0);
   b1 = b1 + (theta - beta);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    
    b1Patches(:,patchCounter) = b1;
    thetaPatches(:,patchCounter) = theta;
    betaPatches(:,patchCounter) = beta;
    
    a1 = K*c;
    a2 = P*beta;
    
    a1 = reshape(a1,[patchsize,patchsize]);
    a2 = reshape(a2,[patchsize,patchsize]);
    newPatch = a1 + a2;
    
    A(xx:xx+patchsize-1, yy:yy+patchsize-1) = A(xx:xx+patchsize-1, yy:yy+patchsize-1) + newPatch;
    B(xx:xx+patchsize-1, yy:yy+patchsize-1) = B(xx:xx+patchsize-1, yy:yy+patchsize-1) + 1; %OVERLAP

    Kd(xx:xx+patchsize-1, yy:yy+patchsize-1) = Kd(xx:xx+patchsize-1, yy:yy+patchsize-1) + a1;
    Pb(xx:xx+patchsize-1, yy:yy+patchsize-1) = Pb(xx:xx+patchsize-1, yy:yy+patchsize-1) + a2;

end
end
oldbetaPatches = betaPatches;

B(B==0) = 1;
J = A./B;
PsiBeta = Pb./B;


g = 1./(1+iota.*(PsiBeta).^2);
    
    
JorIm = J;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% c1 c2 problem - done
c1N = sum(sum(JorIm.*u)); c1N = c1N + c1ProxConst*c1;
c1D = sum(u(:)); c1D = c1D + c1D + c1ProxConst*c1;
c1 = c1N/c1D;

c2N = sum(sum(JorIm.*(1-u))); c2N = c2N + c2ProxConst*c2;
c2D = sum(sum(1-u)); c2D = c2D + c2D + c2ProxConst*c2;
c2 = c2N/c2D;

param.c1 = c1;
param.c2 = c2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% u problem - done
   uRa = rho2.*(lambda1Conj.*fft2(wx+b2x) + lambda2Conj.*fft2(wy+b2y));
   uRb = lambda.*fft2((JorIm-c1).^2 - (JorIm-c2).^2);
   uRc = uProxConst.*fft2(u);
   uR = uRa - uRb + uRc;
   uL = rho2.*lambdaSum + uProxConst;
   u = real(ifft2(uR./uL));
   u = min(max(u,0),1);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %%%% w problem - remains the same.
   [ux,uy] = gradient(u);
   w1a = sign(ux-b2x);
   w2a = sign(uy-b2y);
   w1b = abs(ux-b2x);
   w2b = abs(uy-b2y);
   wx = w1a.*max(w1b - g.*(mu/rho2),0);
   wy = w2a.*max(w2b - g.*(mu/rho2),0);
   %wx = min(max(wx,0),1); wy = min(max(wy,0),1);

   b2x = b2x + wx - ux;
   b2y = b2y + wy - uy;
   
    


   
   e1 = norm(im-J,2); %fitting term
   e20 = abs(wx)+abs(wy);
   e2 = mu*(g(:))'*e20(:); %g'*w
   e30 = (J-c1).^2;
   e3 = lambda*u(:)'*e30(:); %
   e40 = (J-c2).^2;
   e4 = lambda*(1-u(:))'*e40(:);
   E0 = e1+e2+e3+e4;
   E=[E;E0];

end

toc
figure; plot(E);
B(B==0) = 1;
imNew = A./B;

imT = Kd./B;
imP = Pb./B;
end
%% 
function [H, HH] = consP_2d_arbDim(n_gridx, n_gridy, N_gridx, N_gridy, ksi, m)

   % ksi: control the amount of smooth
   % larger ksi, more smoothing
   

 
   %----------coarse grid--------%
    tx = [0:(n_gridx-1)]'/(n_gridx-1);
    ty = [0:(n_gridy-1)]'/(n_gridy-1);
%     tx = [1:n_gridx]'/n_gridx;
%     ty = [1:n_gridy]'/n_gridy;
    
    X = repmat(tx,n_gridy,1);
    Y = [];

    for i = 1:n_gridy
        e = repmat(ty(i), n_gridx,1);
        Y = [Y; e];
    end

    %------fine grid--------------
    ttx = [0:(N_gridx-1)]'/(N_gridx-1);
    tty = [0:(N_gridy-1)]'/(N_gridy-1);
%     ttx = [1:N_gridx]'/N_gridx;
%     tty = [1:N_gridy]'/N_gridy;
    
    XX = repmat(ttx,N_gridy,1);
    YY = [];

    for i = 1:N_gridy
        ee = repmat(tty(i), N_gridx,1);
        YY = [YY; ee];
    end
    
    %====== function: H=1/2+1/pi*arctan(z/ksi) ===============%
    %====== z = cos(theta)*x + sin(theta)*y + c==============%
    %====== c = tx or ttx  =====================%
     n = n_gridx*n_gridy;
     N = N_gridx*N_gridy;
     %ksi = 1e-4;
    %m =4;
    theta = [0: 2*pi/m: 2*pi];
    %theta = [0: pi/m: pi];
    %====== construct H ===========%
%     grid = 19;
%     c = [0:grid]/grid; 
%     c = X;
%     grid = size(c,1);
    %grid = 25;
    grid = n;
    c = [0:(grid-1)]'/(grid-1);
    iter = 1;
    for j = 1:m
       for i = 1:grid
         z = cos(theta(j))*X + sin(theta(j))*Y + repmat(c(i),n,1);
         z = z/ksi;
         r = atan(z);
         H(:,iter) = repmat(1/2, n, 1) + 1/pi*r;
         iter = iter + 1;
       end
    end
    
   %====== construct HH =========%
    iter = 1;
    d = ttx;
      for j = 1:m
       for i = 1:grid
         z = cos(theta(j))*XX + sin(theta(j))*YY + repmat(c(i),N,1);
         z = z/ksi;
         r = atan(z);
         HH(:,iter) = repmat(1/2, N, 1) + 1/pi*r;
         iter = iter + 1;
       end
      end 
    
end
%% 
function K = construct_Gaussian_kernel(n_gridx, n_gridy, cls, sig, coeff1)

   %----------coarse grid--------%
    tx = [0:(n_gridx-1)]'/(n_gridx-1);
    ty = [0:(n_gridy-1)]'/(n_gridy-1);
%     tx = [1:n_gridx]'/n_gridx;
%     ty = [1:n_gridy]'/n_gridy;
    
    X = repmat(tx,n_gridy,1);
    Y = [];

    for i = 1:n_gridy
        e = repmat(ty(i), n_gridx,1);
        Y = [Y; e];
    end
    
    
    type=9; %3,8,9 all right
switch type
    case 1
        %%%%%% This case takes \xi in [0,1], partioned depending on cls,
        %%%%%% and exvaluates gaussian kernel exp( - | x - \xi_i )
     %====== function: K= a * exp{ - | x - \xi_i|^2  / 2sig^2} ===============%
    % ==== a = 1/(sqrt(2*pi)*sig)
    K = zeros(n_gridx*n_gridy,cls);
a = 1/(sqrt(2*pi)*sig);
xi = 0:1/(cls-1):1; % clsX1 vec between 0 and 1
for L = 1:cls
    hhh = exp( -( X - xi(L) ).^2/(2*sig^2) ).*exp( -( Y - xi(L) ).^2/(2*sig^2) ); %xi_i [0,1]
    
    K(:,L) = a*hhh;  
end
    case 2
        %%%% This evaluates thekernel exp{ -(X.^2 +Y.^2)/2sig^2 } on the
        %%%% discretized [0,1]x[0,1] grid
        n = n_gridx*n_gridy;
        
        
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        [repX,repY] = meshgrid(tx,ty);

        a = exp(-(repX.^2)/(2*sig^2));
        b = exp(-(repY.^2)/(2*sig^2));
        mEXP = a.*b;    
        
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff^2*mEXP;        
        
    case 3
        %%%% This evaluates thekernel exp{ -((X-X').^2 + (Y-Y').^2)/2sig^2 } on the
        %%%% discretized [0,1]x[0,1] grid
        n = n_gridx*n_gridy;
        
        [repX,repY] = meshgrid(tx,ty);
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        myX = repX - repX';
        myY = repY - repY';
        a = exp(-(myX.^2)/(2*sig^2));
        b = exp(-(myY.^2)/(2*sig^2));
        mEXP = a.*b;
        
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff^2*mEXP;  
        
    case 4
        %%%This evaluates exp( - (x^2 + y^2)) and repeats the 36x1 vec no
        %%%times
        no = 10;
        K = zeros(n_gridx*n_gridy,no);
        for L=1:no
            K(:,L) = exp(-(X.^2 + Y.^2)/(2*sig^2));
        end
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff*K;
        
    case 5
        %Consider only 6x6 grid
        X6 = reshape(X,[n_gridx,n_gridy]);
        Y6 = reshape(Y,[n_gridx,n_gridy]);
        
        a = exp(-(X6.^2)/(2*sig^2));
        b = exp(-(Y6.^2)/(2*sig^2));
        mEXP = a.*b;
        
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff*mEXP; 
        
    case 6
        %We look at T_ij = K(x_i,y_j) ~ exp( (x_i - y_j)^2 / 2sigma^2 );
        %Use our X and Y 36x1 vecs
        n= n_gridx*n_gridy;
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        
        mExp = exp( ( - (repX-repY).^2) / (2*sig^2) );
        coeff = 1/(sqrt(2*pi)*sig);
        
        K = coeff*mExp;
        
    case 7
        %We look at T_ij = K(x_i,y_j) ~ exp( (x_i - x_j)^2 / 2sigma^2 ) exp( (y_i - y_j)^2 / 2sigma^2 )
        n= n_gridx*n_gridy;
        %K will be 36x36
        for L=1:n
        Kx(L,:) = exp( - (( X(L) - X(:)).^2)/(2*sig^2));
        Ky(L,:) = exp( - (( Y(L) - Y(:)).^2)/(2*sig^2));
        end
        
        mExp = Kx.*Ky;
        coeff = 1/(sqrt(2*pi)*sig);
        coeff = coeff^2;
        K = coeff.*mExp;
        
    case 8
        %We look at now considering advice after 36x2
        XX = [X Y];
        YY = [X Y];
        
        n= n_gridx*n_gridy;
        %K will be 36x36
        %K(i,j) = exp( - ( (XX(i)-YY(i)).^2)/(2*sig^2)) * j part
        %T_12 = exp( - ( (XX(1,1)-YY(2,1
        
        for i=1:n
        for j=1:n
           a = XX(i,1) - YY(j,1);
           b = XX(i,2) - YY(j,2);
           tao = a.^2 + b.^2;
           K(i,j) = exp(-tao./(2*sig^2));
        end
        end
        
%         for i=1:n
%         for j=1:n
%         Kx(i,j) = exp( - (( XX(i,1) - YY(i,1)).^2)/(2*sig^2));
%         Ky(i,j) = exp( - (( XX(j,2) - YY(j,2)).^2)/(2*sig^2));
%         end
%         end
        
        %mExp = Kx.*Ky;
        coeff = 1/(sqrt(2*pi)*sig);
        coeff = coeff^2;
        K = coeff.*K;
        
    case 9
        n = n_gridx*n_gridy;
        
        [repX,repY] = meshgrid(tx,ty);
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        a = repX - repX';
        b = repY - repY';
        tao = sqrt(a.^2 + b.^2);
        tao(find(tao==0)) = eps;
        K = exp(-tao.^2./(2*sig^2));
        if coeff1 == 1
            coeff = 1;
        else
        coeff = 1/(sqrt(2*pi)*sig);
        end
        coeff = coeff^2;
        K = coeff.*K;
end
end