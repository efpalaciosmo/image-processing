% out_test.m for 4D Display ONLY
%                   save Out_data I0 phi phi0 slice_start slice_end
%================ Post view:  slice_start / slice_end known by now

if exist('phi') && exist('phi0') && exist('I')
    save Out_data I phi phi0 slice_start slice_end
    disp('Just Saved DATA = Out_data.mat')
elseif exist('Out_data.mat'); load  Out_data.mat;
    disp('USING the loaded DATA = Out_data.mat')
end
fprintf('[%d:%d] - see phi, I >> >> \n', slice_start, slice_end);
            close all; I0=double(I); figure; disp(' 3D data plotting ...')
[m n p0]=size(phi);  oc2=1; M=max(I0(:)); disp(' ') %%D=zeros(m,n,p);   
  D=zeros(m,n,slice_end+2-slice_start); scrsz = get(0,'ScreenSize'); 
    fprintf(' Check Original Slice Fig.  %%%');organ=0;
for s=(slice_start):(slice_end) %% 1:p0
    fprintf('\b\b\b%3d',s);  whitebg('black')
    c=contour(squeeze(phi(:,:,s)),[-eps 0 eps],'r','Linewidth',1);
    if size(c,2)>0
        organ=organ+1; %  
        posX = 100; posY = 100;  rows=m; columns=m; 
        set(gcf,'Position',[posX posY columns rows]); % To save figure as MATRIX
        set(gca,'units','pixels');
        set(gca,'units','normalized','position',[0 0 1 1]);
        f=getframe(gcf);  new=f.cdata;  
        new=imclearborder(new,8);  new=rgb2gray(new);
        D(:,:,organ)=new; clf 
    end
    cs=c;
end
   close
   
    D=D(:,:,1:organ); f = figure('visible','on'); close %close all, 
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b: %d Slices in range ', organ);
fprintf('[%d:%d] - see phi, I >> \n Done....\n', slice_start,slice_end);
%% ====================================================================
   ss=scrsz(3:4)/1.8; clear p1;
 figure('Position',[9 9 ss]);    %--------------------------------------
   p1(:,:,1,:)=I0;  montage1(p1);  title('Input 3D image I_0');
 figure('Position',[19 90 ss]);  %--------------------------------------
   p1(:,:,1,:)=-phi0; montage1(p1(:,:,1,32:67)); 
   title('\phi_0 on slices 32 - 67: Note Number 33/52/65 are set.')

%% montage Data(:,:,1,:) or Data(:,:,1:3,:)  %--------------------
figure('Position',[29 180 ss]), clear K; K(:,:,1,:)=D; montage1(K)
         title('Segmented boundaries')

%% ABOVE for my_test.m  BELOW for 3D viewer
 figure('Position',[440 9 ss]); Ip=double(I0);  fprintf('\b\b')
             for s=1:p0, 
                 fprintf('\b\b\b%3d',s);
                 Ip(:,:,s)=imclearborder(Ip(:,:,s),8);
              if s<91, Ip(1:75,1:75,s)=0; end
              if s>90, Ip(1:95,1:95,s)=0; end
             end; disp(' |'); clear p
p=patch(isosurface(Ip, 85), 'FaceColor', [1 0 0], 'EdgeColor', 'none');
  view([303 -18]);  daspect([1 1 1]);axis tight
    isonormals(Ip,p)
    camlight, camlight(-80,-10), lighting phong %gouraud %
title('Original Data Normals')
%% ----------
 figure('Position',[459 90 ss]); %--------------------------------------
p=patch(isosurface(phi0, .5), 'FaceColor', [1 0 1], 'EdgeColor', 'none');
 view(3); daspect([1 1 1]); a=[1 150 1 150 1 110]; axis(a)
    isonormals(phi0,p)
    camlight, camlight(-40,-10), lighting gouraud %phong %
title('Constructed \phi_0 level based on markers of 3 slices')
%% 
figure('Position',[479 180 ss]); %--------------------------------------
D1 = subvolume(-5*phi,[nan,nan,nan,nan,nan,nan]);
Ds = smooth3(D1);
hiso = patch(isosurface(Ds,2),...
    'FaceColor',  [1,.75,.65],...
    'EdgeColor','none');
hcap = patch(isocaps(D1,2),...
    'FaceColor','interp',...
    'EdgeColor','none');
view(37, -22); lightangle(45,30); 
    set(gcf,'Renderer','zbuffer'); lighting phong
    isonormals(Ds,hiso);set(hcap,'AmbientStrength',.6)
    set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)
title('Segmented 3D Object')
%======================================================================%%