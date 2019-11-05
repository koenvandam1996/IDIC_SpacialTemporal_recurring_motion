%% help In this file (in order)
%   - All image from folder Images are read ine
%   - The Correlation options are set (e.g. coarse graining steps, margins)
%   - The basis functions and corresponding DOF are created
%   - The globalDIC2D.m function is called
%   - The results are gathered and plotted
close all; clear
%% m-file options
fontsize = 16;
printpdf = true;
load('plotop.colormap.mat');
%% Loading the images for both the front and the back part
Dfront='C:\Users\s149743\Pictures\12-9\backmagnfull\att1 - Copy (2)';%134
Sfront = dir(fullfile(Dfront,'Basler_scA1400-30gm__21502907__20190912_090545803_0*.tiff')); % pattern to match filenames.
Nfront = natsortfiles({Sfront.name});
imstfr = 134;
NoI = 3;
for k = imstfr:imstfr+NoI%length(Nfront)
    File = fullfile(Dfront,char(Nfront(k)));
    I = imread(File);
    if size(I,3)==3
        I = rgb2gray(imread(File));
    end
    Manifoldfront(1:392,:,k-133) = I(9:400,:);
end

Dback ='C:\Users\s149743\Pictures\12-9\front\att1 - Copy (3)';%104
Sback = dir(fullfile(Dback,'Basler_scA1400-30gm__21502907__20190912_085153200_0*.tiff')); % pattern to match filenames.
Nback = natsortfiles({Sback.name});
imstback = 104;
for k = imstback:imstback+NoI%length(Nback)
    File = fullfile(Dback,char(Nback(k)));
    I = imread(File);
    if size(I,3)==4
        I = I(:,:,1:3);
    end
    if size(I,3)==3
        I = rgb2gray(imread(I));
    end
    Manifoldback(1:392,:,k-103)= I(1:392,:);
end

% Create position columns (keep images the same size)
pixelsize = [1 1]; 
[nfront, mfront] = size(Manifoldfront(:,:,1));
x1front = linspace(1,mfront*pixelsize(1),mfront);
y1front = linspace(1,nfront*pixelsize(2),nfront);
[nback, mback] = size(Manifoldback(:,:,1));
x1back = linspace(1,mback*pixelsize(1),mback);
y1back = linspace(1,nback*pixelsize(2),nback);
% center the axes
x1front = x1front - mean(x1front);
y1front = y1front - mean(y1front);
x1back = x1back - mean(x1back);
y1back = y1back - mean(y1back);
x1tot = linspace(1,(mback+mfront)*pixelsize(1),(mback+mfront));

%% GDIC Options
% Region Of Interest and Masking
margin.xleft        = 0*pixelsize(1);
margin.xright       = 50*pixelsize(1);
margin.y            = 0*pixelsize(2);
maskfront=imread(fullfile(Dfront,'mask.png'));
if size(maskfront,3)==3
    maskfront            = rgb2gray(imread(fullfile(Dfront,'mask.png')));
end
maskback=imread(fullfile(Dback,'mask.png'));
if size(maskback,3)==3
    maskback            = rgb2gray(imread(fullfile(Dback,'mask.png')));
end
options.ROIfront(1)      = x1front(1)   + margin.xleft;
options.ROIfront(2)  	 = x1front(end) - margin.xright;
options.ROIfront(3)      = y1front(1)   + margin.y;
options.ROIfront(4)      = y1front(end) - margin.y;
options.ROIback(1)      = x1back(1)   + margin.xleft;
options.ROIback(2)  	= x1back(end) - margin.xright;
options.ROIback(3)      = y1back(1)   + margin.y;
options.ROIback(4)      = y1back(end) - margin.y;
options.maskfront        = find(maskfront(9:400,:)==255);
options.maskback       = find(maskback(1:392,:)==255);
options.coarsesteps = [3 2 1]; 
options.background  = 0;
% options.convcrit   = 8e-4;
options.factorback  = ones(1,NoI*2+1);
options.maxiter     = 400;
% options.list        = true;
options.mc          = 2; %specify at which coarsegrain step to use mc
options.interpmeth  = 'cubic';
%options.basis      = 'chebyshev';
options.gradg       = 'G6';
options.strain      = 'greenlagrange';%'small','membrane','none','logarithmic','greenlagrange','euler-almansi'
options.datatype    = 'single';
options.slider(1)       = 910-margin.xright;
options.slider(2)       = 910-margin.xright;
Manifoldtot =[Manifoldback(:,:,1),Manifoldfront(:,options.slider(1)+margin.xright*2+1:mback,1)];

%% Basis and DOF
% Generate a full basis for polynomial or chebyshev functions (polynomial degree);
phi = dofbuild_poly(4);

% Initial Guess
Ndof = size(phi,1);
u    = zeros(Ndof,1);
u(1) = 26*pixelsize(1);%
u=[26.165085;0.30546516;-0.017617902;-0.94123757;6.5929527;-0.38367406;-1.9603693;-0.36585912;1.6118947;-0.18152377;-5.3541126;-1.6627792;3.1715448;-1.0541872;-0.20105629;-0.093309276;-0.72670704;1.6936055;-4.5462489;0.77148622;10.819376;-0.30850366;-1.7653918;-0.81109941;1.3287251;0.45012948;-2.1772091;0.32886606;3.8423023;1.5706958];
%% Start the digital image correlation
% clear the old logfile
if isfield(options,'logfile')
    delete(options.logfile);
end 

gdic = globalDIC2Dcomb(x1back,y1back,u,phi,Manifoldback,x1front,y1front,Manifoldfront,x1tot,options);

%% Plot the results
% Extract some data from the gdic structure
Ux = gdic.Ux;
Uy = gdic.Uy;
x  = gdic.x;
y  = gdic.y;

% Displacement Fields
Ux(gdic.mask)=NaN;
Uy(gdic.mask)=NaN;
figure(2)
subplot(2,2,1);
imagesc(Manifoldtot);%x1back,y1back,crs.Manf(:,:,1));
colormap(gray(256))
title('Image 1','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
subplot(2,2,2);
imagesc(Manifoldtot);%x1back,y1back,crs.Manf(:,:,2));
colormap(gray(256))
title('Image 2','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax1=subplot(2,2,3);
imagesc(Manifoldtot);
colormap(ax1,gray(256));
title('Displacement Field Ux','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(Ux,'alphadata',~isnan(Ux));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Displacement[\mum]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(2,2,4);
imagesc(Manifoldtot);
colormap(ax1,gray(256));
title('Displacement Field Uy','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(Uy,'alphadata',~isnan(Uy));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Displacement[\mum]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
% globalDICplot(x,y,Ux,Uy,plotop);
if printpdf
    savepdf('dic_DisplacementFields')
end
