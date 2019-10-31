%% help In this file (in order)
%   - All image from folder Images are read ine
%   - The Correlation options are set (e.g. coarse graining steps, margins)
%   - The basis functions and corresponding DOF are created
%   - The IDIC_Calculation.m function is called
%   - The results are gathered and plotted
close all; clear
%% m-file options
fontsize = 16;
printpdf = true;
load('plotop.colormap.mat');
%% Loading the images (from folder)
R='C:\Users\s149743\Documents\globalDIC - v1.0.rc - Release Candidate\Images\new';
imst=8;
NoI=2;
S = dir(fullfile(R,'roll_*.png'));
N = natsortfiles({S.name});
% read images in black and white
for k = imst:imst+NoI
    File = fullfile(R,char(N(k)));
    I = imread(File);
    if size(I,3)==3
        I = rgb2gray(imread(File));
    end
    Manifold(:,:,k-imst+1) = I;
end

% size of images and create pixel position columns
pixelsize = [1 1]; % Pixel dimensions (in um)
[n, m] = size(Manifold(:,:,1));
x1 = linspace(1,m*pixelsize(1),m);
y1 = linspace(1,n*pixelsize(2),n);

% center the axes
x1 = x1 - mean(x1);
y1 = y1 - mean(y1);

%% DIC Options, Working of each function is explained in help.mat
% Region Of Interest and Masking
margin.xleft        = 0*pixelsize(1);
margin.xright       = 150*pixelsize(1);
margin.y            = 0*pixelsize(2);
mask=imread(fullfile(R,'mask.png'));
if size(mask,3)==3
    mask            = rgb2gray(imread(fullfile(R,'mask.png')));
end
options.ROI(1)      = x1(1)   + margin.xleft;
options.ROI(2)  	= x1(end) - margin.xright;
options.ROI(3)      = y1(1)   + margin.y;
options.ROI(4)      = y1(end) - margin.y;
options.mask        = find(mask==255);
options.coarsesteps = [3 2 1]; 
options.background  = 1;
% options.convcrit   = 2e-6;
options.factor      = ones(1,length(Manifold(1,1,:))-1);
options.maxiter     = 400;
% options.list       = true;
options.interpmeth  = 'cubic';
% options.normalized = false;
% options.basis      = 'chebyshev';
options.gradg       = 'G6';
% options.verbose    = 0;
% options.logfile    = 'globalDIC2D.log' ;
% options.comment    = false
%options.debug       = true ;
% options.fulloutput = true;
options.strain      = 'greenlagrange';%'small','membrane','none','logarithmic','greenlagrange','euler-almansi'
options.mc          = 2; %specify at which coarsegrain step to use mc


%% Basis and DOF
% Generate a full basis for polynomial or chebyshev functions (polynomial degree);
phi = dofbuild_poly(4);

% Initial Guess
Ndof = size(phi,1);
u    = zeros(Ndof,1);
u(1) = 81;%

%% Start the digital image correlation
% clear the old logfile
if isfield(options,'logfile')
    delete(options.logfile);
end 
gdic = IDIC_Calculation(x1,y1,u,phi,Manifold,options);


%% Plot the results
% Extract some data from the gdic structure
% Displayed as average fields
Ux = gdic.Ux*mean(gdic.factors); 
Uy = gdic.Uy*mean(gdic.factors); 
x  = gdic.x;
y  = gdic.y;


%% Displacement Fields
Ux(gdic.mask)=NaN;
Uy(gdic.mask)=NaN;

figure
subplot(2,2,1);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(gray(256))
title('Image 1','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
subplot(2,2,2);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(gray(256))
title('Image 2','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax1=subplot(2,2,3);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Displacement Field Ux','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,Ux,'alphadata',~isnan(Ux));
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
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Displacement Field Uy','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,Uy,'alphadata',~isnan(Uy));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Displacement[\mum]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
if printpdf
    savepdf('dic_DisplacementFields')
end

%% plot reflection image
if options.background>0
    figure(5)
    ax1=gca;
    imagesc(x1,y1,crs.Manf(:,:,1));
    colormap(ax1,gray(256));
    title('background','Fontsize',fontsize);
    xlabel('x [\mum]');
    ylabel('y [\mum]');
    ax2 = axes;
    imagesc(x,y,gdic.background,'alphadata',~isnan(gdic.background));
    colormap(ax2,plotop.colormap);
    h=colorbar;
    h.Label.String='Grayscale[-]';
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    linkaxes([ax1,ax2],'xy');
    hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));
end
%% Residual image
figure(6)
ax1=subplot(3,2,1);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Av residual Image: cg 3, it 1','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,residual01(:,:,1),'alphadata',~isnan(residual01(:,:,1)));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Grayscale[-]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(3,2,2);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Av residual Image: cg 3, it end','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,residual01(:,:,2),'alphadata',~isnan(residual01(:,:,2)));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Grayscale[-]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(3,2,3);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Av residual Image: cg 2, it 1','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,residual02(:,:,1),'alphadata',~isnan(residual02(:,:,1)));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Grayscale[-]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(3,2,4);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Av residual Image: cg 2, it end','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,residual02(:,:,2),'alphadata',~isnan(residual02(:,:,2)));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Grayscale[-]';
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(3,2,5);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Av residual Image: cg 1, it 1','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,residual03(:,:,1),'alphadata',~isnan(residual03(:,:,1)));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Grayscale[-]';
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(3,2,6);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Av residual Image: cg 1, it end','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,residual03(:,:,2),'alphadata',~isnan(residual03(:,:,2)));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='Grayscale[-]';
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
set(ax1,'Position',get(ax2, 'position'));

%% Incremental strain calculation 
% ========================

dx = mean(diff(x));
dy = mean(diff(y));
[F.xx, F.xy] = gradient(Ux,dx,dy);
[F.yx, F.yy] = gradient(Uy,dx,dy);

F.xx(isnan(F.xx)) = 0;
F.xy(isnan(F.xy)) = 0;
F.yx(isnan(F.yx)) = 0;
F.yy(isnan(F.yy)) = 0;
[n, m] = size(F.xx);

ZERO = zeros(n,m,'int8');


F.xx = F.xx + 1;
F.yy = F.yy + 1;
F.xx(gdic.mask)=NaN;
F.xy(gdic.mask)=NaN;
F.yx(gdic.mask)=NaN;
F.yy(gdic.mask)=NaN;
% Cauchy Green deformation tensor
C.xx = F.xx .* F.xx + F.yx .* F.yx;
C.xy = F.xx .* F.xy + F.yx .* F.yy;
C.yx = F.yy .* F.yx + F.xy .* F.xx;
C.yy = F.yy .* F.yy + F.xy .* F.xy;

% Finger deformation tensor
B.xx = F.xx .* F.xx + F.xy .* F.xy;
B.xy = F.xx .* F.yx + F.xy .* F.yy;
B.yx = F.yx .* F.xx + F.yy .* F.xy;
B.yy = F.yx .* F.yx + F.yy .* F.yy;

if any(strcmpi(options.strain,{'none'}))
        % Do not compute strain
        E.xx = ZERO;
        E.yy = ZERO;
        E.xy = ZERO;
        E.yx = ZERO;
elseif any(strcmpi(options.strain,{'membrane'}))
    
    % assuming flat initial situation
    [X, Y] = meshgrid(D.cor(inc).xroi,D.cor(inc).yroi);
    [n, m] = size(X);
    Z = zeros(n,m);
    
    [dXx, dXy] = gradient(X+Ux,dx,dy);
    [dYx, dYy] = gradient(Y+Uy,dx,dy);
    [dZx, dZy] = gradient(Z+Uz,dx,dy);

    % this strain definition works for stretched membranes, but is far from
    % universal, test, check, and verify, before using.
    E.xx  = hypot(dXx,dZx) - 1;
    E.yy  = hypot(dYy,dZy) - 1;
    E.xy = hypot(dXy,dZx);
    E.yx = hypot(dYx,dZy);

elseif any(strcmpi(options.strain,{'small'}))
    E.xx = F.xx - 1;
    E.yy = F.yy - 1;
    E.xy = 0.5*(F.xy + F.yx);
    E.yx = 0.5*(F.xy + F.yx);
    
elseif any(strcmpi(options.strain,{'logarithmic'}))
    C = {C.xx,C.xy;C.yx,C.yy};
        [v, d] = eig2d(C);
    E.xx = log(sqrt(d{1})) .* v{1,1} .* v{1,1} ...
        + log(sqrt(d{2})) .* v{2,1} .* v{2,1} ;
    E.xy = log(sqrt(d{1})) .* v{1,1} .* v{1,2} ...
        + log(sqrt(d{2})) .* v{2,1} .* v{2,2} ;
    E.yx = log(sqrt(d{1})) .* v{1,2} .* v{1,1} ...
        + log(sqrt(d{2})) .* v{2,2} .* v{2,1} ;
    E.yy = log(sqrt(d{1})) .* v{1,2} .* v{1,2} ...
        + log(sqrt(d{2})) .* v{2,2} .* v{2,2} ;
elseif any(strcmpi(options.strain,{'greenlagrange'}))
    % Green Lagrange Strain tensor
    E.xx = 0.5 * (C.xx - 1);
    E.xy = 0.5 * (C.xy);
    E.yx = 0.5 * (C.yx);
    E.yy = 0.5 * (C.yy - 1);
elseif any(strcmpi(options.strain,{'euler-almansi'}))
    width = {B.xx,B.xy;B.yx,B.yy};
    [v, d] = eig2d(width);
    E.xx = (1./d{1}) .* v{1,1} .* v{1,1} ...
        + (1./d{2}) .* v{2,1} .* v{2,1} ;
    E.xy = (1./d{1}) .* v{1,1} .* v{1,2} ...
        + (1./d{2}) .* v{2,1} .* v{2,2} ;
    E.yx = (1./d{1}) .* v{1,2} .* v{1,1} ...
        + (1./d{2}) .* v{2,2} .* v{2,1} ;
    E.yy = (1./d{1}) .* v{1,2} .* v{1,2} ...
        + (1./d{2}) .* v{2,2} .* v{2,2} ;

    E.xx = 0.5*(1-E.xx);
    E.xy = 0.5*(0-E.xy);
    E.yx = 0.5*(0-E.yx);
    E.yy = 0.5*(1-E.yy);
end

figure(7)
ax1=subplot(2,2,1);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Incremental strain Field \epsilon_x','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,F.xx,'alphadata',~isnan(Ux));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));

ax1=subplot(2,2,2);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Incremental strain Field \epsilon_y','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,F.yy,'alphadata',~isnan(E.yy));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(2,2,3);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Incremental strain Field \epsilon_{xy}','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,F.xy,'alphadata',~isnan(Ux));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));

ax1=subplot(2,2,4);
imagesc(x1,y1,crs.Manf(:,:,1));
colormap(ax1,gray(256));
title('Incremental strain Field \epsilon_{yx}','Fontsize',fontsize);
xlabel('x [\mum]');
ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,F.yx,'alphadata',~isnan(E.yy));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));
if printpdf
    savepdf('dic_StrainFields')
end
