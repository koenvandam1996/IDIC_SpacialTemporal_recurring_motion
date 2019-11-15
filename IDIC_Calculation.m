function gdic = IDIC_Calculation(x,y,u,phi,Manifold,options)
currentversion = 'version 1(with working Mask and Manifold to obtain one average field)';
tic % start a cpu time counter
NoI  = length(Manifold(1,1,:))-1;

% perform some input checking
options = inputcheckfun(Manifold(:,:,1),x,y,u,phi,options,currentversion,NoI);

% reduce memory footprint with singles
x = single(x); y = single(y); Manifold=single(Manifold);

background=0;
if options.background>0
    [n, m] = size(Manifold(:,:,1));
    background  = zeros(n,m,'single');
    for i=1:length(Manifold(1,1,:))
        background=background+Manifold(:,:,i);
    end
    background=background/length(Manifold(1,1,:));
    background=imgaussfilt(background,options.background);
    for i=1:length(Manifold(1,1,:))
        l=Manifold(:,:,i)-background;
        scaled = l- min(min(l));  
        scaled = scaled / max(max(scaled));
        Manifold(:,:,i)=scaled*255;
    end
end

% open the logfile
if ischar(options.logfile)
    fid = fopen(options.logfile,'a');
    fprintf(fid,'===================================== \n');
    fprintf(fid,'globalDIC %s\n',currentversion);
    if iscell(options.comment)
        for k = 1:length(options.comment)
            fprintf(fid,'%s \n',options.comment{k});
        end
    end
    fprintf(fid,'opening log file: %s \n',datestr(now));
    fprintf(fid,'===== Correlating =================== \n');
else
    fid = -1;
end

if options.verbose > 0
    fprintf('globalDIC %s\n',currentversion);
    fprintf('===================================== \n');
end

% get the number of coarse graining steps (see subfunction)
Ncoarse = coarsestepsfun(Manifold(:,:,1),u,options);
options.Ncoarse = Ncoarse;

if options.list
    coarsesteps = options.coarsesteps;
    Ncoarse = length(coarsesteps);
else
    coarsesteps = Ncoarse:-1:1;
end

% get the number of degrees of freedom
Ndof = size(phi,1);
factor=options.factor;
% for each crs graining step
for kc = 1:Ncoarse
    % store initial guess
    iguess = u;
    
    if options.list
        maxiter = options.maxiter(kc);
        convcrit = options.convcrit(kc);
    else
        maxiter = options.maxiter(1);
        convcrit = options.convcrit(1);
    end
    
    sps = 2^(coarsesteps(kc)-1);
    
    try % skip a coarse graining step if it gives some error
        
        % ============================
        % Coarse Graining
        % ============================
        if options.verbose > 0
            % output to command window
            fprintf('crs:%d/%d, sps:%2d ',kc,Ncoarse,sps)
        end
        if fid ~= -1
            % output to log file
            fprintf(fid,'crs:%d/%d, sps:%2d ==================== \n',kc,Ncoarse,sps);
        end
        crs = coarsegrainfun(x,y,coarsesteps(kc),options,Manifold);

        % ============================
        % Allocate Memory
        % ============================
        Ux = zeros(crs.Npx,options.datatype);
        Uy = zeros(crs.Npx,options.datatype);
        Uz = zeros(crs.Npx,options.datatype);
        
        % progress indicator (1) memory allocated
        if options.verbose > 0; fprintf('='); end
        % ============================
        % update h for initial guess average displacement field
        % ============================
        % loop over the shapefunctions and store the basis function fields              
        Phii=zeros([crs.Npx Ndof],options.datatype);
        for kdof = 1:Ndof
            % call the matlab function which creates the basis for this DOF
            [Phi, dir] = basis(crs.x,crs.y,phi,kdof,options);
            % Store full shape function field
            Phii(:,:,kdof)=Phi;
            % Calculate Displacement Fields
            % ============================
            % separate DOF per direction
            if dir == 1
                Ux  = Ux + u(kdof) * Phi ;
            elseif dir == 2
                Uy  = Uy + u(kdof) * Phi ;
            elseif dir == 3
                Uz  = Uz + u(kdof) * Phi ;
            else
                error('unexpected value for the direction of phi')
            end
        end
        
        % progress indicator (2) displacement field interpolated
        if options.verbose > 0; fprintf('='); end
        
        % convert deformed images with the borders
        for i=1:NoI
            crs.Manh(:,:,i)= interp2(crs.Xfov,crs.Yfov,crs.Mang(:,:,i+1),crs.X+(Ux*factor(i)),crs.Y+(Uy*factor(i)),options.interpmeth,0);
            crs.Manh(:,:,i)=crs.Manh(:,:,i)-Uz;
        end
        
        % progress indicator (3) lsq mass matrix build
        if options.verbose > 0; fprintf('='); end
        % ============================
        % Iterating
        % ============================
        % progress indicator (full) mass matrix build
        if options.verbose > 0
            % output to command window
            fprintf('================== \n')
        end
        
        % initiate some (large) residual and iteration counters
        r       = 1e3;
        conv    = 0;
        it      = 0;
        
        if options.debug
            clear Dbug
        end
        
        % iterating
        while ~conv
            it = it + 1;
            % ============================
            % Build M
            % ============================
            M  = zeros((Ndof+NoI-1),(Ndof+NoI-1),options.datatype);
            b  = zeros((Ndof+NoI-1),1,options.datatype);
            Mc=0;
            % calculate M and b by looping over the image pairs
            for i=1:NoI
                r       =  crs.Manf(:,:,i) - crs.Manh(:,:,i);
                % update m and M
                [m, mc] = buildm(Ndof,phi,options,Ux,Uy,Phii,factor,r,crs.Manf(:,:,i),crs.Manh(:,:,i),i,NoI,crs.Npx,crs.pixelsize(1),crs.pixelsize(2),crs.mask,kc);
                M = M + (m(crs.unmask,:).' * m(crs.unmask,:));
                Mc= Mc + mc;
                % construct the righthand member
                b = b + (m(crs.unmask,:)' * r(crs.unmask));
            end
            % solve for change in degrees of freedom
            du = (M+Mc)\b;
            u  = u + du(1:Ndof);
            % store new factors with first factor fixed at 1
            for i=2:NoI
                factor(i)=factor(i)+(du(Ndof+i-1));
            end
            
            if options.verbose > 0
                % output to command window
                fprintf('it: %2d, ',it);
            end
            
            % update displacement fields
            % ==========================
            % reset to zero
            % ============================
            Ux = zeros(crs.Npx,options.datatype);
            Uy = zeros(crs.Npx,options.datatype);
            Uz = zeros(crs.Npx,options.datatype);
            
            for kdof = 1:Ndof
                % Calculate Displacement Fields
                dir=phi(kdof,3);
                % ============================
                % separate DOF per direction
                if dir == 1
                    Ux  = Ux + u(kdof) * Phii(:,:,kdof) ;
                elseif dir == 2
                    Uy  = Uy + u(kdof) * Phii(:,:,kdof) ;
                elseif dir == 3
                    Uz  = Uz + u(kdof) * Phii(:,:,kdof) ;
                else
                    error('unexpected value for the direction of phi')
                end
            end

            % convert deformed images with the borders
            for i=1:NoI
                crs.Manh(:,:,i)= interp2(crs.Xfov,crs.Yfov,crs.Mang(:,:,i+1),crs.X+(Ux*factor(i)),crs.Y+(Uy*factor(i)),options.interpmeth,0);
                crs.Manh(:,:,i)=crs.Manh(:,:,i)-Uz;
            end
            
            if options.debug
                % store some data for debugging
                Dbug(it).crs  = crs;
                assignin('base',sprintf('Dbug%02d',kc),Dbug)
            end
            r=0;
            for i=1:NoI
                r       = r + crs.Manf(:,:,i) - crs.Manh(:,:,i);
            end
            r       = r/NoI;
            r(crs.mask)=NaN;
            meanr   = mean(abs( r(crs.unmask) ));
            meandu  = mean(mean(abs(du))) / crs.sps;
            maxdu   = max(max(abs(du))) / crs.sps;
            meanalpha = mean(mean(abs(du(Ndof+1:Ndof+NoI-1)))) / crs.sps;
                        
            if options.verbose > 0
                % output to command window
                fprintf('mean(r):%10.3e,  mean(du):%9.2e, mean(alpha):%9.2e,  max(du):%9.2e\n',meanr,meandu, meanalpha, maxdu)
            end
            if fid ~= -1
                % output to log file
                fprintf(fid,'    it: %2d, mean(r):%10.3e,  mean(du):%9.2e,  max(du):%9.2e \n',it,meanr,meandu,maxdu);
            end
            
            if options.debug
                % store some data for debugging
                Dbug(it).du     = du;
                Dbug(it).u      = u;
                Dbug(it).r      = r;
                Dbug(it).factor = factor;
                Dbug(it).meanr  = meanr;
                Dbug(it).meandu = meandu;
                Dbug(it).maxdu  = maxdu;
                Dbug(it).crs    = crs;
                assignin('base',sprintf('Dbug%02d',kc),Dbug)
                 assignin('base','crs',crs)
            end
            
            % test for convergence
            conv = meandu < convcrit;
            if it==1
                residual(:,:,1)=r;
            elseif conv || it==maxiter
                residual(:,:,2)=r;
                assignin('base',sprintf('residual%02d',kc),residual)
                clear residual
                assignin('base','crs',crs)
            end
            
            % test if maximum iterations is reached
            if it >= maxiter
                if options.verbose > 0
                    % output to command window
                    fprintf('    maximum iterations reached\n')
                end
                if fid ~= -1
                    % output to log file
                    fprintf(fid,'    maximum iterations reached\n');
                end
                % terminate the while loop
                break
            end
        end
        
        if options.plot
            % Plot Position Fields
            plotop.name     = sprintf('Position Fields (sps: %02d)',crs.sps);
            plotop.titles     = {'f','h','r'};
            plotop.colorlabel = {'z [\mum]'};
            globalDICplot(x,y,crs.Manf(:,:,1),crs.Manh(:,:,1),r,plotop);
            
            % Plot Displacments
            plotop.name     = sprintf('Displacements (sps: %02d)',crs.sps);
            plotop.titles     = {'Ux','Uy','Uz'};
            plotop.colorlabel = {'U_i [\mum]'};
            globalDICplot(x,y,Ux,Uy,Uz,plotop);
        end
        
    catch ME
        % if an error occured in this crs step, reset u, start the next
        if options.verbose > 0
            % output to command window
            getReport(ME)
            fprintf('    coarsegraining step failed, resetting u \n')
        end
        if fid ~= -1
            % output to log file
            fprintf(fid,'    coarsegraining step failed, resetting u \n');
        end
        % store the error in the workspace
        assignin('base','ME',ME)
        % reset u
        u = iguess;
    end %end try catch
end % forloop over crs graining steps

if options.verbose > 0
    % output to command window
    fprintf('===== Finished in %8.2e sec ====== \n',toc)
end
if fid ~= -1
    % output to log file
    fprintf(fid,'===== Finished in %8.2e sec ====== \n',toc);
end

if exist('r','var')
    % Genererate output structure
    gdic.u   = u;
    gdic.phi = phi;
    gdic.x   = crs.x;
    gdic.y   = crs.y;
    gdic.Ux  = Ux;
    gdic.Uy  = Uy;
    gdic.Uz  = Uz;
    gdic.r   = r;
    gdic.background   = background;
    gdic.mask = crs.mask;
    gdic.factors=factor;
    % Extra output variables
    if options.fulloutput
        gdic.f    = crs.Manf;
        gdic.h    = crs.Manh;
        gdic.g    = crs.g;
        gdic.xfov = crs.xfov;
        gdic.yfov = crs.yfov;
        gdic.Ix   = crs.indx;
        gdic.Iy   = crs.indy;
        gdic.unmask = crs.unmask;
    end
end

gdic.options = options;
gdic.cputime = toc;

if options.verbose > 1
    % output to command window
    logresults(1,gdic)
end
if fid ~= -1
    % output to log file
    logresults(fid,gdic)
    fprintf(fid,'closing log file: %s \n',datestr(now));
    fprintf(fid,'\n\n\n');
    % close the log file
    fclose(fid);
end
end
% ====================================================
%      end of main function
% ====================================================
% ====================================================
%% begin subfunctions
function [m, mc] = buildm(Ndof,phi,options,Ux,Uy,Phi,factor,r,f,g,i,NoI,Npx,dx,dy,mask,kc)
% This function builds the least squares form of the "stiffness" matrix M
m  = zeros(prod(Npx),(Ndof+NoI-1),options.datatype);
mc = zeros((Ndof+NoI-1),(Ndof+NoI-1),options.datatype);

% Gradient of f
% =========================
% Derivative of f
[dfdx, dfdy] = gradient(f,dx,dy);
% Derivative of g
[dhdx, dhdy] = gradient(g,dx,dy);
if any(strcmp(options.gradg,{'G1'}))
    G.x = dfdx(:);
    G.y = dfdy(:);
elseif any(strcmp(options.gradg ,{'G2'}))
    % deformation gradient tensor
    F = defgrad2d(factor(i)*Ux,factor(i)*Uy,dx,dy);
    % inverse
    Fi  = inverse2d(F);
    % transpose
    FiT = transpose2d(Fi);
    % dot product
    G.x = 0.5*FiT{1,1}.*(dhdx) + 0.5*FiT{1,2}.*(dhdy);
    G.y = 0.5*FiT{2,1}.*(dhdx) + 0.5*FiT{2,2}.*(dhdy);
    G.x=G.x(:);
    G.y=G.y(:);
elseif any(strcmp(options.gradg ,{'G3'}))
    % deformation gradient tensor
    F = defgrad2d(factor(i)*Ux,factor(i)*Uy,dx,dy);
    % inverse
    Fi  = inverse2d(F);
    % transpose
    FiT = transpose2d(Fi);
    % dot product
    G.x = 0.5*FiT{1,1}.*(dfdx) + 0.5*FiT{1,2}.*(dfdy);
    G.y = 0.5*FiT{2,1}.*(dfdx) + 0.5*FiT{2,2}.*(dfdy);
    G.x=G.x(:);
    G.y=G.y(:);
elseif any(strcmp(options.gradg ,{'G4'}))
    G.x = dhdx(:);
    G.y = dhdy(:);
elseif any(strcmp(options.gradg,{'G5'}))
    G.x = dfdx(:);
    G.y = dfdy(:);
elseif any(strcmp(options.gradg ,{'G6'}))
    % deformation gradient tensor
    F = defgrad2d(factor(i)*Ux,factor(i)*Uy,dx,dy);
    % inverse
    Fi  = inverse2d(F);
    % transpose
    FiT = transpose2d(Fi);
    % dot product
    G.x = 0.5*FiT{1,1}.*(dfdx+dhdx) + 0.5*FiT{1,2}.*(dfdy+dhdy);
    G.y = 0.5*FiT{2,1}.*(dfdx+dhdx) + 0.5*FiT{2,2}.*(dfdy+dhdy);
    G.x=G.x(:);
    G.y=G.y(:);
elseif any(strcmp(options.gradg ,{'G7'}))
    G.x = 0.5*(dfdx(:) + dhdx(:));
    G.y = 0.5*(dfdy(:) + dhdy(:));
end
G.z = ones(prod(Npx),1);
G.x(mask)=0;
G.y(mask)=0;
Ux(mask)=0;
Uy(mask)=0;
r(mask)=0;
% ============================
% Build m
% ============================
% loop over the shapefunctions (building the stiffness matrix)
for kdof = 1:Ndof
    % ============================
    % Basis Functions
    % ============================
    % call the matlab function which creates the basis for this DOF
    dir = phi(kdof,3);
    %[Phi, dir] = basis(crs.x,crs.y,phi,kdof,options);
    phii=Phi(:,:,kdof);
    % least squares form of the mass matrix (M = m' * m)
    % build one column of m
    if dir == 1
        Gm = G.x .* phii(:);
    elseif dir == 2
        Gm = G.y .* phii(:);
    elseif dir == 3
        Gm = - G.z .* phii(:);
    end
    % store the column to the full matrix in memory
    m( : , kdof ) = m( : , kdof )+Gm*factor(i);
    if kc>=options.mc
        if i~=1
        mc( Ndof+i-1 , kdof ) = Gm'*r(:);
        mc( kdof, Ndof+i-1  ) = mc( Ndof+i-1 , kdof );
        end
    end
end
if i~=1
    m(:,Ndof+i-1)     = ((G.x .* Ux(:))+ (G.y .* Uy(:)));
end % end of loop over Ndof
end
% =======================================================================
function crs = coarsegrainfun(x,y,coarsegrainstep,options,Manifold)
% this function coarse grains the images f and g and crops f to the ROI

%Ncoarse = options.Ncoarse;
[n, m]   = size(Manifold(:,:,1));

% crs graining factor
crs.factor = coarsegrainstep - 1;
% superpixelsize
crs.sps = 2^(crs.factor);

% Building the mask image
% ============================
Imask  = options.mask;
Imask  = [Imask ; find(isnan(Manifold(:,:,1))) ] ;

mask   = zeros(n,m);
mask(Imask) = 1 ;

% Crop to integer superpixel size
% =================================
nmod = mod(n,crs.sps);
In = floor(1+nmod/2):(n-(nmod-floor(nmod/2)));

mmod = mod(m,crs.sps);
Im = floor(1+mmod/2):(m-(mmod-floor(mmod/2)));

% crop
Manifold=Manifold(In,Im,:);
y = y(In);
x = x(Im);
mask = mask(In,Im);

% new size of the image
[n, m] = size(Manifold(:,:,1));

% number of superpixels
if coarsegrainstep ~= 1
    crs.Npx = [n m]./crs.sps ;
else
    crs.Npx = [n m] ;
end

% ============================
% Coarse Graining
% ============================

if coarsegrainstep ~= 1
    % left coarse graining permutation matrix
    % (unity matrix with each row duplicated superpixel size times)
    Pn = repmat(eye(crs.Npx(1)),crs.sps,1);
    Pn = reshape(Pn,crs.Npx(1),n)';
    
    % right coarse graining permutation matrix
    Pm = repmat(eye(crs.Npx(2)),crs.sps,1);
    Pm = reshape(Pm,crs.Npx(2),m)';
    
    % crs grain the mask
    crs.m = ( Pn' * mask * Pm ) ./ crs.sps^2;
    
    % crs grain x and y
    crs.x = ( x * Pm ) ./ crs.sps;
    crs.y = ( y * Pn ) ./ crs.sps;
    
    if any(isnan(Manifold))
        error('there are NaNs in your image')
    end
    % crs grain Manifold 
    for i=1:length(Manifold(1,1,:))
        crs.Mang(:,:,i)= ( Pn' * Manifold(:,:,i) * Pm ) ./ crs.sps^2;
    end
    
else % highest detail
    crs.Mang= Manifold;
    crs.x = x;
    crs.y = y;
    crs.m = mask ;
end

% ============================
% Select the ROI in f
% ============================
% get the size of a superpixel (in mu)
pixelsize(1) = mean(diff(x));
pixelsize(2) = mean(diff(y));
crs.pixelsize = pixelsize*crs.sps;

% the roi needs to contain whole superpixels
ROI(1) = options.ROI(1) - 0.5*pixelsize(1) + 0.5*crs.pixelsize(1);
ROI(2) = options.ROI(2) + 0.5*pixelsize(1) - 0.5*crs.pixelsize(1);
ROI(3) = options.ROI(3) - 0.5*pixelsize(2) + 0.5*crs.pixelsize(2);
ROI(4) = options.ROI(4) + 0.5*pixelsize(2) - 0.5*crs.pixelsize(2);

% the indices which are in the ROI
crs.indx = find(crs.x >= ROI(1) & crs.x <= ROI(2));
crs.indy = find(crs.y >= ROI(3) & crs.y <= ROI(4));

% storing the FOV position vectors
crs.xfov = crs.x;
crs.yfov = crs.y;

% create FOV position fields
[crs.Xfov, crs.Yfov] = meshgrid(crs.xfov,crs.yfov);

% cropping f to the ROI
crs.x = crs.x(crs.indx);
crs.y = crs.y(crs.indy);
for i=1:length(Manifold(1,1,:))
    crs.Manf(:,:,i)= crs.Mang(crs.indy,crs.indx,i);
end

% store the mask as binary image
crs.m      = crs.m(crs.indy,crs.indx);
% store the mask as list of indices
crs.mask   = find(crs.m > 0);
% store the inverse mask index list
crs.unmask = find(crs.m == 0);

% the number of pixels in the ROI
crs.Npx = size(crs.Manf(:,:,1));

% create position fields
[crs.X, crs.Y] = meshgrid(crs.x,crs.y);

if options.debug
    assignin('base','crs',crs)
end
end
% =======================================================================
function [Phi, dir] = basis(x,y,phi,kdof,options)
% this function calls the correct basis-function type build-function, and
% returns the field for this basis-function and direction
% normalize the domain
m = length(x);
n = length(y);
if options.normalized
    x = linspace(-1,1,m);
    y = linspace(-1,1,n);
end
% store the direction in which the basis operates
dir = phi(kdof,3);
% store the order
a = phi(kdof,1);
b = phi(kdof,2);
% call the correct basis-function
if any(strcmpi(options.basis,{'polynomial','p'}))
    Phi = basis_polynomial(x,y,a,b);
elseif any(strcmpi(options.basis,{'chebyshev','c'}))
    Phi = basis_chebyshev(n,m,a,b);
elseif any(strcmpi(options.basis,{'bspline','b'}))
    % for splines also get the number of (unique) knots and the order
    na = phi(kdof,4);
    nb = phi(kdof,5);
    np = phi(kdof,6);
    Phi = basis_bspline(n,m,a,b,na,nb,np);
else
    error('unknown basis: [%s]',options.basis);
end
end
% =======================================================================
function Phi = basis_polynomial(x,y,a,b)
% This function creates a polynomial basisfunction
% which is defined at the pixel level
% this basis creates a two dimensional Nth order field of the type
% Phi = x^a * y^b;
% create the shape in each direction
x = x.^a;
y = y.^b;
% expand each direction to a field
[X, Y] = meshgrid(x,y);
% multiply the two
Phi = X .* Y;
end
% =======================================================================
function Phi = basis_chebyshev(n,m,a,b)
% This function creates a chebyshev basisfunction of the first kind
% which is defined at the pixel level. 
% http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html

% chebyshev functions only make sense on a normalized base
x = linspace(-1,1,m);
y = linspace(-1,1,n);

% For the x direction
% =====================

% initial the first two chebyshev polynomials
Tn = ones(1,m);
Tk = x;

% Recursively build the polynomials from the first two
if a == 0
    Ta = Tn;
elseif a == 1
    Ta = Tk;
else
    for k = 2:a
        
        % calculate the new polynomial
        Ta = 2*x.*Tk - Tn;
        
        % update the previous two
        Tn = Tk;
        Tk = Ta;
    end
end

% For the y direction
% =====================

% initial the first two chebyshev polynomials
Tn = ones(1,n);
Tk = y;

% Recursively build the polynomials from the first two
if b == 0
    Tb = Tn;
elseif b == 1
    Tb = Tk;
else
    for k = 2:b
        
        % calculate the new polynomial
        Tb = 2*y.*Tk - Tn;
        
        % update the previous two
        Tn = Tk;
        Tk = Tb;
    end
end

% expand each direction to a field
[X, Y] = meshgrid(Ta,Tb);

% combine the direction to one basis-function
Phi = X.*Y ;
end
% =======================================================================
function Phi = basis_bspline(n,m,a,b,na,nb,np)
% nk = number of (unique) knots
% np = polynomial degree

% forloop for x and y directions
for xy = 1:2
    % space the knots uniformly
    if xy == 1
        nk = na;
    else
        nk = nb;
    end
    knots = linspace(-1,1,nk);
    
    % Make the B-Spline Open: repeat end knots p times
    knots = [knots(1)*ones(1,np) knots knots(end)*ones(1,np)];
    % new number of knots
    Ni = length(knots);
    
    % Initiate the parametric space
    if xy == 1
        Nzeta = m;
    else
        Nzeta = n;
    end
    zeta = linspace(knots(1),knots(end),Nzeta);
    
    % Zero order
    % =============
    
    % degree
    p = 0;
    
    % number of basis-functions
    Nn = Ni-(p+1);
    
    % initiate matrix for basis functions
    N0 = zeros(Nn,Nzeta);
    
    for i = 1:Nn
        if i == Nn-np
            % only for the right most (single) knot
            I = find(zeta >= knots(i) & zeta <= knots(i+1));
        else
            % for all other knots
            I = find(zeta >= knots(i) & zeta < knots(i+1));
        end
        
        % set the function to zero
        N0(i,I) = 1;
    end
    % copy the zero basis for later use
    N = N0;
    
    % Subsequent orders
    % =============
    for p = 1:np
        % calculate the number of shape functions for this degree
        Nn = Ni-(p+1);
        
        % store the previous N as N1
        N1 = N;
        % initiate the current N
        N  = zeros(Nn,Nzeta);
        
        % forloop over the knots
        for i = 1:Nn
            if knots(i+p) == knots(i)
                % if first term == zero (double knot on right side)
                N(i,:) = N1(i+1,:).* (knots(i+p+1) - zeta ) ./ ( knots(i+p+1) - knots(i+1) );
            elseif knots(i+p+1) == knots(i+1)
                % if second term is zero (double knot on left side)
                N(i,:) = N1(i,:)  .* (   zeta - knots(i)   ) ./ (  knots(i+p)   - knots(i)  ) ;
            else
                % for all other knots
                N(i,:) = N1(i,:)  .* (   zeta - knots(i)   ) ./ (  knots(i+p)   - knots(i)  ) + ...
                    N1(i+1,:).* (knots(i+p+1) - zeta ) ./ ( knots(i+p+1) - knots(i+1) );
            end
        end
        
    end
    
    % Select the combination of shapes
    if xy == 1
        Nx = N(a,:);
    else
        Ny = N(b,:);
    end
    
end

% expand each direction to a field
[X, Y] = meshgrid(Nx,Ny);

% combine the direction to one basis-function
Phi = X.*Y ;
end
% =======================================================================
function Ncoarse = coarsestepsfun(f,u,options)
% this function returns the number of coarse graining steps
% Ndofget size of f
N = length(u);
[n, m] = size(f);

% direction with the least pixels
p = min([n m]);

% coarse image should at least have more pixels then the number of DOF
pmin = max([8 2^nextpow2(N)]);

% the maximum number of coarse graining steps (including full detail)
Ncoarse = floor(log2(p/pmin)) ;

if options.coarsesteps(1) < 1
    Ncoarse = Ncoarse + options.coarsesteps(1);
    Ncoarse = max([1 Ncoarse]);
else
    Ncoarse = options.coarsesteps(1);
end
end
% =======================================================================
function options = inputcheckfun(f,x,y,u,phi,options,currentversion,NoI)
% this function handles all the input checking and sets the default options
% if they are not provided

% check
if length(size(f)) ~= 2
    error('f must be a nxm matrix')
end
if size(f,2) ~= length(x)
    error('x must be a 1xm matrix')
end
if size(f,1) ~= length(y)
    error('y must be a 1xn matrix')
end
if length(u) ~= size(phi,1)
    error('u must have the same length as phi')
end

% ======================
% default settings
% ======================
defaultoptions.ROI         = 'default';
defaultoptions.list        = false ;
defaultoptions.coarsesteps = 0 ;
defaultoptions.convcrit    = 1e-5 ;
defaultoptions.maxiter     = 20 ;
defaultoptions.interpmeth  = 'spline';
defaultoptions.normalized  = true;
defaultoptions.basis       = 'polynomial';
defaultoptions.gradg       = 'G7' ;
defaultoptions.verbose     = 1 ;
defaultoptions.fulloutput  = false;
defaultoptions.plot        = false ;
defaultoptions.debug       = false ;
defaultoptions.datatype    = 'double';
defaultoptions.logfile     = false ;
defaultoptions.comment     = false ;
defaultoptions.mask        = [] ;
defaultoptions.mc          = false ;
defaultoptions.brightness  = false ;
defaultoptions.background  = 0;
defaultoptions.factor      = ones(NoI,1);
% ======================
% Get the input options
% ======================
tempoptions = defaultoptions;
% get the default options list
field = fieldnames(defaultoptions);
% for each field in the default options
for k = 1:length(field)
    if isfield(options,field{k})
        % if this field also exists in the input options --> overwrite
        tempoptions.(field{k}) = options.(field{k});
    end
end
% Set all options
options = tempoptions;
% ======================
% Add u to the options for reference
% ======================
options.init = u;
options.version = currentversion;
% ======================
% Input checking
% ======================
% ROIcheck input
% ======================
if length(options.ROI)~=4
    error('unknown ROI option [%s]',options.ROI)
end
% List
% ======================
if ~islogical(options.list)
    error('unknown list option [%s]',options.list)
end
if options.list && any(options.coarsesteps == 0)
    error('the list option is incompatible with coarsesteps = 0')
end

% Coarse steps
% ======================
if ischar(options.coarsesteps(1))
    error('unknown coarsesteps option [%s]',options.coarsesteps(1))
elseif ~isnumeric(options.coarsesteps(1)) && (mod(options.coarsesteps(1),1) ~= 0)
    error('unknown coarsesteps option [%g]',options.coarsesteps(1))
end
if options.list
    Nlist = length(options.coarsesteps);
else
    options.coarsesteps = options.coarsesteps(1);
end

% Convergence criterea
% ======================
if ~isnumeric(options.convcrit(1))
    error('unknown convcrit option [%s]',options.convcrit(1))
end
if options.list && length(options.convcrit) == 1
    options.convcrit = ones(1,Nlist)*options.convcrit;
elseif options.list && length(options.convcrit) ~= Nlist
    error('the length of convcrit must match the length of coarsesteps')
elseif ~options.list
    options.convcrit = options.convcrit(1);
end

% Maximum number of iterations
% ======================
if ~isnumeric(options.maxiter(1)) && (mod(options.maxiter(1),1) ~= 0)
    error('unknown maxiter option [%s]',options.maxiter(1))
end
if options.list && length(options.maxiter) == 1
    options.maxiter = ones(1,Nlist)*options.maxiter;
elseif options.list && length(options.maxiter) ~= Nlist
    error('the length of maxiter must match the length of coarsesteps')
elseif ~options.list
    options.maxiter = options.maxiter(1);
end

% Interpmeth
% ======================
opt = {'nearest','linear','cubic','spline'};
if ~any(strcmp(options.interpmeth,opt))
    error('unknown interpolation method [%s]',options.interpmeth)
end

% Normalized basis
% ======================
if ~islogical(options.normalized)
    error('unknown normalized option [%s]',options.normalized)
end

% Basis function family
% ======================
if ~any(strcmpi(options.basis,{'polynomial','p','chebyshev','c','bspline','b'}))
    error('unknown basis option [%s]',options.basis)
end

% Gradient of g
% ======================
% if ~islogical(options.gradg)
%     error('unknown gradg option [%s]',options.gradg)
% end

% Verbose
% ======================
if mod(options.verbose,1) ~= 0
    error('unknown verbose option [%s]',options.verbose)
end

% Full Output
% ======================
if ~islogical(options.fulloutput)
    error('unknown fulloutput option [%s]',options.fulloutput)
end

% Plot
% ======================
if ~islogical(options.plot)
    error('unknown plot option [%s]',options.plot)
end

% Debug
% ======================
if ~islogical(options.debug)
    error('unknown debug option [%s]',options.debug)
end

% Datatype
% ======================
opt = {'single','double'};
if ~any(strcmpi(options.datatype,opt))
    error('unsupported datatype [%s]',options.datatype)
end

% Use logfile
% ======================
if islogical(options.logfile)
    if options.logfile
        error('options.logfile must be false or a string pointing to a filename')
    end
elseif ~ischar(options.logfile)
    error('unknown logfile option')
end

% Use Comment
% ======================
if islogical(options.comment)
    if options.comment
        error('options.comment must be false or a string (or a cell with strings)')
    end
elseif ischar(options.comment)
    % convert to cell
    options.comment = {options.comment};
elseif ~iscell(options.comment)
    error('unknown comment option')
end
end
% =======================================================================
function logresults(fid,gdic)
% print some results to the logfile (or screen)

% print a summary of the options structure to the logfile
% =====================
fprintf(fid,'\n===== Options ======================= \n');

% get all fieldnames in the gdic structure
field = fieldnames(gdic.options);

% forloop over each field
for k = 1:length(field)
    % get the content of the field
    content = gdic.options.(field{k});
    % convert to string
    if strcmpi(field{k},'mask')
        content = sprintf('%d',length(content));
    elseif strcmp(field{k},'phi')
        content = sprintf('%d',size(content,1));
    elseif strcmp(field{k},'init')
        content = sprintf('%d',size(content,1));
    elseif iscell(content)
        content = 'cell';
    elseif isnumeric(content) || islogical(content)
        content = strtrim(num2str(content));
    end
    fprintf(fid,'%13s: %s \n',field{k},content);
end

fprintf(fid,'\n===== Results ======================= \n');

% print a summary of the gdic structure to the logfile
% =====================
fprintf(fid,'===== gdic structure \n');

% get all fieldnames in the gdic structure
field = fieldnames(gdic);

% forloop over each field
for k = 1:length(field)
    % get the size of the field
    [n, m] = size(gdic.(field{k}));
    % get the class type of the field
    cls   = class(gdic.(field{k}));
    % if it is a structure, don't calculate the mean
    if isstruct(gdic.(field{k}))
        % print to the logfile
        fprintf(fid,'%10s: [%4dx%4d %8s] \n',field{k},n,m,cls);
    else
        % mean of the field
        tmp   = gdic.(field{k})(:) ;
        mn    = mean(tmp(~isnan(tmp))) ;
        % print to the logfile
        fprintf(fid,'%10s: [%4dx%4d %8s], mean: %8.1e \n',field{k},n,m,cls,mn);
    end
end

% print a the result of each DOF to the logfile
% =====================
fprintf(fid,'===== DOF \n');

Ndof = length(gdic.u);
% the output of each DOF
for k = 1:Ndof
    dirs = 'xyz';
    a = gdic.phi(k,1);
    b = gdic.phi(k,2);
    dir = dirs(gdic.phi(k,3));
    % print the DOF
    fprintf(fid,'     u(%2d): %9.2e, dir: %1s, order: [%d,%d] \n',k,gdic.u(k),dir,a,b);
end
fprintf(fid,'====================================== \n');
end
% =======================================================================
function F = defgrad2d(Ux,Uy,dx,dy)
% compute a matrix of 2D deformation gradient tensors from the x
% displacemetn field (Ux) and the y displacement field(Uy)
%
% inputs
% Ux = (n,m);
% Uy = (n,m);
%
% dx,dy the pixelsize in x and y
%
% output
% F = {2,2}(n,m), i.e. a 2x2 cell array of matrices, where
%    F{1,1}(:,:) is F11
%    F{1,2}(:,:) is F12
%    F{2,1}(:,:) is F21
%    F{2,2}(:,:) is F22
%
% the returned matrix has the same structure

% deformation gradient tensor
[F{1,1}, F{1,2}] = gradient(Ux,dx,dy);
[F{2,1}, F{2,2}] = gradient(Uy,dx,dy);
F{1,1} = F{1,1} + 1;
F{2,2} = F{2,2} + 1;
end
% =======================================================================
function Ainv = inverse2d(A)
% invert a matrix of 2D tensors
% A = {2,2}(n,m), i.e. a 2x2 cell array of matrices, where
%    A{1,1}(:,:) is A11
%    A{1,2}(:,:) is A12
%    A{2,1}(:,:) is A21
%    A{2,2}(:,:) is A22
%
% the returned matrix has the same structure

% using: http://en.wikipedia.org/wiki/Invertible_matrix

% determinant
D = A{1,1}.*A{2,2} - A{1,2}.*A{2,1};

% inverse
Ainv{1,1} =  (1./D) .* A{2,2};
Ainv{1,2} = -(1./D) .* A{1,2};
Ainv{2,1} = -(1./D) .* A{2,1};
Ainv{2,2} =  (1./D) .* A{1,1};
end
% =======================================================================
function At = transpose2d(A)
% transpose a matrix of 2D tensors
% A = {2,2}(n,m), i.e. a 2x2 cell array of matrices, where
%    A{1,1}(:,:) is A11
%    A{1,2}(:,:) is A12
%    A{2,1}(:,:) is A21
%    A{2,2}(:,:) is A22
%
% the returned matrix has the same structure
At{1,1} = A{1,1};
At{1,2} = A{2,1};
At{2,1} = A{1,2};
At{2,2} = A{2,2};
end