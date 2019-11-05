function gdic = globalDIC2Dcomb(xback,yback,u,phi,Manifoldback,xfront,yfront,Manifoldfront,xtot,options)
currentversion = 'version 1(with working Mask and Manifold)';
tic % start a cpu time counter
NoI  = length(Manifoldback(1,1,:))-1;

% perform some input checking
options = inputcheckfun(Manifoldback(:,:,1),xback,yback,u,phi,options,currentversion,NoI);

% reduce memory footprint with singles
xback = single(xback); yback = single(yback); Manifoldback=single(Manifoldback);
xfront = single(xfront); yfront = single(yfront); Manifoldfront=single(Manifoldfront);

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
Ncoarse = coarsestepsfun(Manifoldback(:,:,1),u,options);
options.Ncoarse = Ncoarse;

if options.list
    coarsesteps = options.coarsesteps;
    Ncoarse = length(coarsesteps);
else
    coarsesteps = Ncoarse:-1:1;
end

% get the number of degrees of freedom
Ndof = size(phi,1);
factorback=options.factorback;
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
    options.slider(1)=round(options.slider(2)/sps);
    
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
        crsback = coarsegrainfun(xback,yback,coarsesteps(kc),Manifoldback,options.maskback,options.ROIback);
        crsfront = coarsegrainfun(xfront,yfront,coarsesteps(kc),Manifoldfront,options.maskfront,options.ROIfront);
       
        % ============================
        % Allocate Memory
        % ============================
        Ux = zeros([crsback.Npx(1) (crsback.Npx(2)-options.slider(1)+crsfront.Npx(2))],options.datatype);
        Uy = zeros([crsback.Npx(1) (crsback.Npx(2)-options.slider(1)+crsfront.Npx(2))],options.datatype);
        % progress indicator (1) memory allocated
        if options.verbose > 0; fprintf('='); end
        % ============================
        % update h for initial guess average displacement field
        % ============================
        % loop over the shapefunctions
        
        xtot = linspace(1,((crsback.Npx(2)-options.slider(1)+crsfront.Npx(2))),((crsback.Npx(2)-options.slider(1)+crsfront.Npx(2))));
        xtot = xtot - mean(xtot);
        Phii=zeros([crsback.Npx(1) (crsback.Npx(2)-options.slider(1)+crsfront.Npx(2)) Ndof],options.datatype);
        for kdof = 1:Ndof
            % call the matlab function which creates the basis for this DOF
            [Phi, dir] = basis(xtot,crsback.y,phi,kdof,options);
            % Store full shape function field
            Phii(:,:,kdof)=Phi;
            % Calculate Displacement Fields
            % ============================
            % separate DOF per direction
            if dir == 1
                Ux  = Ux + u(kdof) * Phi ;
            elseif dir == 2
                Uy  = Uy + u(kdof) * Phi ;
            else
                error('unexpected value for the direction of phi')
            end
        end
        
        % progress indicator (2) displacement field interpolated
        if options.verbose > 0; fprintf('='); end
        
        % convert deformed images with the borders
        for i=1:NoI
            crsback.Manh(:,:,i)= interp2(crsback.Xfov,crsback.Yfov,crsback.Mang(:,:,i+1),crsback.X+(Ux(:,1:crsback.Npx(2))*factorback(i)),crsback.Y+(Uy(:,1:crsback.Npx(2))*factorback(i)),options.interpmeth,0);
            crsfront.Manh(:,:,i)= interp2(crsfront.Xfov,crsfront.Yfov,crsfront.Mang(:,:,i+1),crsfront.X+(Ux(:,(crsback.Npx(2)-options.slider(1)+1):(crsback.Npx(2)-options.slider(1)+crsfront.Npx(2)))*factorback(NoI+i)),crsfront.Y+(Uy(:,(crsback.Npx(2)-options.slider(1)+1):(crsback.Npx(2)-options.slider(1)+crsfront.Npx(2)))*factorback(NoI+i)),options.interpmeth,0);
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
        unmaskfront=(crsfront.unmask+prod(size(crsback.Manh(:,(1:crsback.Npx(2)),1))));
        unmasktot = vertcat(crsback.unmask,unmaskfront);
        maskfront=(crsfront.mask+prod(size(crsback.Manh(:,(1:crsback.Npx(2)),1))));
        masktot= vertcat(crsback.mask,maskfront);
        if options.debug
            clear Dbug
        end
        
        % iterating
        while ~conv
            it = it + 1;
            % ============================
            % Build M
            % ============================
            M  = zeros((Ndof+(2*NoI)-1),(Ndof+(2*NoI)-1),options.datatype);
            b  = zeros((Ndof+(2*NoI)-1),1,options.datatype);
            Mc=0;
            % construct the average residual
            for i=1:NoI
                r(1:crsback.Npx(1),1:(crsback.Npx(2))) =  crsback.Manf(:,1:(crsback.Npx(2)),i) - crsback.Manh(:,1:(crsback.Npx(2)),i);
                r(1:crsback.Npx(1),(crsback.Npx(2)+1):(crsback.Npx(2)+crsfront.Npx(2))) =  crsfront.Manf(:,1:(crsfront.Npx(2)),i) - crsfront.Manh(:,1:(crsfront.Npx(2)),i);
                for j=1:options.slider(1)
                    r(:,crsback.Npx(2)+j)=r(:,crsback.Npx(2)+j)*((options.slider(1)-j)/options.slider(1));
                    r(:,crsback.Npx(2)-options.slider(1)+j)=r(:,crsback.Npx(2)-options.slider(1)+j)*(j/options.slider(1));
                end
                % update m and M
                [m, mc] = buildm(Ndof,phi,options,Ux,Uy,Phii,factorback,r,crsback.Manf(:,:,i),crsback.Manh(:,:,i),i,NoI,crsback.Npx,crsback.pixelsize(1),crsback.pixelsize(2),masktot,kc,crsfront.Manf(:,:,i),crsfront.Manh(:,:,i),unmasktot,crsback.mask,crsfront.mask);
                M = M + (m(unmasktot,:).' * m(unmasktot,:)); 
                Mc= Mc + mc;
                % construct the righthand member
                b = b + (m(unmasktot,:)' * r(unmasktot));
            end
            % solve for change in degrees of freedom
            du = (M+mc)\b;
            u  = u + du(1:Ndof);
       
            for i=2:NoI*2
                factorback(i)=factorback(i)+(du(Ndof+i-1));
%                 if kc>2 && (factorback(i)>1.5 | factorback(i)<0.55)
%                    factorback(i)=1; 
%                 end
            end
            
            if options.verbose > 0
                % output to command window
                fprintf('it: %2d, ',it);
            end
            
            % update displacement fields
            % ==========================
            % reset to zero
            % ============================
            Ux = zeros([crsback.Npx(1) (crsback.Npx(2)-options.slider(1)+crsfront.Npx(2))],options.datatype);
            Uy = zeros([crsback.Npx(1) (crsback.Npx(2)-options.slider(1)+crsfront.Npx(2))],options.datatype);
        
            for kdof = 1:Ndof
                % Calculate Displacement Fields
                dir=phi(kdof,3);
                % ============================
                % separate DOF per direction
                if dir == 1
                    Ux  = Ux + u(kdof) * Phii(:,:,kdof) ;
                elseif dir == 2
                    Uy  = Uy + u(kdof) * Phii(:,:,kdof) ;
                else
                    error('unexpected value for the direction of phi')
                end
            end

            % convert deformed images with the borders
             for i=1:NoI
            crsback.Manh(:,:,i)= interp2(crsback.Xfov,crsback.Yfov,crsback.Mang(:,:,i+1),crsback.X+(Ux(:,1:crsback.Npx(2))*factorback(i)),crsback.Y+(Uy(:,1:crsback.Npx(2))*factorback(i)),options.interpmeth,0);
            crsfront.Manh(:,:,i)= interp2(crsfront.Xfov,crsfront.Yfov,crsfront.Mang(:,:,i+1),crsfront.X+(Ux(:,(crsback.Npx(2)-options.slider(1)+1):(crsback.Npx(2)-options.slider(1)+crsfront.Npx(2)))*factorback(NoI+i)),crsfront.Y+(Uy(:,(crsback.Npx(2)-options.slider(1)+1):(crsback.Npx(2)-options.slider(1)+crsfront.Npx(2)))*factorback(NoI+i)),options.interpmeth,0);
            end
            
            if options.debug
                % store some data for debugging
                Dbug(it).crs  = crsback;
                assignin('base',sprintf('Dbug%02d',kc),Dbug)
            end
            r = zeros([crsback.Npx(1) (crsback.Npx(2)*2)],options.datatype);
            for i=1:NoI
                r(:,1:(crsback.Npx(2)))       = r(:,1:(crsback.Npx(2))) + crsback.Manf(:,:,i) - crsback.Manh(:,:,i);
                r(:,(crsback.Npx(2)+1):crsback.Npx(2)+crsfront.Npx(2))     = r(:,(crsback.Npx(2)+1):crsback.Npx(2)+crsfront.Npx(2)) + crsfront.Manf(:,:,i) - crsfront.Manh(:,:,i);
            end
            r       = r/NoI;
            %r(crsback.mask)=NaN;
            meanr   = mean(abs( r(crsback.unmask) ));
            meandu  = mean(mean(abs(du))) / crsback.sps;
            maxdu   = max(max(abs(du))) / crsback.sps;
            meanalpha = mean(mean(abs(du(Ndof+1:Ndof+NoI-1)))) / crsback.sps;
                        
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
                Dbug(it).factor = factorback;
                Dbug(it).meanr  = meanr;
                Dbug(it).meandu = meandu;
                Dbug(it).maxdu  = maxdu;
                Dbug(it).crs    = crsback;
                assignin('base',sprintf('Dbug%02d',kc),Dbug)
                 assignin('base','crs',crsback)
            end
            
            % test for convergence
            conv = meandu < convcrit;
            if it==1
                residual(:,:,1)=r;
            elseif conv || it==maxiter
                residual(:,:,2)=r;
                assignin('base',sprintf('residual%02d',kc),residual)
                clear residual
                assignin('base','crs',crsback)
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
            plotop.name     = sprintf('Position Fields (sps: %02d)',crsback.sps);
            plotop.titles     = {'f','h','r'};
            plotop.colorlabel = {'z [\mum]'};
            globalDICplot(xback,yback,crsback.Manf(:,:,1),crsback.Manh(:,:,1),r,plotop);
            
            % Plot Displacments
            plotop.name     = sprintf('Displacements (sps: %02d)',crsback.sps);
            plotop.titles     = {'Ux','Uy'};
            plotop.colorlabel = {'U_i [\mum]'};
            globalDICplot(xback,yback,Ux,Uy,plotop);
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
    masksec=(crsfront.mask+prod(size(crsback.Manh(:,(1:crsback.Npx(2)-options.slider(1))))));
    masktotsec = vertcat(crsback.mask,masksec);
if exist('r','var')
    % Genererate output structure
    gdic.u   = u;
    gdic.phi = phi;
    gdic.x   = crsback.x;
    gdic.y   = crsback.y;
    gdic.Ux  = Ux;
    gdic.Uy  = Uy;
    gdic.r   = r;
    gdic.mask = masktotsec;
    gdic.factors=factorback;
    % Extra output variables
    if options.fulloutput
        gdic.f    = crsback.Manf;
        gdic.h    = crsback.Manh;
        gdic.g    = crsback.g;
        gdic.xfov = crsback.xfov;
        gdic.yfov = crsback.yfov;
        gdic.Ix   = crsback.indx;
        gdic.Iy   = crsback.indy;
        gdic.unmask = crsback.unmask;
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
function [m, mc] = buildm(Ndof,phi,options,Ux,Uy,Phi,factor,r,f,g,i,NoI,Npx,dx,dy,mask,kc,ffront,gfront,unmasktot,maskback,maskfront)
% This function builds the least squares form of the "stiffness" matrix M
Npxtot=size(r);
m  = zeros(prod(Npxtot),(Ndof+(NoI*2)-1),options.datatype);
mc = zeros((Ndof+(NoI*2)-1),(Ndof+(NoI*2)-1),options.datatype);

% Gradient of f
% =========================
% Derivative of f
[dfdx, dfdy] = gradient(f,dx,dy);
[dffrontdx, dffrontdy] = gradient(ffront,dx,dy);
% Derivative of g
[dhdx, dhdy] = gradient(g,dx,dy);
[dhfrontdx, dhfrontdy] = gradient(gfront,dx,dy);
Uxb=Ux(:,1:Npx(2))*factor(i);
Uxfr=Ux(:,Npx(2)-options.slider(1)+1:end)*factor(i+NoI);
Uyb=Uy(:,1:Npx(2))*factor(i);
Uyfr=Uy(:,Npx(2)-options.slider(1)+1:end)*factor(i+NoI);

if any(strcmp(options.gradg,{'G1'}))
    G.x = dfdx(:);
    G.y = dfdy(:);
    G.xfr = dffrontdx(:);
    G.yfr = dffrontdy(:);
elseif any(strcmp(options.gradg ,{'G6'}))
    % deformation gradient tensor
    Fb = defgrad2d(Uxb,Uyb,dx,dy);
    Ffr = defgrad2d(Uxfr,Uyfr,dx,dy);
    % inverse
    Fib  = inverse2d(Fb);
    Fifr  = inverse2d(Ffr);
    % transpose
    FiTb = transpose2d(Fib);
    FiTfr = transpose2d(Fifr);
    % dot product
    G.xb = 0.5*FiTb{1,1}.*(dfdx+dhdx) + 0.5*FiTb{1,2}.*(dfdy+dhdy);
    G.yb = 0.5*FiTb{2,1}.*(dfdx+dhdx) + 0.5*FiTb{2,2}.*(dfdy+dhdy);
    G.xfr = 0.5*FiTfr{1,1}.*(dffrontdx+dhfrontdx) + 0.5*FiTfr{1,2}.*(dffrontdy+dhfrontdy);
    G.yfr = 0.5*FiTfr{2,1}.*(dffrontdx+dhfrontdx) + 0.5*FiTfr{2,2}.*(dffrontdy+dhfrontdy);
    G.xb=G.xb(:);
    G.yb=G.yb(:);
    G.xfr=G.xfr(:);
    G.yfr=G.yfr(:);
elseif any(strcmp(options.gradg ,{'G7'}))
    G.xb = 0.5*(dfdx(:) + dhdx(:));
    G.yb = 0.5*(dfdy(:) + dhdy(:));
    G.xfr = 0.5*(dffrontdx(:) + dhfrontdx(:));
    G.yfr = 0.5*(dffrontdy(:) + dhfrontdy(:));
    G.xb=G.xb(:);
    G.yb=G.yb(:);
    G.xfr=G.xfr(:);
    G.yfr=G.yfr(:);
end
G.z = ones(prod(Npx),1);
G.xb(maskback)=0;
G.yb(maskback)=0;
G.xfr(maskfront)=0;
G.yfr(maskfront)=0;
r(mask)=0;
Uxb(maskback)=0;
Uxfr(maskfront)=0;
Uyb(maskback)=0;
Uyfr(maskfront)=0;

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
    phiib=phii(:,1:Npx(2));
    phiifr=phii(:,Npx(2)-options.slider(1)+1:end);
    % least squares form of the mass matrix (M = m' * m)
    % build one column of m
    if dir == 1
        Gmb = G.xb .* phiib(:);
        Gmfr = G.xfr .* phiifr(:);
    elseif dir == 2
        Gmb = G.yb .* phiib(:);
        Gmfr = G.yfr .* phiifr(:);
    end
    % store the column to the full matrix in memory
    m( 1:prod(Npx) , kdof ) = m( 1:prod(Npx) , kdof )+(Gmb)*factor(i);
    m( prod(Npx)+1:prod(Npxtot) , kdof ) = m(prod(Npx)+1:prod(Npxtot) , kdof )+(Gmfr)*factor(i+NoI);
r=r(:);
if kc>=options.mc
        if i~=1
        mc( Ndof+i-1 , kdof ) = Gmb'*r(( 1:prod(Npx)));
        mc( kdof, Ndof+i-1  ) = mc( Ndof+i-1 , kdof );
        mc( Ndof+NoI+i-1 , kdof ) = Gmfr'*r(( prod(Npx)+1:prod(Npxtot)));
        mc( kdof, Ndof+NoI+i-1  ) = mc( Ndof+i-1 , kdof );
        end
end
end

if i~=1
    m(1:prod(Npx),Ndof+i-1)     = ((G.xb .* Uxb(:))+ (G.yb .* Uyb(:)));
end % end of loop over Ndof
    m(prod(Npx)+1:prod(Npxtot),Ndof+NoI+i-1)     = ((G.xfr .* Uxfr(:))+ (G.yfr .* Uyfr(:)));
end
% =======================================================================
function crs = coarsegrainfun(x,y,coarsegrainstep,Manifold,opmask,oproi)
% this function coarse grains the images f and g and crops f to the ROI

%Ncoarse = options.Ncoarse;
[n, m]   = size(Manifold(:,:,1));

% crs graining factor
crs.factor = coarsegrainstep - 1;
% superpixelsize
crs.sps = 2^(crs.factor);

% Building the mask image
% ============================
Imask  = opmask;
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
ROI(1) = oproi(1) - 0.5*pixelsize(1) + 0.5*crs.pixelsize(1);
ROI(2) = oproi(2) + 0.5*pixelsize(1) - 0.5*crs.pixelsize(1);
ROI(3) = oproi(3) - 0.5*pixelsize(2) + 0.5*crs.pixelsize(2);
ROI(4) = oproi(4) + 0.5*pixelsize(2) - 0.5*crs.pixelsize(2);

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

%if options.debug
%    assignin('base','crs',crs)
%end
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
defaultoptions.ROIback         = 'default';
defaultoptions.ROIfront        = 'default';
defaultoptions.slider          = 0;
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
defaultoptions.maskback        = [] ;
defaultoptions.maskfront        = [] ;
defaultoptions.mc          = false ;
defaultoptions.brightness  = false ;
defaultoptions.factorback  = ones(NoI,1);
defaultoptions.factorfront  = ones(NoI,1);
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
if length(options.ROIback)~=4
    error('unknown ROI option [%s]',options.ROIback)
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