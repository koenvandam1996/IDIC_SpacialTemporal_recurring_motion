load('C:\Users\s149743\Documents\globalDIC - v1.0.rc - Release Candidate\15-15backgr.mat')
[n, m] = size(Ux);
dx = mean(diff(x));
dy = mean(diff(y));
Ux(gdic.mask)=NaN;
Uy(gdic.mask)=NaN;
edgetop=0;
edgebot=0;
Ux((gdic.mask(1:end-edgetop)+edgetop))=NaN;
Uy((gdic.mask(1:end-edgetop)+edgetop))=NaN;
Ux((gdic.mask(edgebot+1:end)-edgebot))=NaN;
Uy((gdic.mask(edgebot+1:end)-edgebot))=NaN;
[F.xx, F.xy] = gradient(Ux,dx,dy);
[F.yx, F.yy] = gradient(Uy,dx,dy);
ZERO = zeros(n,m,'int8');
F.xx = F.xx + 1;
F.yy = F.yy + 1;
F.xx(gdic.mask)=NaN;
F.xy(gdic.mask)=NaN;
F.yx(gdic.mask)=NaN;
F.yy(gdic.mask)=NaN;
F.xxtot=zeros(n,m);
F.xytot=zeros(n,m);
F.yxtot=zeros(n,m);
F.yytot=zeros(n,m);
for i=1:n
    for j=1:400
        ifirst=i;
        jfirst=j;
        displx=Ux(round(ifirst),round(jfirst));
        disply=Uy(round(ifirst),round(jfirst));
        while ~isnan(displx) && ~isnan(disply) && jfirst<m-displx %&& ifirst<(n+disply)
%             if F.xxtot(round(ifirst),round(jfirst))==0
%                 F.xxtot(round(ifirst),round(jfirst))=((F.xx(round(ifirst),round(jfirst))));
%                 F.xytot(round(ifirst),round(jfirst))=(F.xy(round(ifirst),round(jfirst)));
%                 F.yxtot(round(ifirst),round(jfirst))=(F.yx(round(ifirst),round(jfirst)));
%                 F.yytot(round(ifirst),round(jfirst))=((F.yy(round(ifirst),round(jfirst))));
%             end
            if F.xxtot(round(ifirst),round(jfirst))==0 && jfirst<=200
                F.xxtot(round(ifirst),round(jfirst))=((F.xx(round(ifirst),round(jfirst))-1)*((jfirst)/22))+1;
                F.xytot(round(ifirst),round(jfirst))=(F.xy(round(ifirst),round(jfirst))*(jfirst))/22;
                F.yxtot(round(ifirst),round(jfirst))=(F.yx(round(ifirst),round(jfirst))*(jfirst))/22;
                F.yytot(round(ifirst),round(jfirst))=((F.yy(round(ifirst),round(jfirst))-1)*((jfirst)/22))+1;
            elseif F.xxtot(round(ifirst),round(jfirst))==0 && jfirst>200
                F.xxtot(round(ifirst),round(jfirst))=F.xx(round(ifirst),round(jfirst-1));
                F.xytot(round(ifirst),round(jfirst))=F.xy(round(ifirst),round(jfirst-1));
                F.yxtot(round(ifirst),round(jfirst))=F.yx(round(ifirst),round(jfirst-1));
                F.yytot(round(ifirst),round(jfirst))=F.yy(round(ifirst),round(jfirst-1));
            end
            
            isec=ifirst+disply;
            jsec=jfirst+displx;
            if F.xxtot(round(isec),round(jsec))==0
            F.xxtot(round(isec),round(jsec))=F.xxtot(round(ifirst),round(jfirst))*F.xx(round(isec),round(jsec))+F.xytot(round(ifirst),round(jfirst))*F.yx(round(isec),round(jsec));
            F.xytot(round(isec),round(jsec))=F.xxtot(round(ifirst),round(jfirst))*F.xy(round(isec),round(jsec))+F.xytot(round(ifirst),round(jfirst))*F.yy(round(isec),round(jsec));
            F.yxtot(round(isec),round(jsec))=F.yxtot(round(ifirst),round(jfirst))*F.xx(round(isec),round(jsec))+F.yytot(round(ifirst),round(jfirst))*F.yx(round(isec),round(jsec));
            F.yytot(round(isec),round(jsec))=F.yxtot(round(ifirst),round(jfirst))*F.xy(round(isec),round(jsec))+F.yytot(round(ifirst),round(jfirst))*F.yy(round(isec),round(jsec));
            else 
                F.xxtot(round(isec),round(jsec))=(F.xxtot(round(isec),round(jsec))+(F.xxtot(round(ifirst),round(jfirst))*F.xx(round(isec),round(jsec))+F.xytot(round(ifirst),round(jfirst))*F.yx(round(isec),round(jsec))))/2;
                F.xytot(round(isec),round(jsec))=(F.xytot(round(isec),round(jsec))+(F.xxtot(round(ifirst),round(jfirst))*F.xy(round(isec),round(jsec))+F.xytot(round(ifirst),round(jfirst))*F.yy(round(isec),round(jsec))))/2;
                F.yxtot(round(isec),round(jsec))=(F.yxtot(round(isec),round(jsec))+(F.yxtot(round(ifirst),round(jfirst))*F.xx(round(isec),round(jsec))+F.yytot(round(ifirst),round(jfirst))*F.yx(round(isec),round(jsec))))/2;
                F.yytot(round(isec),round(jsec))=(F.yytot(round(isec),round(jsec))+(F.yxtot(round(ifirst),round(jfirst))*F.xy(round(isec),round(jsec))+F.yytot(round(ifirst),round(jfirst))*F.yy(round(isec),round(jsec))))/2;
            end
            ifirst=isec;
            jfirst=jsec;
            displx=Ux(round(ifirst),round(jfirst));
            disply=Uy(round(ifirst),round(jfirst));
              
        end
    end
end
F.xxtot(F.xxtot == 0) = NaN;
F.xytot(F.xytot == 0) = NaN;
F.yxtot(F.yxtot == 0) = NaN;
F.yytot(F.yytot == 0) = NaN;

% Cauchy Green deformation tensor
C.xx = F.xxtot .* F.xxtot + F.yxtot .* F.yxtot;
C.xy = F.xxtot .* F.xytot + F.yxtot .* F.yytot;
C.yx = F.yytot .* F.yxtot + F.xytot .* F.xxtot;
C.yy = F.yytot .* F.yytot + F.xytot .* F.xytot;
E.xxtot = 0.5 * (C.xx - 1);
E.xytot = 0.5 * (C.xy);
E.yxtot = 0.5 * (C.yx);
E.yytot = 0.5 * (C.yy - 1);

figure
ax1=subplot(3,1,1);
imagesc(x1,y1,crs.Manf(:,:,2));
colormap(ax1,gray(256));
title('Strain Fields \epsilon_x','Fontsize',fontsize);
%xlabel('x [\mum]');
%ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,1e2*E.xxtot,'alphadata',~isnan(E.xxtot));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
caxis(ax2,[-20 20]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));

ax1=subplot(3,1,2);
imagesc(x1,y1,crs.Manf(:,:,2));
colormap(ax1,gray(256));
title('Strain Fields \epsilon_y','Fontsize',fontsize);
%xlabel('x [\mum]');
%ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,1e2*E.yytot,'alphadata',~isnan(E.yytot));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));
ax1=subplot(3,1,3);
imagesc(x1,y1,crs.Manf(:,:,2));
colormap(ax1,gray(256));
title('Strain Fields \epsilon_{xy} and \epsilon_{yx}','Fontsize',fontsize);
%xlabel('x [\mum]');
%ylabel('y [\mum]');
ax2 = axes;
imagesc(x,y,1e2*E.xytot,'alphadata',~isnan(E.xytot));
colormap(ax2,plotop.colormap);
h=colorbar;
h.Label.String='\epsilon [%]';
% caxis(ax2,[-10 10]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
linkaxes([ax1,ax2],'xy');
hFig.ResizeFcn = @(varargin) ...
    set(ax1,'Position',get(ax2, 'position'));
