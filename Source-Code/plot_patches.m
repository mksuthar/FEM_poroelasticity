function [k] = plot_patches(N,C,defor, xq,yq,xdis,ydis,dis_mag, pres_all,Sigma_x,Sigma_y, Sigma_xy)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dt = datestr(now,'mm_dd_yyyy HH-MM-SS AM');
dir_name = ['results_',dt];
mkdir(dir_name)

order = [1 5 2 6 3 7 4 8];

ix = 0;
iy = 0;
if sum(N(:,2)) > 0 && sum(N(:,3)) > 0
    ix = 2;
    iy = 3;
end

if sum(N(:,3)) > 0 && sum(N(:,4)) > 0
    ix = 3;
    iy = 4;
end

if sum(N(:,2)) > 0 && sum(N(:,4)) > 0
    ix = 2;
    iy = 4;
end

figure
plot(N(C(1,order),ix),N(C(1,order),iy),'-ok')
hold on
plot(N(C(1,9),ix),N(C(1,9),iy),'ok')
plot(N(C(1,[1 2 3 4]),ix),N(C(1,[1 2 3 4]),iy),'^k')
for i=2:length(C(:,1))
    plot(N(C(i,order),ix),N(C(i,order),iy),'-ok')
    plot(N(C(i,9),ix),N(C(i,9),iy),'ok')
    plot(N(C(i,[1 2 3 4]),ix),N(C(i,[1 2 3 4]),iy),'^k')
end
hold off
xlabel('\bfX')
ylabel('\bfY')
title('Geometry of the problem (circle = displacement node, triangle = pressure node)')
grid on
k = save_pdf(gcf,'geometry.pdf',dir_name)

N = [N(:,ix), N(:,iy)];

figure
subplot(2,1,1);
plot(N(:,1),N(:,2),'ok')
xlabel('\bfX')
ylabel('\bfY')
title('Initial condition')
grid on
subplot(2,1,2);
plot(defor(:,1),defor(:,2),'ok')
xlabel('\bfX')
ylabel('\bfY')
title('Simulated Displacement (x1 magnification)')
grid on
k = save_pdf(gcf,'dispRAW.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,xdis,11)
grid on
title('Displacement in X')
colorbar
set(ch,'edgecolor','none');
shading flat
xlabel('\bfX')
ylabel('\bfY')
k = save_pdf(gcf,'dispX.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,ydis,11)
title('Displacement in Y')
grid on
set(ch,'edgecolor','none');
shading flat
xlabel('\bfX')
ylabel('\bfY')
colorbar
k = save_pdf(gcf,'dispY.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,dis_mag,11)
title('Total Displacement magnitude')
grid on
xlabel('\bfX')
ylabel('\bfY')
set(ch,'edgecolor','none');
shading flat
colorbar
k = save_pdf(gcf,'dispM.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,pres_all,11)
xlabel('\bfX')
ylabel('\bfY')
title('Pressure magnitude')
grid on
set(ch,'edgecolor','none');
shading flat
colorbar
k = save_pdf(gcf,'pressM.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,Sigma_x,11)
title('\sigma_{xx}')
xlabel('\bfX')
ylabel('\bfY')
set(ch,'edgecolor','none');
shading flat
grid on
colorbar
k = save_pdf(gcf,'sigma_xx.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,Sigma_y,11)
xlabel('\bfX')
ylabel('\bfY')
title('\sigma_{yy}')
set(ch,'edgecolor','none');
shading flat
grid on
colorbar
k = save_pdf(gcf,'sigma_yy.pdf',dir_name)

figure
[ch,ch]=contourf(xq,yq,Sigma_xy,11)
xlabel('\bfX')
ylabel('\bfY')
title('\sigma_{xy}')
set(ch,'edgecolor','none');
shading flat
grid on
colorbar
k = save_pdf(gcf,'sigma_xy.pdf',dir_name)

%close all

end

function [k] = save_pdf(gcf,name,dir_name)
    h = gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    set(gcf,'Renderer','Zbuffer')
    print(gcf, '-dpdf', [dir_name,'/',name]);
    k = 0;
end

