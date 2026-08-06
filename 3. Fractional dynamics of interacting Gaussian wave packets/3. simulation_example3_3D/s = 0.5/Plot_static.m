% This program uses the data computed by equation3D_diffusion
% (MT_TSM_equation3D_static).m to plot the contour maps and x2
% cross-sections of the corresponding numerical solutions.

clear, close all;
load('simulation3D_diffusion_s=0.5.mat');
M = 13.5; % Plotting area parameter

%% Compute the values of numerical solution over [-M,M]^3
% (length(pt)=101, plot=0.4s; lenght(pt)=151, plot=0.9s; lenght(pt)=201, plot=2s)
% (length(pt)=301, plot=6.5s; length(pt)=401, plot=15.1s; lenght(pt)=501, plot=32.8s)
pt = -M:M/250:M; % the points [pt1,pt2,pt3] = ndgrid(pt,pt,pt) used to plot the isosurface figure
[pt_plot1, pt_plot2, pt_plot3] = meshgrid(pt,pt,pt); % the node values ​​used for plotting (must use "meshgrid")

[NN,ppt] = meshgrid(-N:N-1,pt); MTmatrix = MT(NN,mu*ppt); % common matrix, which is used to compute the value of sum(j=-N:N-1)sum(k=-N:N-1)sum(l=-N:N-1) U(j,k,l)*MT(j,mu*x)*MT(k,mu*y)*MT(l,mu*z) on points [pt1,pt2,pt3] based on U
unum_node_pt = zeros(length(pt),length(pt),length(pt),inter+1); % unum_nodept(:,:,:,q) denotes the values of \psi(x,y,z,(q-1)*(Tstep/inter)*dt) on points [pt1,pt2,pt3]
for m = 1:inter+1
    unum_node_pt_temp = tensorprod(unum_coeff(:,:,:,m),MTmatrix,3,2);
    unum_node_pt(:,:,:,m) = pagemtimes(pagemtimes(MTmatrix,unum_node_pt_temp),MTmatrix.');
end

%% Compute the values of numerical solution over [-M,M]×0×[-M,M]
pt_cross = -M:M/1000:M; [pt_cross2,pt_cross1] = meshgrid(pt_cross); % the points used to plot the cross section x2=0

B = (1i).^(-N:N-1)*sqrt(2/pi); % B = (MT(-N),...,MT(N-1))
[NN,ppt_cross] = meshgrid(-N:N-1,pt_cross); MTmatrix = MT(NN,mu*ppt_cross);
unum_node_pt_cross = zeros(length(pt_cross),length(pt_cross),inter+1); % unum_node_pt_cross(:,:,q) denotes the values of \psi(x,0,z,(q-1)*(Tstep/inter)*dt) on points [pt_cross1,pt_cross2]
for m = 1:inter+1
    % \psi(x,0,z) = sum(p=-N:N-1)sum(q=-N:N-1)sum(r=-N:N-1) U(p,q,r)*MT(p,mu*x)*MT(q,0)*MT(r,mu*z) 
    % = sum(p=-N:N-1)sum(r=-N:N-1) (sum(q=-N:N-1)U(p,q,r)*MT(q,0))*MT(p,mu*x)*MT(r,mu*z)
    % = sum(p=-N:N-1)sum(r=-N:N-1) A(p,r)*MT(p,mu*x)*MT(r,mu*z)
    A = squeeze(sum(unum_coeff(:,:,:,m).*B,2));
    unum_node_pt_cross(:,:,m) = MTmatrix*A*MTmatrix.';
end

%% Plot the isosurface figure and the cross section x2=0
% isosurface level
level = 1.5*1e-2; % plot |u|=level

% standardize the coordinate axis scale
xticks = -10:5:10;
yticks = -10:10:10;
zticks = -10:5:10;

% standardize the view
v = [160,10];

figure('Position', [0, 50, 1800, 800])
plot_par = tight_subplot(2, 4, [0.065 0.05], [0.015 0.035], [0.04 0.015]);
for m = 1:8
    axes(plot_par(m))
    if m<= 4
        %% plot isosurfaces
        data = abs(unum_node_pt(:,:,:,m)); % original data, its arrangement corresponds to ndgrid(pt,pt,pt).
        data_plot = permute(data, [2,1,3]); % transform the original data into a data arrangement adapted for "meshgrid"
        pc = patch(isosurface(pt_plot1, pt_plot2, pt_plot3, data_plot, level));
        isonormals(pt_plot1, pt_plot2, pt_plot3, data_plot, pc); % plot isosurface

        % set the color and transparency of the isosurface
        pc.FaceColor = [0.4, 0.7, 1];
        pc.FaceAlpha = 0.35; % set the transparency level: the higher the value, the less transparent (设置透明度, 越大越不透明)
        pc.EdgeColor = 'none'; % hide border (隐藏边线)

        % lighting effect (光照效果)
        camlight headlight,
        lighting phong

        % grid, box, view
        grid off, box on, view(v)

        % label and axis
        hx = xlabel('$x_1$','Interpreter','latex','FontSize',14);
        hy = ylabel('$x_2$','Interpreter','latex','FontSize',14);
        hz = zlabel('$x_3$','Interpreter','latex','FontSize',14);
        set(gca,'XTick',xticks); set(gca,'XTickLabel',num2str(xticks'));
        set(gca,'YTick',yticks); set(gca,'YTickLabel',num2str(yticks'));
        set(gca,'ZTick',zticks); set(gca,'ZTickLabel',num2str(zticks'));
        set(gca,'FontSize',14);
        axis([-M, M, -M, M, -M, M]); daspect([1, 1, 1])

        drawnow;
        hy.Units = 'data';
        pos = hy.Position;
        pos(1) = pos(1) - 5.0;
        pos(2) = pos(2) - 11.0;
        hy.Position = pos; % modify the position of ylabel

        % subtitle
        title(sprintf('$t = %d$', (m-1)*(Tstep/inter)*dt),'Interpreter','latex', 'FontSize', 20);
    else
        %% plot the cross section x2=0
        surf(pt_cross1, zeros(size(pt_cross1)), pt_cross2, abs(unum_node_pt_cross(:,:,m-4)));
        shading interp; % smooth interpolation
        colorbar

        % grid, box, view
        grid off, box on, view(v)

        % label and axis
        hx = xlabel('$x_1$','Interpreter','latex','FontSize',14);
        hy = ylabel('$x_2$','Interpreter','latex','FontSize',14);
        hz = zlabel('$x_3$','Interpreter','latex','FontSize',14);
        set(gca,'XTick',xticks); set(gca,'XTickLabel',num2str(xticks'));
        set(gca,'YTick',0); set(gca,'YTickLabel',0);
        set(gca,'ZTick',zticks); set(gca,'ZTickLabel',num2str(zticks'));
        set(gca,'FontSize',14);
        axis([-M, M, -3*M/5, 3*M/5, -M, M]); daspect([1, 1, 1])

        drawnow;
        hy.Units = 'data';
        pos = hy.Position;
        pos(1) = pos(1) - 3.2;
        pos(2) = pos(2) - 7.0;
        hy.Position = pos; % modify the position of ylabel
    end
end