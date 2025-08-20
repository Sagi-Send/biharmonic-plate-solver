%% Distributed force on a clamped plate: numerical implementation (multi-mode + whole-plate plots)
clear; clc;

%% Parameters (numeric)
k  = 1.0;          % exponent in w(x) = w0*exp(k*x/l)
l  = 10.0;          % span
h  = 1.0;          % half-thickness (y in [-h,h])
E  = 1.0;          % Young's modulus
nu = 0.20;         % Poisson's ratio
w0 = 1.0;          % load amplitude

% grid
Nx = 600;          % x points
Ny = 241;          % y points

% which truncations to plot (line at y=0) and which to show on whole plate
mode_list       = [0 1 2 4 8 16 32 50];
Nmode_wholeplot = mode_list(end);

%% Figure 1: v(x,0) for multiple truncations
xv = linspace(0,l,Nx);
figure('Name','v(x,0) vs N'); hold on;
leg = strings(1,numel(mode_list));

for i = 1:numel(mode_list)
    Nmode = mode_list(i);
    [~,v,xv_plot,yv_plot] = solve_plate(Nmode, k,l,h,E,nu,w0, Nx,Ny);
    [~,iy0] = min(abs(yv_plot-0));
    plot(xv_plot, v(iy0,:), 'LineWidth', 1.4);
    leg(i) = sprintf('N = %d', Nmode);
end
grid off; xlabel('x'); ylabel('v(x,0)');
title(sprintf('v(x,0) for multiple truncations (Nx=%d, Ny=%d)', Nx, Ny));
legend(leg, 'Location','best'); hold off;

%% Figure 2: whole-plate plots for a selected N
[u,v,xv,yv] = solve_plate(Nmode_wholeplot, k,l,h,E,nu,w0, Nx,Ny);
[X,Y] = meshgrid(xv,yv);

figure('Name','Whole-plate fields');
% 1) filled contour of v
subplot(1,2,1);
contourf(xv, yv, v, 30, 'LineColor','none');
axis equal tight;
colorbar; xlabel('x'); ylabel('y'); title('v(x,y) contour');

% 3) in-plane deformed mesh (u,v)
subplot(1,2,2); hold on;
mag = sqrt(u.^2 + v.^2);
sf  = 0.15 * max(l,h) / max(mag(:) + eps);   % deformation scale
Xd = X + sf*u;   Yd = Y + sf*v;
plot(X,  Y,  'k:', X',  Y',  'k:');          % original grid
plot(Xd, Yd, 'b-', Xd', Yd', 'b-');          % deformed grid
axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('Deformed mesh (scale=%.2g), N=%d', sf, Nmode_wholeplot));
hold off;

%% ===== solver: returns u(x,y), v(x,y), and grids =====
function [u,v,xv,yv] = solve_plate(Nmode, k,l,h,E,nu,w0, Nx,Ny)

    % constants and grids
    G = E/(2*(1+nu));
    I = (2/3)*h^3;
    xv = linspace(0,l,Nx);
    yv = linspace(-h,h,Ny);
    [X,Y] = meshgrid(xv,yv);

    dx = xv(2)-xv(1);
    dy = yv(2)-yv(1);

    [~,iy0] = min(abs(yv-0));   % y=0 index
    ixL = Nx;                   % x=l index

    % Fourier coefficients for w(x) = w0*exp(k*x/l)
    a0 = w0*(exp(k)-1)/k;
    n  = 1:Nmode;
    an = 2*w0*k*(exp(k)-1) ./ (k^2 + (2*pi*n).^2);
    bn = -(exp(k)-1)*w0*(4*pi*n) ./ (k^2 + (2*pi*n).^2);

    % n=0 stresses
    sigx  = (a0/(2*I))*( (l - X).*X.*Y ) + (a0/(2*I))*((2/3)*Y.^3 - (2/5)*h^2.*Y);
    sigy  = -(a0/(2*I))*((1/3)*Y.^3 - h^2.*Y + (2/3)*h^3);
    tauxy = -(a0/(2*I))*((h^2 - Y.^2).*(X - l/2));

    % add modal stresses
    for m = 1:Nmode
        alpha = 2*pi*m/l;
        ah = alpha*h; ch = cosh(ah); sh = sinh(ah);
        S2 = sinh(2*ah);
        Dp = S2 + 2*alpha*h;
        Dm = S2 - 2*alpha*h;

        cy = cosh(alpha*Y);
        sy = sinh(alpha*Y);
        cX = cos(alpha*X);
        sX = sin(alpha*X);

        A11 = ( (ah*ch - sh).*cy - alpha*Y.*sy.*sh )/Dp ...
            - ( (ah*sh - ch).*sy - alpha*Y.*cy.*ch )/Dm;

        A12 = ( (ah*sh + ch).*sy - alpha*Y.*cy.*ch )/Dm ...
            - ( (ah*ch + sh).*cy - alpha*Y.*sy.*sh )/Dp;

        % a_n (cos family)
        sigx  = sigx  + an(m).*cX .* A11;
        sigy  = sigy  + an(m).*cX .* A12;
        tauxy = tauxy + an(m).*sX .* ( (ah*ch).*sy - alpha*Y.*cy.*sh )/Dp ...
                      - an(m).*sX .* ( (ah*sh).*cy - alpha*Y.*sy.*ch )/Dm;

        % b_n (sin family)
        sigx  = sigx  + bn(m).*sX .* A11;
        sigy  = sigy  + bn(m).*sX .* A12;
        tauxy = tauxy + bn(m).*cX .* ( (ah*sh).*cy - alpha*Y.*sy.*ch )/Dm ...
                      - bn(m).*cX .* ( (ah*ch).*sy - alpha*Y.*cy.*sh )/Dp;
    end

    % particular displacements by line integration
    epsx = (sigx - nu*sigy)/E;
    epsy = (sigy - nu*sigx)/E;

    u0 = zeros(Ny,Nx);
    for j = 1:Ny
        u0(j,:) = cumtrapz(xv, epsx(j,:));
    end

    v0 = zeros(Ny,Nx);
    for i = 1:Nx
        col = epsy(:,i).';
        vtmp = cumtrapz(yv, col);
        v0(:,i) = vtmp - vtmp(iy0);  % set v0(x,0)=0
    end

    % Section-averaged split: fit R ≈ F(x) + G(y) - K in LS sense
    [Uy,~] = gradient(u0, dy, dx);
    [~,Vx] = gradient(v0, dy, dx);
    R = Uy + Vx - tauxy./G;

    K     = mean(R(:));
    Fbase = mean(R,1) - K;      % 1xNx
    Gbase = mean(R,2).' - K;    % 1xNy

    % f2'(x) = -F(x), f1'(y) = -G(y)
    f2 = -cumtrapz(xv, Fbase);      f2 = f2(:).';     % 1xNx
    f1 = -cumtrapz(yv, Gbase.');    f1 = f1(:);       % Nyx1
    f1 = f1 - f1(iy0);                              % set f1(0)=0

    % add integration functions
    u = u0 + repmat(f1, 1, Nx);
    v = v0 + repmat(f2, Ny, 1);

    % clamp at (x=l,y=0): zero translations
    u = u - u(iy0,ixL);
    v = v - v(iy0,ixL);

    % enforce zero rotation at clamp with backward difference
    S = (v(iy0,ixL) - v(iy0,ixL-1))/dx;    % current slope at x=l
    u = u + S * (yv.' * ones(1,Nx));       % add S*y to u
    v = v - S * (ones(Ny,1) * (xv - l));   % subtract S*(x-l) from v
end
