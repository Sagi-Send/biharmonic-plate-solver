classdef Plate
    properties (Constant)
        h  = 0.5;   % half-thickness (total = 2h)
        E  = 20.0;
        nu = 0.20;
        w0 = 1.0;
    end
    properties (Dependent)
        l                  % span; computed from S and h
    end
    properties
        k  = 1.0;              % exponential load exponent
        S  = 1.0;              % slenderness S = l/(2h)
        Nx = 600;   Ny = 241;  % grid

    end

    methods
        function obj = Plate(k, slenderness, Nx, Ny)
            obj.k  = k;     obj.S  = slenderness;
            obj.Nx = Nx;    obj.Ny = Ny;
        end

        function L = get.l(obj)
            L = 2*obj.h * obj.S;
        end

        % ---- solver: returns u,v, grids, and stresses ----
        function [u,v,xv,yv,sigx,sigy,tauxy] = solve_plate(obj, Nmode)
            k = obj.k;  l = obj.l;  E = obj.E;  nu = obj.nu;  w0 = obj.w0;
            Nx = obj.Nx; Ny = obj.Ny; h = obj.h;

            G = E/(2*(1+nu));
            I = (2/3)*h^3;
            xv = linspace(0,l,Nx);
            yv = linspace(-h,h,Ny);
            [X,Y] = meshgrid(xv,yv);

            dx = xv(2)-xv(1);
            dy = yv(2)-yv(1);

            [~,iy0] = min(abs(yv-0));
            ixL = Nx;

            % Fourier of w(x) = w0 * exp(k*x/l)
            a0 = -w0*(exp(k)-1)/k;
            n  = 1:Nmode;
            an = -2*w0*k*(exp(k)-1) ./ (k^2 + (2*pi*n).^2);
            bn = (exp(k)-1)*w0*(4*pi*n) ./ (k^2 + (2*pi*n).^2);

            % n=0 stresses
            sigx  = (a0/(2*I))*((l - X).*X.*Y) + (a0/(2*I))*((2/3)*Y.^3 - (2/5)*h^2.*Y);
            sigy  = -(a0/(2*I))*((1/3)*Y.^3 - h^2.*Y + (2/3)*h^3);
            tauxy = -(a0/(2*I))*((h^2 - Y.^2).*(X - l/2));

            % modal stresses
            for m = 1:Nmode
                alpha = 2*pi*m/l;
                ah = alpha*h; ch = cosh(ah); sh = sinh(ah);
                S2 = sinh(2*ah);
                Dp = S2 + 2*alpha*h;   % D_+
                Dm = S2 - 2*alpha*h;   % D_-

                cy = cosh(alpha*Y);  sy = sinh(alpha*Y);
                cX = cos(alpha*X);   sX = sin(alpha*X);

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

            % strains (plane stress) -> integrate to u0,v0
            epsx = (sigx - nu*sigy)/E;
            epsy = (sigy - nu*sigx)/E;

            u0 = zeros(Ny,Nx);
            for j = 1:Ny
                u0(j,:) = cumtrapz(xv, epsx(j,:));
            end

            v0 = zeros(Ny,Nx);
            for i = 1:Nx
                col  = epsy(:,i);
                vtmp = cumtrapz(yv, col);
                v0(:,i) = vtmp - vtmp(iy0);  % enforce v0(x,0)=0
            end

            % shear compatibility correction via separated fit
            [~,Uy] = gradient(u0, dx, dy);   % Uy = du/dy
            [Vx,~] = gradient(v0, dx, dy);   % Vx = dv/dx
            R = Uy + Vx - tauxy./G;

            K     = mean(R(:));
            Fbase = mean(R,1) - K;          % function of x
            Gbase = mean(R,2).' - K;        % function of y

            f2 = -cumtrapz(xv, Fbase);  f2 = f2(:).';
            f1 = -cumtrapz(yv, Gbase.'); f1 = f1(:);
            f1 = f1 - f1(iy0);           % anchor f1(0)=0

            u = u0 + repmat(f1, 1, Nx);
            v = v0 + repmat(f2, Ny, 1);

            % clamp (x=l,y=0): zero translations + zero slope v_x
            u = u - u(iy0,ixL);
            v = v - v(iy0,ixL);
            Srot = (v(iy0,ixL) - v(iy0,ixL-1))/dx;
            u = u + Srot * (yv.' * ones(1,Nx));
            v = v - Srot * (ones(Ny,1) * (xv - l));
        end
    end
end