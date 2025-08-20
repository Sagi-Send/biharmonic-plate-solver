classdef Plotter
    properties
        plate
        mode_list
    end

    methods
        function obj = Plotter(plate, mode_list)
            obj.plate     = plate;
            obj.mode_list = mode_list;
        end

        function plot_displacement(obj)
            plate = obj.plate;  mode_list = obj.mode_list;
            Nmode_wholeplot = mode_list(end);

            % --- Figure 1: v(x,0) vs truncation ---
            figure('Name','v(x,0) vs N'); hold on;
            leg = strings(1,numel(mode_list));
            for i = 1:numel(mode_list)
                Nmode = mode_list(i);
                [~,v,xv_plot,yv_plot] = plate.solve_plate(Nmode);
                [~,iy0] = min(abs(yv_plot-0));
                plot(xv_plot, v(iy0,:), 'LineWidth', 1.4);
                leg(i) = sprintf('N = %d', Nmode);
            end
            grid off; xlabel('x'); ylabel('v(x,0)');
            title(sprintf('v(x,0) for multiple truncations (Nx=%d, Ny=%d)',...
                plate.Nx, plate.Ny));
            legend(leg, 'Location','best'); hold off;
            set(gca,'YDir','reverse')

            % --- Figure 2: deformed field colored by v(x,y) ---
            [u,v,xv,yv] = plate.solve_plate(Nmode_wholeplot);
            [X,Y] = meshgrid(xv,yv);

            mag = sqrt(u.^2 + v.^2);
            sf  = 0.15 * max(plate.l,plate.h) / max(mag(:) + eps);
            Xd = X + sf*u;   Yd = Y + sf*v;

            figure('Name','Deformed field (v)');
            surf(Xd, Yd, zeros(size(v)), v, 'EdgeColor','none');
            view(2); axis equal tight; colormap(parula); colorbar;
            cmaxV = max(abs(v(:))); if cmaxV>0, caxis([-cmaxV cmaxV]); end
            xlabel('x'); ylabel('y');
            title(sprintf('Deformed mesh + v-coloring (scale=%.3g), N=%d',...
                sf, Nmode_wholeplot));
            set(gca,'YDir','reverse')
        end

        function plot_stresses(obj)
            plate = obj.plate;
            Nmode_wholeplot = obj.mode_list(end);

            % solve and build deformed coords
            [u,v,xv,yv,sigx,sigy,tauxy] = plate.solve_plate...
                (Nmode_wholeplot);
            [X,Y] = meshgrid(xv,yv);

            mag = sqrt(u.^2 + v.^2);
            sf  = 0.15 * max(plate.l,plate.h) / max(mag(:) + eps);
            Xd = X + sf*u;   Yd = Y + sf*v;

            figure('Name','Stress distributions (deformed grid)');

            % σx
            subplot(1,3,1);
            surf(Xd, Yd, 0*sigx, sigx, 'EdgeColor','none');
            view(2); axis equal tight; colorbar;
            cmax = max(abs(sigx(:))); if cmax>0, caxis([-cmax cmax]); end
            xlabel('x'); ylabel('y'); title('\sigma_x');
            set(gca,'YDir','reverse')

            % σy
            subplot(1,3,2);
            surf(Xd, Yd, 0*sigy, sigy, 'EdgeColor','none');
            view(2); axis equal tight; colorbar;
            cmax = max(abs(sigy(:))); if cmax>0, caxis([-cmax cmax]); end
            xlabel('x'); ylabel('y'); title('\sigma_y');
            set(gca,'YDir','reverse')

            % τxy
            subplot(1,3,3);
            surf(Xd, Yd, 0*tauxy, tauxy, 'EdgeColor','none');
            view(2); axis equal tight; colorbar;
            cmax = max(abs(tauxy(:))); if cmax>0, caxis([-cmax cmax]); end
            xlabel('x'); ylabel('y'); title('\tau_{xy}');
            set(gca,'YDir','reverse')
        end
    end
end
