classdef Plotter
    properties
        plates
        mode_list
    end
    properties (Constant, Access=private)
        CB_LABEL_FS = 14;     % colorbar title size
        CB_TICK_FS  = 11;     % colorbar tick size
        AX_LABEL_FS = 12;     % axes label size
    end

    methods
        function obj = Plotter(plates, mode_list)
            obj.plates    = plates;
            obj.mode_list = mode_list;
        end

        function plot_displacement(obj)
            plates = obj.plates;  mode_list = obj.mode_list;
            Nmode  = mode_list(end);
            nP     = numel(plates);

            % --- Figure 1: v(x,0) vs N ---
            figure('Name','v(x,0) vs N'); hold on;
            leg = strings(1,numel(mode_list));
            [~,v,xv,yv] = plates(1).solve_plate(mode_list(1));
            [~,iy0] = min(abs(yv));  % y=0 index
            for j = 1:numel(mode_list)
                [~,v,xv,~] = plates(1).solve_plate(mode_list(j));
                plot(xv, v(iy0,:), 'LineWidth', 1.4);
                leg(j) = sprintf('N = %d', mode_list(j));
            end
            grid off; xlabel('x'); ylabel('v(x,0)');
            % title('v(x,0) Midline Deflection for Multiple Truncations');
            legend(leg, 'Location','best'); hold off;
            set(gca,'YDir','reverse');

            % --- Figure 2: deformed field colored by v/u for each plate ---
            figure('Name','Deformed field (v and u)');
            t = tiledlayout(nP,2,'Padding','compact','TileSpacing','compact');

            for i = 1:nP
                [u,v,xv,yv] = plates(i).solve_plate(Nmode);
                [X,Y] = meshgrid(xv,yv);
                mag = sqrt(u.^2 + v.^2);
                sf  = 0.15 * max(plates(i).l, plates(i).h) / max(mag(:)+eps);
                Xd  = X + sf*u;   Yd = Y + sf*v;

                % v (left)
                ax = nexttile(t,(i-1)*2+1);
                surf(ax, Xd, Yd, 0*v, v, 'EdgeColor','none');
                set(ax,'YDir','reverse'); view(ax,2); axis(ax,'equal','tight');
                colormap(ax, parula); obj.setCbar(ax,'v');
                obj.symClim(ax, v); grid(ax,'off');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                % title(ax, sprintf('h/l=%.3g, N=%d (v)', plates(i).S^-1, Nmode));

                % u (right)
                ax = nexttile(t,(i-1)*2+2);
                surf(ax, Xd, Yd, 0*u, u, 'EdgeColor','none');
                set(ax,'YDir','reverse'); view(ax,2); axis(ax,'equal','tight');
                colormap(ax, parula); obj.setCbar(ax,'u');
                obj.symClim(ax, u); grid(ax,'off');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                % title(ax, sprintf('h/l=%.3g, N=%d (u)', plates(i).S^-1, Nmode));
            end
        end

        function plot_stresses(obj)
            plates = obj.plates;  nP = numel(plates);
            Nmode  = obj.mode_list(end);

            figure('Name','Stress distributions (deformed grid)');
            t = tiledlayout(nP,3,'Padding','compact','TileSpacing','compact');

            for i = 1:nP
                [u,v,xv,yv,sigx,sigy,tauxy] = plates(i).solve_plate(Nmode);
                [X,Y] = meshgrid(xv,yv);
                mag = sqrt(u.^2 + v.^2);
                sf  = 0.15 * max(plates(i).l, plates(i).h) / max(mag(:)+eps);
                Xd  = X + sf*u;   Yd = Y + sf*v;

                % σx
                ax = nexttile(t,(i-1)*3+1);
                surf(ax, Xd, Yd, 0*sigx, sigx, 'EdgeColor','none');
                set(ax,'YDir','reverse'); view(ax,2); axis(ax,'equal','tight');
                obj.symClim(ax, sigx); grid(ax,'off'); obj.setCbar(ax,'\sigma_x');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                % title(ax, sprintf('\\sigma_x (S=%.3g)', plates(i).S^-1));

                % σy
                ax = nexttile(t,(i-1)*3+2);
                surf(ax, Xd, Yd, 0*sigy, sigy, 'EdgeColor','none');
                set(ax,'YDir','reverse'); view(ax,2); axis(ax,'equal','tight');
                obj.symClim(ax, sigy); grid(ax,'off'); obj.setCbar(ax,'\sigma_y');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                % title(ax, sprintf('\\sigma_y'));

                % τxy
                ax = nexttile(t,(i-1)*3+3);
                surf(ax, Xd, Yd, 0*tauxy, tauxy, 'EdgeColor','none');
                set(ax,'YDir','reverse'); view(ax,2); axis(ax,'equal','tight');
                obj.symClim(ax, tauxy); grid(ax,'off'); obj.setCbar(ax,'\tau_{xy}');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                % title(ax, sprintf('\\tau_{xy}'));
            end
        end
    end

    methods (Access=private)
        function setCbar(obj, ax, label)
            cb = colorbar(ax);
            cb.Label.String     = label;
            cb.Label.FontWeight = 'bold';
            cb.Label.FontSize   = obj.CB_LABEL_FS;   % larger “ruler title”
            cb.FontWeight       = 'bold';            % bold tick labels
            cb.FontSize         = obj.CB_TICK_FS;    % tick size
        end

        function symClim(~, ax, Z)
            cmax = max(abs(Z(:)));
            if cmax>0, clim(ax,[-cmax cmax]); end
        end
    end
end
