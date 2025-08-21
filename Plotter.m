classdef Plotter
    properties
        plates
        mode_list
    end
    properties (Constant, Access=private)
        CB_LABEL_FS = 16;     % colorbar title size
        CB_TICK_FS  = 13;     % colorbar tick size
        AX_LABEL_FS = 14;     % axes label size
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

        function plot_displacement_exp(obj)
            plates = obj.plates;
            nP     = numel(plates);
        
            figure('Name','Deformed field (v and u)');
            t = tiledlayout(nP,2,'Padding','compact','TileSpacing','compact');
        
            for i = 1:nP
                [u,v,xv,yv] = plates(i).solve_plate_exp();
                [X,Y] = meshgrid(xv,yv);
                mag = sqrt(u.^2 + v.^2);
                sf  = 0.15 * max(plates(i).l, plates(i).h) / max(mag(:)+eps);
                Xd  = X + sf*u;   Yd = Y + sf*v;
        
                % v (left)
                ax = nexttile(t,(i-1)*2+1);
                surf(ax, Xd, Yd, 0*v, v, 'EdgeColor','none');
                view(ax,2); axis(ax,'equal','tight');
                colormap(ax, parula); obj.setCbar(ax,'v');
                obj.symClim(ax, v); grid(ax,'off');
                set(ax,'FontSize',obj.CB_TICK_FS);
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                set(gca,'YDir','reverse');
        
                % u (right)
                ax = nexttile(t,(i-1)*2+2);
                surf(ax, Xd, Yd, 0*u, u, 'EdgeColor','none');
                view(ax,2); axis(ax,'equal','tight');
                colormap(ax, parula); obj.setCbar(ax,'u');
                obj.symClim(ax, u); grid(ax,'off');
                set(ax,'FontSize',obj.CB_TICK_FS);
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                set(gca,'YDir','reverse');
            end
        end

        function plot_stresses_exp(obj)
            plates = obj.plates;  nP = numel(plates);
        
            figure('Name','Stress distributions (deformed grid)');
            t = tiledlayout(nP,3,'Padding','compact','TileSpacing','compact');
        
            for i = 1:nP
                [u,v,xv,yv,sigx,sigy,tauxy] = plates(i).solve_plate_exp();
                [X,Y] = meshgrid(xv,yv);
                mag = sqrt(u.^2 + v.^2);
                sf  = 0.15 * max(plates(i).l, plates(i).h) / max(mag(:)+eps);
                Xd  = X + sf*u;   Yd = Y + sf*v;
        
                % σx
                ax = nexttile(t,(i-1)*3+1);
                surf(ax, Xd, Yd, 0*sigx, sigx, 'EdgeColor','none');
                view(ax,2); axis(ax,'equal','tight');
                obj.symClim(ax, sigx); grid(ax,'off'); obj.setCbar(ax,'\sigma_x');
                set(ax,'FontSize',obj.CB_TICK_FS);
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                set(gca,'YDir','reverse');
        
                % σy
                ax = nexttile(t,(i-1)*3+2);
                surf(ax, Xd, Yd, 0*sigy, sigy, 'EdgeColor','none');
                view(ax,2); axis(ax,'equal','tight');
                obj.symClim(ax, sigy); grid(ax,'off'); obj.setCbar(ax,'\sigma_y');
                set(ax,'FontSize',obj.CB_TICK_FS);
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                set(gca,'YDir','reverse');
        
                % τxy
                ax = nexttile(t,(i-1)*3+3);
                surf(ax, Xd, Yd, 0*tauxy, tauxy, 'EdgeColor','none');
                view(ax,2); axis(ax,'equal','tight');
                obj.symClim(ax, tauxy); grid(ax,'off'); obj.setCbar(ax,'\tau_{xy}');
                set(ax,'FontSize',obj.CB_TICK_FS);
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                set(gca,'YDir','reverse');
            end
        end

        function plot_disp_concentrated(obj)
            % Plot u and v (uy0 clamp) for concentrated load, side-by-side per plate
            plates = obj.plates;
            nP     = numel(plates);
        
            % ==== Figure A: v only, one row x two columns (Plate 1 vs Plate 2) ====
            figure('Name','v (uy0) on displaced shape — concentrated load: Plate 1 vs Plate 2');
            t = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
            
            nShow = min(nP,2);                          % show first two plates
            U = cell(1,nShow); V = cell(1,nShow);
            XV = cell(1,nShow); YV = cell(1,nShow);
            MAG = cell(1,nShow);
            
            % Solve once for the plates to get a shared color scale
            for k = 1:nShow
                [U{k}, V{k}, XV{k}, YV{k}] = plates(k).solve_plate_concentrated("uy0");
                MAG{k} = hypot(U{k}, V{k});
            end
            vabs = max( abs([V{1}(:); V{min(2,nShow)}(:)]) );
            
            for k = 1:nShow
                [Xk,Yk] = meshgrid(XV{k}, YV{k});
                sf  = 0.15 * max(plates(k).l, plates(k).h) / max([MAG{k}(:); eps]);
            
                ax = nexttile(t,k);
                Xd = Xk + sf*U{k};  Yd = Yk + sf*V{k};
                surf(ax, Xd, Yd, zeros(size(Xd)), V{k}, 'EdgeColor','none'); hold(ax,'on');
                set(ax,'YDir','reverse'); view(ax,2); axis(ax,'equal','tight');
                colormap(ax, parula); grid(ax,'off'); clim(ax,[-vabs vabs]);
                obj.setCbar(ax,'v');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                title(ax, sprintf('Plate %d', k));
            end
            
            
            % ==== Single-row figure: v(x,0) overlays for Plates 1 & 2 ====
            figure('Name','v(x,0): vx0 vs uy0 — Plates 1 & 2');
            t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
            
            % -------- Tile 1: Plate 1 (v at y=0 vs x) --------
            [u1,v1,xv1,yv1] = plates(1).solve_plate_concentrated("vx0");
            [u2,v2,~,   ~ ] = plates(1).solve_plate_concentrated("uy0");
            [~,iy0] = min(abs(yv1));  % y = 0 index
            
            ax = nexttile(t,1);
            plot(ax, xv1, v1(iy0,:), 'LineWidth',1.5, 'DisplayName','v,x=0'); hold(ax,'on');
            plot(ax, xv1, v2(iy0,:), 'LineWidth',1.5, 'DisplayName','u,y=0');
            grid(ax,'off'); xlim(ax,[xv1(1) xv1(end)]);
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'v(x,0)','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            legend(ax,'Location','best');
            
            % -------- Tile 2: Plate 1 (u along y at x = l) --------
            x0 = plates(1).l;                     % choose the section (e.g., right end)
            [~,ix0] = min(abs(xv1 - x0));         % x index for x0
            
            ax = nexttile(t,2);
            plot(ax, u1(:,ix0), yv1,  'LineWidth',1.5, 'DisplayName','v,x=0'); hold(ax,'on');
            plot(ax, u2(:,ix0), yv1,  'LineWidth',1.5, 'DisplayName','u,y=0');
            grid(ax,'off'); xlim(ax,[u1(1,ix0) u1(end,ix0)]);
            xlabel(ax, sprintf('u(y, x=%.3g)', x0),'FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax, 'y', 'FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            legend(ax,'Location','best');
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
