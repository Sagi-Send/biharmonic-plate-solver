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
        
                % Decide if this plate is "slender"
                lambda = plates(i).l / plates(i).h;      % slenderness (l/h)
                make_thicker = (lambda >= 6);            % threshold — tune if needed
                pb_ratio = [2 1 1];                      % 3:1 width:height plot box
        
                % σx
                ax = nexttile(t,(i-1)*3+1);
                surf(ax, Xd, Yd, 0*sigx, sigx, 'EdgeColor','none');
                view(ax,2);
                if make_thicker
                    axis(ax,'tight');                    % allow visual thickening
                    pbaspect(ax,pb_ratio);
                else
                    axis(ax,'equal','tight');            % preserve geometry when not slender
                end
                obj.symClim(ax, sigx); grid(ax,'off'); obj.setCbar(ax,'\sigma_x');
                set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
                % σy
                ax = nexttile(t,(i-1)*3+2);
                surf(ax, Xd, Yd, 0*sigy, sigy, 'EdgeColor','none');
                view(ax,2);
                if make_thicker
                    axis(ax,'tight'); pbaspect(ax,pb_ratio);
                else
                    axis(ax,'equal','tight');
                end
                obj.symClim(ax, sigy); grid(ax,'off'); obj.setCbar(ax,'\sigma_y');
                set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
                % τxy
                ax = nexttile(t,(i-1)*3+3);
                surf(ax, Xd, Yd, 0*tauxy, tauxy, 'EdgeColor','none');
                view(ax,2);
                if make_thicker
                    axis(ax,'tight'); pbaspect(ax,pb_ratio);
                else
                    axis(ax,'equal','tight');
                end
                obj.symClim(ax, tauxy); grid(ax,'off'); obj.setCbar(ax,'\tau_{xy}');
                set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
                xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
                ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            end
        end

        function plot_disp_concentrated(obj)
            % Plot u and v (uy0 clamp) for concentrated load, side-by-side per plate
            plates = obj.plates;
            nP     = numel(plates);
        
            % ==== Figure A: v only, one row x two columns (Plate 1 vs Plate 2) ====
            figure('Name','v(x,0): vx0 vs uy0 — Plates 1 & 2');
            t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
            
            % -------- Tile 1: Plate 1 (v at y=0 vs x) --------
            [u1,v1,xv1,yv1] = plates(1).solve_plate_concentrated("vx0");
            [u2,v2,~,   ~ ] = plates(1).solve_plate_concentrated("uy0");
            [~,iy0] = min(abs(yv1));  % y = 0 index
            
            ax = nexttile(t,1);
            set(ax,'FontSize', max(12, obj.AX_LABEL_FS+2), 'FontName','Arial');  % <-- ticks size
            plot(ax, xv1, v1(iy0,:), 'LineWidth',2.0, 'DisplayName','v_{x}=0'); hold(ax,'on');
            plot(ax, xv1, v2(iy0,:), 'LineWidth',2.0, 'DisplayName','u_{y}=0');
            grid(ax,'off'); xlim(ax,[xv1(1) xv1(end)]);
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS+2);
            ylabel(ax,'v(x,0)','FontWeight','bold','FontSize',obj.AX_LABEL_FS+2);
            legend(ax,'Location','best','FontSize',obj.AX_LABEL_FS+2,'Box','on');
            
            % -------- Tile 2: Plate 1 (u along y at x = l) --------
            x0 = plates(1).l;  [~,ix0] = min(abs(xv1 - x0));
            
            ax = nexttile(t,2);
            set(ax,'FontSize', max(12, obj.AX_LABEL_FS+2), 'FontName','Arial');  % <-- ticks size
            plot(ax, u1(:,ix0), yv1, 'LineWidth',2.0, 'DisplayName','v_x=0'); hold(ax,'on');
            plot(ax, u2(:,ix0), yv1, 'LineWidth',2.0, 'DisplayName','u_y=0');
            grid(ax,'off'); ylim(ax,[yv1(1) yv1(end)]);
            xlabel(ax, sprintf('u(y, x=%.3g)', x0),'FontWeight','bold','FontSize',obj.AX_LABEL_FS+2);
            ylabel(ax, 'y', 'FontWeight','bold','FontSize',obj.AX_LABEL_FS+2);
            legend(ax,'Location','best','FontSize',obj.AX_LABEL_FS+2,'Box','on');
        end

        function plot_displacement_superposed(obj)
            plates = obj.plates;
            nP     = numel(plates);
        
            figure('Name','Deformed field (v and u)');
            t = tiledlayout(nP,2,'Padding','compact','TileSpacing','compact');
        
            for i = 1:nP
                [u,v,xv,yv] = plates(i).solve_superposed_zeroV();
                [X,Y] = meshgrid(xv,yv);
                mag = sqrt(u.^2 + v.^2);
                sf  = 0.1 * max(plates(i).l, plates(i).h) / max(mag(:)+eps);
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

        function plot_superposed(obj)
            % One plate: top row = u, v, rotation θ ; bottom row = σx, σy, τxy
            plate = obj.plates(2);
        
            figure('Name','Superposed response — deformations (top) & stresses (bottom)');
            t = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
        
            % --- Solve superposed response (exponential + concentrated, v(0,0)=0) ---
            [u,v,xv,yv,sigx,sigy,tauxy] = plate.solve_superposed_zeroV();
            [X,Y] = meshgrid(xv,yv);
        
            % Deformation scale
            mag = hypot(u,v);
            sf  = 0.07 * max(plate.l, plate.h) / max([mag(:); eps]);
            Xd  = X + sf*u;   Yd = Y + sf*v;
        
            % Slenderness-based aspect
            lambda = plate.l / plate.h;
            make_thicker = (lambda >= 6);
            pb_ratio = [3 1 1];
        
            % Symmetric limits for u,v
            uabs = max(abs(u(:)));
            vabs = max(abs(v(:)));
        
            % ---- NEW: small-rotation (rigid-body) angle θ = 1/2 (v_x - u_y) ----
            dx = xv(2)-xv(1);  dy = yv(2)-yv(1);
            [du_dy, du_dx] = gradient(u, dy, dx);   % returns [∂u/∂y, ∂u/∂x]
            [dv_dy, dv_dx] = gradient(v, dy, dx);   % returns [∂v/∂y, ∂v/∂x]
            theta = 0.5*(dv_dx - du_dy);            % radians
            tabs  = max(abs(theta(:)));
        
            % ---------- Top row: deformations ----------
            % u(x,y)
            ax = nexttile(t,1);
            surf(ax, Xd, Yd, zeros(size(Xd)), u, 'EdgeColor','none');
            view(ax,2); if make_thicker, axis(ax,'tight'); pbaspect(ax,pb_ratio); else, axis(ax,'equal','tight'); end
            obj.setCbar(ax,'u'); grid(ax,'off'); clim(ax,[-uabs uabs]);
            set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
            % v(x,y)
            ax = nexttile(t,2);
            surf(ax, Xd, Yd, zeros(size(Xd)), v, 'EdgeColor','none');
            view(ax,2); if make_thicker, axis(ax,'tight'); pbaspect(ax,pb_ratio); else, axis(ax,'equal','tight'); end
            obj.setCbar(ax,'v'); grid(ax,'off'); clim(ax,[-vabs vabs]);
            set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
            % θ(x,y) rotation
            ax = nexttile(t,3);
            surf(ax, Xd, Yd, zeros(size(Xd)), theta, 'EdgeColor','none');
            view(ax,2); if make_thicker, axis(ax,'tight'); pbaspect(ax,pb_ratio); else, axis(ax,'equal','tight'); end
            grid(ax,'off'); clim(ax,[-tabs tabs]);
            cb = colorbar(ax); cb.Label.String = '\theta (rad)'; cb.Label.FontSize = obj.CB_LABEL_FS;
            set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
            % ---------- Bottom row: stresses ----------
            % σx
            ax = nexttile(t,4);
            surf(ax, Xd, Yd, zeros(size(Xd)), sigx, 'EdgeColor','none');
            view(ax,2); if make_thicker, axis(ax,'tight'); pbaspect(ax,pb_ratio); else, axis(ax,'equal','tight'); end
            obj.symClim(ax, sigx); grid(ax,'off'); obj.setCbar(ax,'\sigma_x');
            set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
            % σy
            ax = nexttile(t,5);
            surf(ax, Xd, Yd, zeros(size(Xd)), sigy, 'EdgeColor','none');
            view(ax,2); if make_thicker, axis(ax,'tight'); pbaspect(ax,pb_ratio); else, axis(ax,'equal','tight'); end
            obj.symClim(ax, sigy); grid(ax,'off'); obj.setCbar(ax,'\sigma_y');
            set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        
            % τxy
            ax = nexttile(t,6);
            surf(ax, Xd, Yd, zeros(size(Xd)), tauxy, 'EdgeColor','none');
            view(ax,2); if make_thicker, axis(ax,'tight');
            pbaspect(ax,pb_ratio); else, axis(ax,'equal','tight'); end
            obj.symClim(ax, tauxy); grid(ax,'off');
            obj.setCbar(ax,'\tau_{xy}');
            set(ax,'FontSize',obj.CB_TICK_FS,'YDir','reverse');
            xlabel(ax,'x','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
            ylabel(ax,'y','FontWeight','bold','FontSize',obj.AX_LABEL_FS);
        end

        function plot_midlines(obj)
            plates = obj.plates;
            for i = 1:numel(plates)
                plates(i).l = 1;
            end
            nP     = numel(plates);
            
            figure('Name','midline_1');
            t = tiledlayout(2, nP,'Padding','compact','TileSpacing',...
                'compact');
            
            for i = 1:nP
                [up,vp,xv,yv]   = plates(i).solve_plate_k0();    
                [ub,vb,~,~]     = plates(i).solve_beam_k0();
                [~,iy0] = min(abs(yv));  % y = 0 index
                
                ax = nexttile(t,i);
                set(ax,'FontSize', max(12, obj.AX_LABEL_FS+2),...
                    'FontName','Arial');
                plot(ax, xv, vp(iy0,:), 'LineWidth',2.0,...
                    'DisplayName','plate'); hold(ax,'on');
                plot(ax, xv, vb(iy0,:), 'LineWidth',2.0,...
                    'DisplayName','beam');
                grid(ax,'off'); xlim(ax,[xv(1) xv(end)]);
                xlabel(ax,'x','FontWeight','bold','FontSize',...
                    obj.AX_LABEL_FS+2);
                ylabel(ax,'v(x,0)','FontWeight','bold','FontSize',...
                    obj.AX_LABEL_FS+2);
                legend(ax,'Location','best','FontSize',obj.AX_LABEL_FS+2,...
                    'Box','on');
                set(gca,'box','off','YDir','reverse');
                % title(2*plates(i).h/plates(i).l)
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
