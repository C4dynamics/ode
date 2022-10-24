classdef rocket_dynamics < matlab.mixin.Copyable
    % class rocket_dynamics represents an airframe response of a rocket by the short period dynamics. 
    % input arguments:
    %   system type:
    %       1 open loop 
    %       2 gatr (damping feedback)
    %       3 pole placement controller 
    %       4 lqr 
    %       5 noisy open loop with full state estimation (lqe)
    %       6 lqg
    %       the state vector of each system is 
    %               |delta|                 wing deflection 
    %       x   =   |  q  |                 pitch rate 
    %               |alpha|                 angle of atack 
    %       and the air conditions are given: mach = 0.6, h = 1000ft
    %   x0  initial conditions (wrt xe). (degrees, degrees per second). 
    %   cmd (acceleration): 
    %       []          stability problem (initial condition)
    %       otherwise   tracking problem: 
    %           1 scalar    constant input (scalar gain)
    %           2 scalars   velocity input (y intersection and slope)
    %           3 scalars   acceleration input (quadratic polynomial coefficients)
    %   an  trim angle of attack about which to linearize the system. (deg)        
    % out arguments:
    %   y   full state vector 
    %
    % user object-functions: 
    %       set.mach
    %       set.xe
    %       set.x0
    %       run_sim
    %       drawstate
    % static functions 
    %       windtunnel
    %       phase_portait
    %%

    properties (Constant)
        hft     =   1000;           % ft            
        taud 	=	0.1;            % s             
        iyy  	=	1406.1689;      % kg*m^2        
        m    	=	453.592;        % kg            
        s    	=	0.073;          % m^2           
        d    	=	0.3048;         % m             
        g       =   9.8;            % m/s^2                 
        dt   	=	.5e-3;          % 0.1;  
    end
    
    properties 
        system
        ifig
        % parameters
%         mach   % --
        v = 914.4   % m/s             3000    ft/s   
        Qs          %      N/m^2 * m^2     

        % aero db
        DeltaDB
        AlphaDB 
        CmDB
        CzDB 
%         CaDB 
        Cmq 
        Cma
        Cza
        Cmd
        Czd
        za
        zd
         
        % state variables  
        t
        q
        alpha
        delta
        am
        x0 = [0; 0; 0]; % rad, rad/sec
        xe = [0; 0; 0]; % relevant only for the linearized system
        deltacmd = 0;
        
        % linear variables
        A
        b
    end
    
    
    methods 
        function obj = rocket_dynamics(system)%, mach

            if ~exist('system', 'var') || isempty(system) %
                system = 'nonlinear';
            end
%             if ~exist('mach', 'var') || isempty(mach) 
%                 mach = 3;
%             end

            obj.system = system;

            rho = 1.225;                    % kg/m^3. 
%             v_sound = obj.v / mach;
            Q = 1/2 * rho * obj.v^2; % 512128 N/m^2. finally it's pressure.
            obj.Qs = Q * obj.s;

            obj.loadaerodata;
            
            % find free figure number - without opening a window. 
            openfigs = findobj('Type', 'figure');
            if ~isempty(openfigs)
                obj.ifig = find(~ismember(1 : 100, openfigs), 1);
            else
                obj.ifig = 1;
            end
        end

        function set.xe(obj, val)
            obj.xe = val;
            obj.loadaerodata;
        end
        function set.x0(obj, val)
%            if strcmp(obj.system, 'linear')
%                val = obj.xe + val;
%            end
            obj.x0 = val;
        end
        

        function loadaerodata(obj)
            
            load('aerodynamic_coef')            
            
%             s6d.aero.ParasiticRollMethod = 1;
%             s6d = AeroDB(s6d);

            obj.DeltaDB = z_db.Delta;% -15 : 15; % [-15.0  -10.0   -7.5   -5.0   -2.5    0.    +2.5   +5.0   +7.5  +10.0   15.0];% s6d.aero.Delta2';% [-15.0  -10.0   -7.5   -5.0   -2.5    0.    +2.5   +5.0   +7.5  +10.0   15.0];
            obj.AlphaDB = z_db.Alpha; % 0 : 20;% s6d.aero.Alpha2';

            % [~, iM] = min(abs(s6d.aero.Mach2 - obj.mach));

%             DX = -0.0800162;
%             Dref = 80 / 1000; % diameter
%             DXD = DX / Dref;

%             Cmof     = squeeze(s6d.aero.CMOF(iM, :, :)); % r = mach. c = delta, d3 = alpha. 
            obj.CzDB = z_db.Cz; % squeeze(s6d.aero.CZOF(iM, :, :));
%             obj.CaDB = squeeze(s6d.aero.CAOF(iM, :, :));
            obj.CmDB = z_db.Cm; % Cmof + obj.CzDB * DXD; % dimensionless. a correction of the wind tunnel measures for the normal force contribution. 
            % sec/rad:                        1/rad                    m            m/s   (*q = rad/sec -> m*sec*rad/(rad*m*sec)                                       
            obj.Cmq = 0;% intrp1v(s6d.aero.Mach4, s6d.aero.CMQ, obj.mach) * obj.d / (2 * obj.v);
            
            %
            % stability derivatives 
            %%
            
            alphatrim = obj.xe(2) * 180 / pi;
            deltatrim = obj.xe(3) * 180 / pi;
            
            % interpolate for aoa in delta trim:
            cmdb_d0 = interp1(obj.DeltaDB, obj.CmDB, deltatrim);
            czdb_d0 = interp1(obj.DeltaDB, obj.CzDB, deltatrim);
            if isnan(cmdb_d0(1))
                warning('enter delta between -15 and 15 deg')
                return;
            end
            
            % make a linear approximation about alpha (delta = deltatrim).
            alphadb = (-20 : 0.1 : 20) * pi / 180;
            iaoa = find(alphadb >= alphatrim - 2 & alphadb <= alphatrim + 2);
            cmdb_d0al = interp1(obj.AlphaDB, cmdb_d0, alphadb);
            czdb_d0al = interp1(obj.AlphaDB, czdb_d0, alphadb);
            
            CMPa  = polyfit(alphadb(iaoa) * pi / 180, cmdb_d0al(iaoa), 1);    % dimensionless
            CZPa  = polyfit(alphadb(iaoa) * pi / 180, czdb_d0al(iaoa), 1);
           
            
            % interpolate for delta in aoa trim:
            cmdb_a0 = interp1(obj.AlphaDB, obj.CmDB', alphatrim);
            czdb_a0 = interp1(obj.AlphaDB, obj.CzDB', alphatrim);
            if isnan(cmdb_a0(1))
                warning('enter alpha between 0 and 40 deg')
                return;
            end
            
            % make a linear approximation about delta (aoa = alphatrim).
            deltadb = (-15 : 0.1 : 15) * pi / 180;
            id    = find(deltadb >= deltatrim - 2 & deltadb <= deltatrim + 2);
            cmdb_a0d = interp1(obj.DeltaDB, cmdb_a0, deltadb);
            czdb_a0d = interp1(obj.DeltaDB, czdb_a0, deltadb);
            
            CMPd  = polyfit(deltadb(id) * pi / 180, cmdb_a0d(id), 1);
            CZPd  = polyfit(deltadb(id) * pi / 180, czdb_a0d(id), 1);

            
            % drop x0
            obj.Cma = CMPa(1);
            obj.Cza = CZPa(1);
            obj.Cmd = CMPd(1);
            obj.Czd = CZPd(1);
        
        
             % build_linear_matrix           
            
            mq = obj.Qs * obj.d / obj.iyy * obj.Cmq; % sec/rad
            ma = obj.Qs * obj.d / obj.iyy * obj.Cma; % 1/rad
            md = obj.Qs * obj.d / obj.iyy * obj.Cmd; % 1/rad
            
%             san = sin(obj.xe(2)); % sin(alpha trim)
%             can = cos(obj.xe(2)); % cos(alpha trim)
%             mv = obj.m * obj.v;
%             Cxq = 0; Cxa = 0; Cxd = 0; CxXe = 0; 
%             Czq = 0;
%             CzXe = intrp2v(obj.AlphaDB, obj.DeltaDB, obj.CzDB, abs(obj.xe(2) * 180 / pi), obj.xe(3) * rocket_dynamics.copysign(obj.xe(2)) * 180 / pi) * rocket_dynamics.copysign(obj.xe(2)); 
%             z21 = 1 - obj.Qs * (can * Czq - san * Cxq) / mv;
%             z22 = -obj.Qs * (can * (obj.Cza - CxXe) - san * (CzXe + Cxa)) / mv;
%             z23 = -obj.Qs * (can * obj.Czd - san * Cxd) / mv;

            obj.za = obj.Qs * cos(obj.xe(2)) * obj.Cza / obj.m;
            obj.zd = obj.Qs * cos(obj.xe(2)) * obj.Czd / obj.m;
            
            obj.A = [mq     ma                  md              ...                  
                    ; 1    -obj.za / obj.v     -obj.zd / obj.v  ... 
                    ; 0     0                  -1 / obj.taud    ];
                
            obj.b =  [0;    0;      1 / obj.taud];
        end

        function x = run_sim(obj, tsim, theta0)

            if ~exist('tsim', 'var') || isempty(tsim) %
                tsim = 3;
            end
            if ~exist('theta0', 'var') || isempty(theta0) %
%                 theta0 = obj.x0(2);
                a0 = obj.xe(2);
                d0 = obj.xe(3);
%                 as = rocket_dynamics.copysign(a0);
                
                    % interp2(X,Y,V,Xq,Yq)
                Cz0 = interp2(obj.AlphaDB, obj.DeltaDB, obj.CzDB, a0, d0);                                            
                theta0 = a0 - acos(obj.Qs * Cz0 * cos(a0) / (obj.m * obj.g));
%                 theta0 = c_shortperiod.get_theta_trim(a0, d0)
            end

            obj.t = 0 : obj.dt : tsim;

            switch obj.system
                case 'nonlinear'
                    th0 = theta0; % -pi/2;% 
                    [~, x] = ode45(@(k1, k2) obj.dx_nonlinear(k1, k2), obj.t, [th0; obj.x0]); % 
                    
                    % x: 
                    %     1) theta, 
                    %     2) q, 
                    %     3) alpha, 
                    %     4) delta 
                    % extract am = gammadot / vm:
                    %       gammadot = -obj.g * cos(alpha - theta) / v + faz * cos(alpha) / (m * v)
                    %%
                    Cz = zeros(size(x, 1), 1);
                    for i = 1 : size(x, 1)
                        alpha_i = x(i, 3);
                        delta_i = x(i, 4);
%                         as = rocket_dynamics.copysign(alpha_i);
                        Cz(i) = interp2(obj.AlphaDB, obj.DeltaDB, obj.CzDB, alpha_i, delta_i);
                    end
%                     gamma0 = x(1, 1) - x(1, 3); % alpha0 - theta0
                    gammadot = -obj.g * cos(x(:, 3) - x(:, 1)) / obj.v ...
                                    + obj.Qs * Cz .* cos(x(:, 3)) / (obj.m * obj.v);
                    obj.am = obj.v * gammadot;
                    y = x(:, 2 : end);
                    % y: 
                    %     1) q, 
                    %     2) alpha, 
                    %     3) delta 
                case 'linear'
                    [~, y] = ode45(@(k1, k2) obj.dx_linear(k1, k2), obj.t, obj.x0); % 
                    alphadot = sum(obj.A(2, :) .* y, 2);
                    gammadot = y(:, 1) - alphadot; % thetadot - alphadot 
                    obj.am = gammadot * obj.v;
                    % y: 
                    %     1) q, 
                    %     2) alpha, 
                    %     3) delta 
                    y = obj.xe' + y;
                    % am = fN / m = -fzb * cos(alpha) / m
                    % linear fzb = Qs * (cza * z + czd * d)
%                     fzb = -obj.Qs * (obj.Cza * y(:, 2) + obj.Czd * y(:, 3));
%                     obj.am = -fzb .* cos(y(:, 2)) / obj.m;
%                   or:
%                     obj.am = obj.za .* y(:, 2) + obj.zd .* y(:, 3);
                    
                    x = [];
            end

            obj.q = y(:, 1);
            obj.alpha = y(:, 2);
            obj.delta = y(:, 3);            
        end

        function dx = dx_nonlinear(obj, t, x)
            % here an additional variable, theta, is necessary, since g is
            % oriented by it.
            % x = [theta; q; alpha; delta]
                        
            theta = x(1);
            qi = x(2);
            alpha_i = x(3);
            delta_i = x(4);
            
%         Cm(-alpha, delta) = -Cm(alpha, -delta)
%             as = rocket_dynamics.copysign(alpha_i);
            Cz = interp2(obj.AlphaDB, obj.DeltaDB, obj.CzDB, alpha_i, delta_i);                                            
            Cm = interp2(obj.AlphaDB, obj.DeltaDB, obj.CmDB, alpha_i, delta_i);% dimensionless
            Cx = 0;% intrp2v(obj.AlphaDB, obj.DeltaDB, obj.CaDB, abs(alpha_i * 180 / pi),                  0) * as;

            Cm = Cm + obj.Cmq * qi; % dimensionless + sec/rad * rad/sec = dimensionless

            dth = qi;
            dq = 1 / obj.iyy * obj.Qs * obj.d * Cm;
            da = qi + obj.g * cos(alpha_i - theta) / obj.v ...
                    + obj.Qs * (-Cz * cos(alpha_i) + Cx * sin(alpha_i)) / (obj.m * obj.v);
            dd = -1 / obj.taud * delta_i + 1 / obj.taud * obj.deltacmd;

            dx = [dth; dq; da; dd];
        end
        
        function drawstate(obj)
%             figure(obj.ifig)
%             close(obj.ifig)
            figure(obj.ifig)
            
            subplot(2, 2, 1)
            plot(obj.t, obj.delta * 180 / pi, 'm', 'linewidth', 2)
            hold on
            plot(obj.t, obj.deltacmd * ones(size(obj.t)) * 180 / pi, 'c', 'linewidth', 1.5)
            legend('\delta_{measure}', '\delta_{command}')
            rocket_dynamics.plotdefaults(gca, 'Fin Deflection', 't', 'deg')
            hold off
            y1 = ylim;
            ylim([y1(1) - 1, y1(2) + 1])

            subplot(2, 2, 3)
            plot(obj.t, obj.q * 180 / pi, 'm', 'linewidth', 2)
            rocket_dynamics.plotdefaults(gca, 'Pitch Rate', 't', 'deg/sec')
            y1 = ylim;
            ylim([y1(1) - 1, y1(2) + 1])

            subplot(2, 2, 2)
            plot(obj.t, obj.alpha * 180 / pi, 'm', 'linewidth', 2)
            rocket_dynamics.plotdefaults(gca, 'Angle of Attack', 't', 'deg')
            y1 = ylim;
            ylim([y1(1) - 1, y1(2) + 1])
            
            subplot(2, 2, 4)
            plot(obj.t, obj.am / obj.g, 'm', 'linewidth', 2)
            rocket_dynamics.plotdefaults(gca, 'Acceleration', 't', 'g')
            y1 = ylim;
            ylim([y1(1) - 1, y1(2) + 1])
        end
        
        function ranimate(obj)
            figure('color', [0.4 0 0.4])
%             axis off
            theta = cumtrapz(obj.q) * obj.dt;
            filename = 'rocketsim.gif';
            lr =  3 / 4; 
            bklr = 1 / 4;
            displaylim = [-0.5 1.5 -0.8 0.8];
            rtail = [0 0];
            ax = gca;
            set(ax, 'visible', 'off');
            
            for t = 1 : 10 : length(theta) %400% 
                plot(displaylim(1 : 2), [0 0], 'm', 'LineWidth', 1)
                hold(ax, 'on')
                set(ax, 'visible', 'off');
                
                %
                % edges
                %% 
                rhead = lr * [cos(theta(t)) sin(theta(t))];
                theta_tail = theta(t) + pi / 2;
                rtail2 = bklr * [-sin(theta_tail) cos(theta_tail)];

                % 
                % canard
                %%
                rcanard = lr * 3 / 4;
                canard_center = rcanard * [cos(theta(t)) sin(theta(t))];
                th_canardl = theta(t) + 20 * pi / 180;
                th_canardr = theta(t) - 20 * pi / 180;
                canard_left = rcanard * 1.05 * [cos(th_canardl) sin(th_canardl)];
                canard_right = rcanard * 1.05 * [cos(th_canardr) sin(th_canardr)];
                plot(ax, [canard_center(1), canard_left(1)], [canard_center(2), canard_left(2)], 'g', 'LineWidth', 5)
                plot(ax, [canard_center(1), canard_right(1)], [canard_center(2), canard_right(2)], 'g', 'LineWidth', 5)
                 
                % 
                % tail fins
                %%
                rfins = bklr * 4 / 5;
                fins_center = rfins * [-sin(theta_tail) cos(theta_tail)];
                th_finsl = theta_tail + 45 * pi / 180;
                th_finsr = theta_tail - 45 * pi / 180;
                fins_left = rfins * 2 * [-sin(th_finsl) cos(th_finsl)];
                fins_right = rfins * 2 * [-sin(th_finsr) cos(th_finsr)];
                plot(ax, [fins_center(1), fins_left(1)], [fins_center(2), fins_left(2)], 'b', 'LineWidth', 12)
                plot(ax, [fins_center(1), fins_right(1)], [fins_center(2), fins_right(2)], 'b', 'LineWidth', 12)

                %
                % edges
                %%
                plot(ax, [rtail(1) rhead(1)], [rtail(2) rhead(2)], 'k', 'LineWidth', 50, 'Color',[0.960784316062927 0.921568632125854 0.921568632125854])
                plot(ax, rhead(1), rhead(2), 'ok', 'LineWidth', 10, 'MarkerFaceColor', [1 0.5 0], 'MarkerSize', 40, 'Color',[0.960784316062927 0.921568632125854 0.921568632125854]);
                plot(ax, [rtail2(1) rtail(1)], [rtail2(2) rtail(2)], 'k', 'LineWidth', 50, 'Color',[0.960784316062927 0.921568632125854 0.921568632125854])

                
                axis(displaylim);
                drawnow
                
                %%
                  frame = getframe(1);
                  im = frame2im(frame);
                  [imind, cm] = rgb2ind(im, 256);
                  if t == 1;
                      imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
                  else
                      imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
                  end
                %%
                
                hold(ax, 'off')
%                 hold off
            end
            
        end
    end

    methods (Static)
        function draw_db()
            rocket = rocket_dynamics();
            
            f1 = figure;
            hold on 
            f2 = figure;
            hold on 

            for d = 1 : 3 : length(rocket.DeltaDB)
                c = rand(1, 3);
                figure(f1)
                plot(rocket.AlphaDB * 180 / pi, rocket.CzDB(d, :), 'color', c)
                figure(f2)
                plot(rocket.AlphaDB * 180 / pi, rocket.CmDB(d, :), 'color', c)
            end

            figure(f1)
            hold off
            grid
            box
            title('Normal Force Coefficient', 'fontname', 'times', 'fontsize', 12)
            xlabel('angle of attack (deg)', 'fontname', 'times', 'fontsize', 12)
            ylabel('C_N', 'fontname', 'times', 'fontsize', 12)
            lg = legend(strsplit(num2str(rocket.DeltaDB(1 : 3 : end)*180/pi)), 'location', 'best');
            lg_ttl = get(lg, 'title');
            set(lg_ttl, 'string', '\delta')

            figure(f2)
            hold off
            grid
            box
            title('Pitch Moment Coefficient', 'fontname', 'times', 'fontsize', 12)
            xlabel('angle of attack (deg)', 'fontname', 'times', 'fontsize', 12)
            ylabel('C_M', 'fontname', 'times', 'fontsize', 12)
            lg = legend(strsplit(num2str(rocket.DeltaDB(1 : 3 : end)*180/pi)), 'location', 'best');
            lg_ttl = get(lg, 'title');
            set(lg_ttl, 'string', '\delta')


        end
        function signed_val = copysign(xsign, xval)
        % function copysign(xval, xsign) cretaes a parameter by the sign of one 
        %   parameter, xsign, and the value of another, xval
        % xsign zero is regarded positive 
        % xval empty is regarded one 
        %%

            if nargin < 2
                xval = 1;
            end

            xsign0 = sign(xsign);
            xsign0(xsign0 == 0) = 1;
            signed_val = xsign0 * abs(xval);

        end
        function plotdefaults(ax, ptitle, pxlabel, pylabel)
            title(ptitle, 'fontname', 'times', 'fontsize', 12)
            xlabel(pxlabel, 'fontname', 'times', 'fontsize', 12)
            ylabel(pylabel, 'fontname', 'times', 'fontsize', 12)
            grid on
            box on
            set(ax, 'fontname', 'times', 'fontsize', 12)
        end
    end
end




