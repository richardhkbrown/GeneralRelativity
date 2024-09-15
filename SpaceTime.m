classdef SpaceTime < handle
    
    % Bequeath to derived classes
    properties (SetAccess=protected)
        
        type = 'Unknown';
        v_tot;
        t_correction = true;
        hasMass = false;
        
    end
    
    % Unstored
    properties (Dependent)
        
        r_max;
        G;
        c;
        
    end
    
    % Bequeath to derived classes
    properties (SetAccess=private,GetAccess=public)
                
        % y = [v_t, v_r, v_theta, v_phi, t, r, theta, phi]
        y;
        
    end

    % Same for all objects
    methods (Static)
        
        % Way to get static parameters
        function out = fakeStatic(data)
          
            persistent var;
            if ( isempty(var) )
                var.r_max = 10;
                var.G = 1;
                var.c = 10;
            end
            if ( nargin )
                var = data;
            end
            out = var;
            
        end
        
    end
    
    % Unique for this class
    methods
        
        % Initialize
        function obj = SpaceTime(t,r,theta,phi,v_t,v_r,v_theta,v_phi,varargin)
            
            if ( nargin>8 && ~isempty(varargin{1}) )
                obj.v_tot = varargin{1};
            else
                obj.v_tot = obj.c;
            end
            if ( nargin>9 && ~isempty(varargin{2}) )
                obj.type = varargin{2};
            end
            if ( isempty(v_t) )
                v_t = nan*v_r;
            end
            
            obj.y = [v_t v_r v_theta v_phi t r theta phi];
            
            obj.hasMass = any(strcmp(fieldnames(obj),'r_s'));
            
        end
        
        % Integrate derivative
        function obj = integrate(obj,h)
            
            idxVt = isnan(obj.y(:,1)) & ~isnan(obj.y(:,2));
            if ( any(idxVt) )
                
                s = obj.y2states(obj.y);
                
                % v_t is not specified at init and needs to be computed
                [g__t_t,           ~,              ~,            ~, ...
                              g__r_r,              ~,            ~, ...
                                      g__theta_theta,            ~, ...
                                                        g__phi_phi] = obj.metric(s.r,s.theta);
                % vTot^2 = A*vA^2 + B*vB^2 + C*vC^2 + D*vD^2
                % vTot^2 - B*vB^2 - C*vC^2 - D*vD^2 =  A*vA^2
                % vA^2 = (vTot^2 - B*vB^2 - C*vC^2 - D*vD^2) / A
                v_t2 = (obj.v_tot^2 - g__r_r.*s.v_r.^2 - g__phi_phi.*s.v_phi.^2 - g__theta_theta.*s.v_theta.^2) ./ g__t_t;
                s.v_t(idxVt) = sqrt( v_t2(idxVt) );
                obj.y = obj.states2y(s);
                
            end
                
            k1 = obj.f(obj.y);
            k2 = obj.f(obj.y+0.5*h*k1);
            k3 = obj.f(obj.y+0.5*h*k2);
            k4 = obj.f(obj.y+h*k3);
            y = obj.y + h/6*(k1+2*k2+2*k3+k4);
            
            s = obj.y2states(y);
            if ( obj.hasMass )
                s.r(s.r < obj.r_s) = obj.r_s; % clip to Schwartzchild radius r_s
            end
            [g__t_t,           ~,              ~,            ~, ...
                          g__r_r,              ~,            ~, ...
                                  g__theta_theta,            ~, ...
                                                    g__phi_phi] = obj.metric(s.r,s.theta);
            if ( obj.t_correction )
                % vTot^2 = A*(vA*s)^2 + B*vB^2 + C*vC^2 + D*vD^2
                % vTot^2 - B*vB^2 - C*vC^2 - D*vD^2 = s^2(A*vA^2)
                % s^2 = (vTot^2 - B*vB^2 - C*vC^2 - D*vD^2) / A*vA^2 
                s2 = (obj.v_tot^2 - g__r_r.*s.v_r.^2 - g__theta_theta.*s.v_theta.^2 - g__phi_phi.*s.v_phi.^2) ./ (g__t_t.*s.v_t.^2);
                s1 = sqrt(s2);
                if ( ~isreal(s1) )
                    disp('what!');
                end
                obj.t_correction = false;
            else
                % vTot^2 = A*vA^2 + B*(s*vB)^2 + C*(s*vC)^2 + D*(s*vD)^2
                % vTot^2 - A*vA^2 = s^2(B*vB^2 + C*vC^2 + D*vD^2)
                % s^2 = (vTot^2 - A*vA^2 ) / (B*vB^2 + C*vC^2 + D*vD^2) 
                s2 = (obj.v_tot^2 - g__t_t.*s.v_t.^2) ./ (g__r_r.*s.v_r.^2 + g__theta_theta.*s.v_theta.^2 + g__phi_phi.*s.v_phi.^2);
                s1 = sqrt(s2);
                if ( ~isreal(s1) )
                    disp('wha!');
                end
                obj.t_correction = true;
            end
            s.v_t = s1.*s.v_t;
            s.v_r = s1.*s.v_r;
            s.v_theta = s1.*s.v_theta;
            s.v_phi = s1.*s.v_phi;

            obj.y = obj.states2y(s);
 
        end

        % Mutators
        function set.type(obj,val)
             obj.type = val;
        end
        function set.v_tot(obj,val)
             obj.v_tot = val;
        end
        
        % Mutators and Accessores for constants
        function set.r_max(obj,value)
            temp = obj.fakeStatic;
            temp.r_max = value;
            obj.fakeStatic(temp)
        end        
        function out = get.r_max(obj)
            out = obj.fakeStatic.r_max;
        end
        function set.G(obj,value)
            temp = obj.fakeStatic;
            temp.G = value;
            obj.fakeStatic(temp)
        end        
        function out = get.G(obj)
            out = obj.fakeStatic.G;
        end
        function set.c(obj,value)
            temp = obj.fakeStatic;
            temp.c = value;
            obj.fakeStatic(temp)
        end        
        function out = get.c(obj)
            out = obj.fakeStatic.c;
        end
        
        % Compute derivative
        function dY = f(obj,y)

            s = SpaceTime.y2states(y);

            v_v = [    s.v_t.*s.v_t     s.v_t.*s.v_r     s.v_t.*s.v_theta     s.v_t.*s.v_phi ...
                       s.v_r.*s.v_t     s.v_r.*s.v_r     s.v_r.*s.v_theta     s.v_r.*s.v_phi ...
                   s.v_theta.*s.v_t s.v_theta.*s.v_r s.v_theta.*s.v_theta s.v_theta.*s.v_phi ...
                     s.v_phi.*s.v_t   s.v_phi.*s.v_r   s.v_phi.*s.v_theta   s.v_phi.*s.v_phi]';

            [Gamma_t,Gamma_r,Gamma_theta,Gamma_phi] = obj.christoffel(s.r,s.theta);
            dvt_dlamda = -sum(Gamma_t.*v_v',2);
            dvr_dlamda = -sum(Gamma_r.*v_v',2);
            dvtheta_dlamda = -sum(Gamma_theta.*v_v',2);
            dvphi_dlamda = -sum(Gamma_phi.*v_v',2);

            dY = [dvt_dlamda dvr_dlamda dvtheta_dlamda dvphi_dlamda s.v_t s.v_r s.v_theta s.v_phi];

        end
        
    end
    
    % No member variables
    methods (Static, Access = public)

        % Create states from y 
        function s = y2states(y)

            s.v_t = y(:,1);
            s.v_r = y(:,2);
            s.v_theta = y(:,3);
            s.v_phi = y(:,4);
            s.t = y(:,5);
            s.r = y(:,6);
            s.theta = y(:,7);
            s.phi = y(:,8);
              
        end
        
        function y = states2y(s)

            % y = [v_t, v_r, v_theta, v_phi, t, y, theta, phi]
            y = [s.v_t s.v_r s.v_theta s.v_phi s.t s.r s.theta s.phi];
              
        end
        
    end
        
    % No member variables
    methods (Abstract, Static, Access=protected)
        
        [g__t_t, g__t_r,     g__t_theta,     g__t_phi, ...
                          g__r_r,     g__r_theta,     g__r_phi, ...
                                  g__theta_theta, g__theta_phi, ...
                                                    g__phi_phi] = metric(r,theta)

        % Compute Christoffel symbols
        [Gamma_t,Gamma_r,Gamma_theta,Gamma_phi] = christoffel(r,theta)

    end
   
end