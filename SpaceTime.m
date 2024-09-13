classdef SpaceTime < handle
    
    % Bequeath to derived classes
    properties (SetAccess=protected)
        
        type = 'Unknown';
        v_tot;
        
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
            
            if (nargin>8 && ~isempty(varargin{1}))
                obj.v_tot = varargin{1};
            else
                obj.v_tot = obj.c;
            end
            if ( isempty(v_t) )
                % v_t is not specified at init and needs to be computed
                [g__t_t,           ~,              ~,            ~, ...
                              g__r_r,              ~,            ~, ...
                                      g__theta_theta,            ~, ...
                                                        g__phi_phi] = obj.metric(r,theta);
                v_t = sqrt( (obj.v_tot^2 - g__r_r.*v_r.*v_r - g__phi_phi.*v_phi.*v_phi - g__theta_theta.*v_theta.*v_theta)./g__t_t );
            end
            
            obj.y = [v_t v_r v_theta v_phi t r theta phi];
            
        end
        
        % Integrate derivative
        function obj = integrate(obj,h)
            
            k1 = obj.f(obj.y);
            k2 = obj.f(obj.y+0.5*h*k1);
            k3 = obj.f(obj.y+0.5*h*k2);
            k4 = obj.f(obj.y+h*k3);
            y = obj.y + h/6*(k1+2*k2+2*k3+k4);
            
            s = obj.y2states(y);
            [g__t_t,           ~,              ~,            ~, ...
                          g__r_r,              ~,            ~, ...
                                  g__theta_theta,            ~, ...
                                                    g__phi_phi] = obj.metric(s.r,s.theta);
            % vTot^2 = A*(vA*s)^2 + B*(vB*s)^2 + C*(vC*s)^2 + D*(vD*s)^2
            % vTot^2 = s^2(A*vA^2 + B*vB^2 + C*vC^2 + D*vD^2)
            % s^2(A*vA^2 + B*vB^2 + C*vC^2 + D*vD^2) - vTot^2 = 0
            v2 = g__t_t.*s.v_t.^2 + g__r_r.*s.v_r.^2 +  g__theta_theta.*s.v_theta.^2 + g__phi_phi.*s.v_phi.^2;
            s2a = 1-1e-9;
            s2b = 1+1e-9;
            v_tot2 = obj.v_tot^2;
            arr = [s2a*v2-v_tot2 s2b*v2-v_tot2];
            frac = -arr(:,1)./(arr(:,2)-arr(:,1));
            s2 = s2a + frac.*(s2b-s2a);
            idx = s2>0;
            s1 = sqrt(s2(idx));
            s.v_t(idx) = s1.*s.v_t(idx);
            s.v_r(idx) = s1.*s.v_r(idx);
            s.v_theta(idx) = s1.*s.v_theta(idx);
            s.v_phi(idx) = s1.*s.v_phi(idx);
            s.v_t(~idx) = nan;
            s.v_r(~idx) = nan;
            s.v_theta(~idx) = nan;
            s.v_phi(~idx) = nan;
            
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

            % y = [v_r, v_t, v_phi, v_theta, r, t, phi, theta]
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

        % Create strcut from states 
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