classdef SpaceTimeSchwarzs < SpaceTimeMinkowski

    % Class members
    properties
        
        r_s = 1;
        
    end
    
    % Public class methods
    methods
        
        % Initialize
        function obj = SpaceTimeSchwarzs(varargin)
            
            obj@SpaceTimeMinkowski(varargin{:});
            obj.type = 'Schwarzschild';
            
        end

    end
    
    % Bequethed class methods
    methods (Access=protected)
        
        function [g__t_t, g__t_r,     g__t_theta,     g__t_phi, ...
                          g__r_r,     g__r_theta,     g__r_phi, ...
                                  g__theta_theta, g__theta_phi, ...
                                                    g__phi_phi] = metric(obj,r,theta)
            
            c = obj.c;
            r_s = obj.r_s;

            % Schwarzschild
            
            g__t_t = (1-r_s./r)*c^2;
            g__t_r = 0*r;
            g__t_theta = 0*r;
            g__t_phi = 0*r;
            g__r_r = -1./(1-r_s./r);
            g__r_theta = 0*r;
            g__r_phi = 0*r;
            g__theta_theta = -r.^2;
            g__theta_phi = 0*r;
            g__phi_phi = -r.^2.*sin(theta).^2;

        end

        % Compute Christoffel symbols
        function [Gamma_t,Gamma_r,Gamma_theta,Gamma_phi] = christoffel(obj,r,theta)
            
            c = obj.c;
            r_s = obj.r_s;
            
            % Schwarzschild

            % t
            Gamma_t__t_t = 0.*r;
            Gamma_t__t_r = r_s./(2*r.*(r-r_s));
            Gamma_t__t_theta = 0.*r;
            Gamma_t__t_phi = 0.*r;
            Gamma_t__r_r = 0.*r;
            Gamma_t__r_theta = 0.*r;
            Gamma_t__r_phi = 0.*r;
            Gamma_t__theta_theta = 0.*r;
            Gamma_t__theta_phi = 0.*r;
            Gamma_t__phi_phi = 0.*r;

            % r
            Gamma_r__t_t = c^2*r_s*(r-r_s)./(2*r.^3);
            Gamma_r__t_r = 0.*r;
            Gamma_r__t_theta = 0.*r;
            Gamma_r__t_phi = 0.*r;
            Gamma_r__r_r = -r_s./(2*r.*(r-r_s));
            Gamma_r__r_theta = 0.*r;
            Gamma_r__r_phi = 0.*r;
            Gamma_r__theta_theta = r_s-r;
            Gamma_r__theta_phi = 0.*r;
            Gamma_r__phi_phi = (r_s-r).*sin(theta).^2;

            % theta
            Gamma_theta__t_t = 0.*r;
            Gamma_theta__t_r = 0.*r;
            Gamma_theta__t_theta = 0.*r;
            Gamma_theta__t_phi = 0.*r;
            Gamma_theta__r_r = 0.*r;
            Gamma_theta__r_theta = 1./r;
            Gamma_theta__r_phi = 0.*r;
            Gamma_theta__theta_theta = 0.*r;
            Gamma_theta__theta_phi = 0.*r;
            Gamma_theta__phi_phi = -sin(theta).*cos(theta);

            % phi
            Gamma_phi__t_t = 0.*r;
            Gamma_phi__t_r = 0.*r;
            Gamma_phi__t_theta = 0.*r;
            Gamma_phi__t_phi = 0.*r;
            Gamma_phi__r_r = 0.*r;
            Gamma_phi__r_theta = 0.*r;
            Gamma_phi__r_phi = 1./r;
            Gamma_phi__theta_theta = 0.*r;
            Gamma_phi__theta_phi = cot(theta);
            Gamma_phi__phi_phi = 0.*r;

            Gamma_t = [Gamma_t__t_t     Gamma_t__t_r     Gamma_t__t_theta     Gamma_t__t_phi ...
                       Gamma_t__t_r     Gamma_t__r_r     Gamma_t__r_theta     Gamma_t__r_phi ...
                       Gamma_t__t_theta Gamma_t__r_theta Gamma_t__theta_theta Gamma_t__theta_phi ...
                       Gamma_t__t_phi   Gamma_t__r_phi   Gamma_t__theta_phi   Gamma_t__phi_phi];

            Gamma_r = [Gamma_r__t_t     Gamma_r__t_r     Gamma_r__t_theta     Gamma_r__t_phi ...
                       Gamma_r__t_r     Gamma_r__r_r     Gamma_r__r_theta     Gamma_r__r_phi ...
                       Gamma_r__t_theta Gamma_r__r_theta Gamma_r__theta_theta Gamma_r__theta_phi ...
                       Gamma_r__t_phi   Gamma_r__r_phi   Gamma_r__theta_phi   Gamma_r__phi_phi];

            Gamma_theta = [Gamma_theta__t_t     Gamma_theta__t_r     Gamma_theta__t_theta     Gamma_theta__t_phi ...
                           Gamma_theta__t_r     Gamma_theta__r_r     Gamma_theta__r_theta     Gamma_theta__r_phi ...
                           Gamma_theta__t_theta Gamma_theta__r_theta Gamma_theta__theta_theta Gamma_theta__theta_phi ...
                           Gamma_theta__t_phi   Gamma_theta__r_phi   Gamma_theta__theta_phi   Gamma_theta__phi_phi];

            Gamma_phi = [Gamma_phi__t_t     Gamma_phi__t_r     Gamma_phi__t_theta     Gamma_phi__t_phi ...
                         Gamma_phi__t_r     Gamma_phi__r_r     Gamma_phi__r_theta     Gamma_phi__r_phi ...
                         Gamma_phi__t_theta Gamma_phi__r_theta Gamma_phi__theta_theta Gamma_phi__theta_phi ...
                         Gamma_phi__t_phi   Gamma_phi__r_phi   Gamma_phi__theta_phi   Gamma_phi__phi_phi];
            
        end

    end
    
end