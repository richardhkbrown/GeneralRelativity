classdef SpaceTimeKerr < SpaceTimeSchwarzs
    
    % Class members
    properties
        
        a = 1;
        
    end
    
    % Not shared
    methods
        
        % Initialize
        function obj = SpaceTimeKerr(varargin)

            if ( nargin<10 )
                varargin{10} = 'Kerr';
            end
            obj@SpaceTimeSchwarzs(varargin{:});
            
        end
        
    end
    
    % Not shared and not public
    methods (Access=protected)
        
        function [g__t_t, g__t_r,     g__t_theta,     g__t_phi, ...
                          g__r_r,     g__r_theta,     g__r_phi, ...
                                  g__theta_theta, g__theta_phi, ...
                                                    g__phi_phi] = metric(obj,r,theta)
            
            c = obj.c;
            r_s = obj.r_s;
            a = obj.a;
            Sigma = r.^2 + a^2*cos(theta).^2;
            Delta = r.^2 - r_s*r + a^2;
            
            % Kerr
            
            g__t_t = (1-r_s*r./Sigma)*c^2;
            g__t_r = 0*r;
            g__t_theta = 0*r;
            g__t_phi = 2*r_s*a*r.*sin(theta).^2*c./Sigma;
            g__r_r = -Sigma./Delta;
            g__r_theta = 0*r;
            g__r_phi = 0*r;
            g__theta_theta = -Sigma;
            g__theta_phi = 0*r;
            g__phi_phi = -(r.^2 + a^2 + r_s*a^2*r.*sin(theta).^2./Sigma).*sin(theta).^2;

        end

        % Compute Christoffel symbols
        function [Gamma_t,Gamma_r,Gamma_theta,Gamma_phi] = christoffel(obj,r,theta)
            
            c = obj.c;
            r_s = obj.r_s;
            a = obj.a;
            Sigma = r.^2 + a^2*cos(theta).^2;
            Delta = r.^2 - r_s*r + a^2;
            mu = r.^2 + a^2;
            A = (r.^2+a^2).^2-a.^2.*Delta.*sin(theta).^2;
    
            % Kerr

            % t
            Gamma_t__t_t = 0*r;
            Gamma_t__t_r = r_s*(r.^2+a^2).*(r.^2-a^2*cos(theta).^2)./(2*Sigma.^2.*Delta);
            Gamma_t__t_theta = -r_s*a^2*r.*sin(theta).*cos(theta)./Sigma.^2;
            Gamma_t__t_phi = 0*r;
            Gamma_t__r_r = 0*r;
            Gamma_t__r_theta = 0*r;
            Gamma_t__r_phi = r_s*a*sin(theta).^2.*(a^2*cos(theta).^2.*(a^2-r.^2)-r.^2.*(a^2+3*r.^2))./(2*c*Sigma.^2.*Delta);
            Gamma_t__theta_theta = 0*r;
            Gamma_t__theta_phi = r_s*a^3*r.*sin(theta).^3.*cos(theta)./(c*Sigma.^2);
            Gamma_t__phi_phi = 0*r;

            % r
            Gamma_r__t_t = c^2*r_s*Delta.*(r.^2-a^2*cos(theta).^2)./(2*Sigma.^3);
            Gamma_r__t_r = 0*r;
            Gamma_r__t_theta = 0*r;
            Gamma_r__t_phi = -c*Delta*r_s*a.*sin(theta).^2.*(r.^2-a^2.*cos(theta).^2)./(2.*Sigma.^3);
            Gamma_r__r_r = (2*r*a^2.*sin(theta).^2-r_s*(r.^2-a^2*cos(theta).^2))./(2.*Sigma.*Delta);
            Gamma_r__r_theta = -a^2*sin(theta).*cos(theta)./Sigma;
            Gamma_r__r_phi = 0*r;
            Gamma_r__theta_theta = -r.*Delta./Sigma;
            Gamma_r__theta_phi = 0*r;
            Gamma_r__phi_phi = Delta.*sin(theta).^2./(2*Sigma.^3).*(-2*r.*Sigma.^2+r_s*a.^2.*sin(theta).^2.*(r.^2-a^2*cos(theta).^2));

            % theta
            Gamma_theta__t_t = -c^2*r_s*a^2*sin(theta).*cos(theta)./Sigma.^3;
            Gamma_theta__t_r = 0*r;
            Gamma_theta__t_theta = 0*r;
            Gamma_theta__t_phi = c*r_s*a*r.*(r.^2+a^2).*sin(theta).*cos(theta)./Sigma.^3;
            Gamma_theta__t_phi = r_s*a*r.*(r.^2+a^2).*sin(theta).*cos(theta)./Sigma.^3;
            Gamma_theta__r_r = a^2*sin(theta).*cos(theta)./(Sigma.*Delta);
            Gamma_theta__r_theta = r./Sigma;
            Gamma_theta__r_phi = 0*r;
            Gamma_theta__theta_theta = -a^2*sin(theta).*cos(theta)./Sigma;
            Gamma_theta__theta_phi = 0*r;
            Gamma_theta__phi_phi = -sin(theta).*cos(theta)./Sigma.^3.*(A.*Sigma+(r.^2+a^2).*r_s*a^2.*r.*sin(theta).^2);

            % phi
            Gamma_phi__t_t = 0*r;
            Gamma_phi__t_r = c*r_s*a.*(r.^2-a^2*cos(theta).^2)./(2*Sigma.^2.*Delta);
            Gamma_phi__t_theta = -c*r_s*a*r.*cot(theta)./Sigma.^2;
            Gamma_phi__t_phi = 0*r;
            Gamma_phi__r_r = 0*r;
            Gamma_phi__r_theta = 0*r;
            Gamma_phi__r_phi = (2*r.*Sigma.^2+r_s*(a^4*sin(theta).^2.*cos(theta).^2-r.^2.*(Sigma+r.^2+a^2)))./(2*Sigma.^2.*Delta);
            Gamma_phi__theta_theta = 0*r;
            Gamma_phi__theta_phi = cot(theta)./Sigma.^2.*(Sigma.^2+r_s*a^2*r.*sin(theta).^2);
            Gamma_phi__phi_phi = 0*r;

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