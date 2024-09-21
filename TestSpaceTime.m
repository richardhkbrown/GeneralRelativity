close all;
clear all;

% set empty for time, set 0 for light
v_tot = [];

if ( isempty(v_tot) )
    r = [8 10 12]';
else
    r = [1 1.0666]';
end
t = 0*ones(size(r));
theta = pi/2*ones(size(r));
phi = 0*ones(size(r));
v_r = 0*ones(size(r));
v_theta = 0*ones(size(r));
v = 0.1;
v_phi = v./r;

% Minkowski
minkST = SpaceTimeMinkowski(t,r,theta,phi,[],v_r,v_theta,v_phi,v_tot,'Mnk');

mink.y = minkST.y;
mink.s = SpaceTime.y2states(mink.y);
x = mink.s.r.*cos(mink.s.phi);
y = mink.s.r.*sin(mink.s.phi);
mink.H = plot([x x]',[y y]','bo-');
axis equal;

% Schwarzschild
if ( isempty(v_tot) )
    r_s_s = [0 .002 0.01];
    a_s = [-20 0 20];
    iLoop = 200000;
    h = 0.002;
else
    % p_s = 3 * r_s / 2
    % r_s = 2 * p_s / 3
    r_s_s = [0 1*r(1)/3 2*r(1)/3];
    a_s = [-0.1 0 0.1];
    iLoop = 200000;
    h = 0.002;
end
schw = [];
hold on;
for iN = 1:length(r_s_s)
    schwST(iN) = SpaceTimeSchwarzs(t,r,theta,phi,[],v_r,v_theta,v_phi,v_tot,sprintf('Sch%02d',iN));
    schwST(iN).r_s = r_s_s(iN);
    schw(iN).y = schwST(iN).y;
    schw(iN).s = SpaceTime.y2states(schw(iN).y);
    schw(iN).H = plot([x x]',[y y]','g+-');
end
hold off;

% Kerr
kerr = [];
hold on;
for iN = 1:length(a_s)
    kerrST(iN) = SpaceTimeKerr(t,r,theta,phi,[],v_r,v_theta,v_phi,v_tot,sprintf('Ker%02d',iN));
    if ( isempty(v_tot) )
        kerrST(iN).r_s = r_s_s(2);
    else
        kerrST(iN).r_s = r_s_s(3);
    end
    kerrST(iN).a = a_s(iN);
    kerr(iN).y = kerrST(iN).y;
    kerr(iN).s = SpaceTime.y2states(kerr(iN).y);
    kerr(iN).H = plot([x x]',[y y]','r.');
end
hold off;

iSkip = 1000;
for iLoop = 1:iLoop
    
    minkST = minkST.integrate(h);
    for iN = 1:length(r_s_s)
        schwST(iN) = schwST(iN).integrate(h);
    end
    for iN = 1:length(a_s)
        kerrST(iN) = kerrST(iN).integrate(h);
    end    
    
    if ( mod(iLoop,iSkip)==1 )
        mink.s = [SpaceTime.y2states(minkST.y)];     
        x = [vertcat(mink.H.XData) mink.s.r.*cos(mink.s.phi)];
        y = [vertcat(mink.H.YData) mink.s.r.*sin(mink.s.phi)];
        x = mat2cell(x,ones(1,size(x,1)),size(x,2));
        y = mat2cell(y,ones(1,size(y,1)),size(y,2));
        [mink.H.XData] = deal(x{:});
        [mink.H.YData] = deal(y{:});
        for iN = 1:length(r_s_s)
            schw(iN).s = SpaceTime.y2states(schwST(iN).y);
            x = [vertcat(schw(iN).H.XData) schw(iN).s.r.*cos(schw(iN).s.phi)];
            y = [vertcat(schw(iN).H.YData) schw(iN).s.r.*sin(schw(iN).s.phi)];
            x = mat2cell(x,ones(1,size(x,1)),size(x,2));
            y = mat2cell(y,ones(1,size(y,1)),size(y,2));
            [schw(iN).H.XData] = deal(x{:});
            [schw(iN).H.YData] = deal(y{:});
        end
        for iN = 1:length(a_s)
            kerr(iN).s = SpaceTime.y2states(kerrST(iN).y);
            x = [vertcat(kerr(iN).H.XData) kerr(iN).s.r.*cos(kerr(iN).s.phi)];
            y = [vertcat(kerr(iN).H.YData) kerr(iN).s.r.*sin(kerr(iN).s.phi)];
            x = mat2cell(x,ones(1,size(x,1)),size(x,2));
            y = mat2cell(y,ones(1,size(y,1)),size(y,2));
            [kerr(iN).H.XData] = deal(x{:});
            [kerr(iN).H.YData] = deal(y{:});
        end
        drawnow;
    end

end
