function rmr()
    %% Constant definitios.
    % Link lengths
    l(1) = 1.0;
    l(2) = 1.0;
    l(3) = 1.1;
    l(4) = 1.1;
    % Sampling period
    ts = 0.05;
    % Redundant solution gain
    kh = 0.05;
    
    %% Desired end-effector line equation (ax + by + c = 0).
    a = 1.0;
    b = 0;
    c = -3.0;
    
    %% Initial joint variable values:
    % q1 = -60 deg, q2 = 0 deg, q3 = 50 deg, q4 = 20 deg.
    q = zeros([4 1]);
    q(1) = -pi / 3.0;
    q(2) = 0;
    q(3) = pi / 3.6;
    q(4) = pi / 9.0;

    %% Obstacle related definitions
    % Obstacle radius
    obr = 0.2;
    % Distance between obstacles (y-direction).
    obdist = 1.0;
    % Initial obstacle positions.
    % ob(1) = x-coordinate of obstacle 1.
    % ob(2) = y-coordinate of obstacle 1.
    % ob(3) = x-coordinate of obstacle 2.
    % ob(4) = x-coordinate of obstacle 2.
    ob(1) = 1.2;
    ob(2) = -2.2;
    ob(3) = 1.2;
    ob(4) = ob(2) + obdist;

    %% Initialize figure and buttons.
    fhandle = figure('Position', [200, 200, 800, 600],...
                     'Resize', 'off');
    upbutton = uicontrol(fhandle,...
                         'Style', 'pushbutton',...
                         'String', 'Up',...
                         'Position', [20, 100, 60, 20]);
    set(upbutton, 'Callback', 'upbutton_pressed()');
    downbutton = uicontrol(fhandle,...
                           'Style', 'pushbutton',...
                           'String', 'Down',...
                           'Position', [20, 60, 60, 20]);
    closebutton = uicontrol(fhandle,...
                            'Style', 'pushbutton',...
                            'String', 'Close',...
                            'Position', [20, 20, 60, 20]);
    
    %% Initialize plot.
    axeshandle = axes('Parent', fhandle,...
                      'Position', [0.15, 0.05, 0.80, 0.90]);
    axis(axeshandle, [-1 5 -4 2]);
    
    
    %% Resolved motion rate control loop.
    while true
        % Clear old plot.
        cla(axeshandle);
        % Draw manipulator.
        draw(q, ob);
        pause(ts);
        % Compute current px.
        px = l(1) * cos(q(1)) +...
             l(2) * cos(q(1) + q(2)) +...
             l(3) * cos(q(1) + q(2) + q(3)) +...
             l(4) * cos(q(1) + q(2) + q(3) + q(4));
        % (Desired speed) = ((desired px) - (current px)) / ts
        % XXX: This only works for a vertical desired line for now.
        pdot = (-c - px) / ts;
        % Find redundant inverse solution.
        J = jac(q);
        PJ = pinv(J);
        qdot = PJ * pdot + kh * (eye([4 4]) - PJ * J) * hgrad(q, ob);
        q = q + qdot * ts;
    end
    
    %% Gradient of the cost function regarding the obstacles.
    function h = hgrad(q, ob)
       h = [0 0 0 0]';
    end
    
    %% Jacobian matrix.
     % Just a 1x4 vector in our case, because the desired task only
     % requires 1 dof.
    function j = jac(q)
        j1 = - l(1) * sin(q(1)) -...
               l(2) * sin(q(1) + q(2)) -...
               l(3) * sin(q(1) + q(2) + q(3)) -...
               l(4) * sin(q(1) + q(2) + q(3) + q(4));
        j2 = - l(2) * sin(q(1) + q(2)) -...
               l(3) * sin(q(1) + q(2) + q(3)) -...
               l(4) * sin(q(1) + q(2) + q(3) + q(4));
        j3 = - l(3) * sin(q(1) + q(2) + q(3)) -...
               l(4) * sin(q(1) + q(2) + q(3) + q(4));
        j4 = - l(4) * sin(q(1) + q(2) + q(3) + q(4));
        j = [j1 j2 j3 j4];
    end
    
    %% Drawing function.
    function draw(q, ob)
        % Constants (colors, sizes).
        eeEdgeColor = [0.8 0 0];
        eeFaceColor = [0.8 0 0];
        jEdgeColor = [0 0 1.0];
        jFaceColor = [0 0 1.0];
        obEdgeColor = [0 0.8 0];
        obFaceColor = [0 0.8 0];
        lWidth = 3;
        eWidth = 2;
        jSize = 10;
        obSize = 10 * obr;
        % Draw the desired end-effector position line (crashes for
        % horizontal lines).
        line([(4*b-c)/a, (-4*b-c)/a], [-3.5, 1.5],...
            'LineWidth', eWidth,...
            'LineStyle', '--');
        set(axeshandle, 'NextPlot', 'add');
        % Draw first joint.
        plot(axeshandle, 0, 0, 'o',...
            'MarkerSize', 15,...
            'MarkerEdgeColor', jEdgeColor,...
            'MarkerFaceColor', jFaceColor);
        % Draw first link.
        sx = 0;
        sy = 0;
        ex = l(1)*cos(q(1));
        ey = l(1)*sin(q(1));
        line([sx, ex], [sy, ey], 'LineWidth', lWidth);
        % Draw second joint.
        plot(axeshandle, ex, ey, 'o',...
            'MarkerSize', jSize,...
            'MarkerEdgeColor', jEdgeColor,...
            'MarkerFaceColor', jFaceColor);
        % Draw second link.
        sx = ex;
        sy = ey;
        ex = sx + l(2)*cos(q(1) + q(2));
        ey = sy + l(2)*sin(q(1) + q(2));
        line([sx, ex], [sy, ey], 'LineWidth', lWidth);
        % Draw third joint.
        plot(axeshandle, ex, ey, 'o',...
            'MarkerSize', jSize,...
            'MarkerEdgeColor', jEdgeColor,...
            'MarkerFaceColor', jFaceColor);
        % Draw third link.
        sx = ex;
        sy = ey;
        ex = sx + l(3)*cos(q(1) + q(2) + q(3));
        ey = sy + l(3)*sin(q(1) + q(2) + q(3));
        line([sx, ex], [sy, ey], 'LineWidth', lWidth);
        % Draw fourth joint.
        plot(axeshandle, ex, ey, 'o',...
            'MarkerSize', jSize,...
            'MarkerEdgeColor', jEdgeColor,...
            'MarkerFaceColor', jFaceColor);
        % Draw fourth link.
        sx = ex;
        sy = ey;
        ex = sx + l(4)*cos(q(1) + q(2) + q(3) + q(4));
        ey = sy + l(4)*sin(q(1) + q(2) + q(3) + q(4));
        line([sx, ex], [sy, ey], 'LineWidth', lWidth);
        % Draw end effector.
        plot(axeshandle, ex, ey, 'o',...
            'MarkerSize', jSize,...
            'MarkerEdgeColor', eeEdgeColor,...
            'MarkerFaceColor', eeFaceColor);
        % Draw obstacles.
        plot(axeshandle, ob(1), ob(2), 'o',...
            'MarkerSize', obSize,...
            'MarkerEdgeColor', obEdgeColor,...
            'MarkerFaceColor', obFaceColor);
        plot(axeshandle, ob(3), ob(4), 'o',...
            'MarkerSize', obSize,...
            'MarkerEdgeColor', obEdgeColor,...
            'MarkerFaceColor', obFaceColor);
    end

    %% Up-button callback function.
    function upbutton_pressed()
       if(ob(4) < 1.5)
           
    end
end