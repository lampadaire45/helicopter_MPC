function visualize_helicopter_trajectory_rotating(states_trajectory,reference_trajectory,pause_duration)
    %% VISUALIZE QUADROTOR TRAJECTORY
    % Higly inspired by example MPC...
    %
    %
    % plots the dynamics of the quadrotor with inverted pendulum given the 
    % provided states_trajectory which must be as specified below
    % INPUTS
    % - time:               Nx1 vector of time indices
    %       
    % - states_trajectory:  Nx8 vector of quadrotor and pendulum states
    %       1) x        displacement of COM in x-direction
    %       2) y        displacement of COM in y-direction
    %       3) z        displacement of COM in z-direction
    %       4) roll     rotation of quadrotor around x-axis
    %       5) pitch    rotation of quadrotor around y-axis
    %       6) yaw      rotation of quadrotor around z-axis
    % OUTPUS
    % - none 
    

    
    %% INIT
    if (nargin == 2)
        pause_duration = 0;
    end
    
    % X is the 9-states of the helicopter
    x = states_trajectory(:,1);
    y = -states_trajectory(:,2); %negative signs because of different reference frame
    z = -states_trajectory(:,3); %negative signs because of different reference frame
    
    roll = states_trajectory(:,7);
    pitch = -states_trajectory(:,8); %negative signs because of different reference frame
    yaw = -states_trajectory(:,9);   %negative signs because of different reference frame
    
    X = [x y z roll pitch yaw];
    
    [N,~] = size(X);

    % reference set point (center of plot)
    x_r = 0;
    y_r = 0;
    z_r = 0;

    % helicopter frame, main rotor and tail rotor
    r_MR = 1; % Main rotor size
    h_MR = 0.5; % Height of main rotor
    rx_MR = r_MR*cos(linspace(0,2*pi,20));
    rx_MR = [rx_MR rx_MR(1)];
    ry_MR = r_MR*sin(linspace(0,2*pi,20));
    ry_MR = [ry_MR ry_MR(1)];

    r_TR = r_MR*0.5; % Tail rotor size
    rx_TR = r_TR*cos(linspace(0,2*pi,20));
    rx_TR = [rx_TR rx_TR(1)];
    rz_TR = r_TR*sin(linspace(0,2*pi,20));
    rz_TR = [rz_TR rz_TR(1)];
    
    l =  r_MR*2; % Length of fuselage
    pl = 2*r_MR; % Length of vertical arrow

    % init figure
    gca42 = figure(42);
    set(gca42,'color','w');
    clf;
    
    % define plot axes limits
    w = 10;
    wz = 10;
    ic_buf = 1.0;
    
    Ax = [-w+x_r w+x_r+ic_buf -w+y_r w+y_r+ic_buf -wz+z_r wz+z_r];

    % loop through trajectory inputs
    for j = [1:N-1 1]

        % obtain rotational matrix for current RPY angles
        Rt = R(X(j,4:6)); 
        
        % define each quadrotor circle
        R1 = Rt*([ rx_MR ; ry_MR  ; zeros(size(rx_MR))+h_MR ]) + X(j,1:3)';
        R2 = Rt*([ rx_TR-l ; zeros(size(rx_MR)); rz_TR ]) + X(j,1:3)';
        
        % define black arms
        A1 = Rt*([0 -l;0 0;0 0]) + X(j,1:3)';
        
        % verticle reference position of helicopter
        Pv = ([0 0;0 0; 0 pl]) + X(j,1:3)';
        
        % plot coordinate reference
        plot3( x_r,y_r,z_r,'r.');
        hold on
        
        % plot quadrotor frame cross and circles
        plot3( A1(1,:),A1(2,:),A1(3,:),'k','LineWidth',1.5);
        plot3( R1(1,:),R1(2,:),R1(3,:),'r',R2(1,:),R2(2,:),R2(3,:),'b','LineWidth',1.5);
        
        % plot equilibrium equilibrium (verticle upwards)
        quiver3(Pv(1,1),Pv(2,1),Pv(3,1),Pv(1,2)-Pv(1,1),Pv(2,2)-Pv(2,1),Pv(3,2)-Pv(3,1),'g','LineWidth',1.5);
        
        hold off
        grid();
        
        % set axes
        axis(Ax);
        view(3);
        % view([15 25]);
        % view([0 0]);
        
        set(gca,'box','on')
        drawnow   
        if (pause_duration > 0) 
            pause(pause_duration);
        end
    end
    
    % plot the trajectory at the end
    figure(42);
    hold on
    plot3(x,y,z,'.k','MarkerSize',0.2);
    %plot(reference_trajectory(1,:),reference_trajectory(2,:),'c-');
    grid on
    
    % Function to determine the rotation matrix for plotting
    function y = R(Xrot)
        
        phi = Xrot(1);
        theta = Xrot(2);
        psi = Xrot(3);

        Rpsi = [cos(psi)  -sin(psi)   0;
                sin(psi)  cos(psi)    0;
                0         0           1];

        % rotation around y with theta
        Rtheta = [cos(theta)    0       sin(theta);
                  0             1       0;
                 -sin(theta)    0       cos(theta)];

        % rotation around x with phi 
        Rphi = [1       0           0;
                0       cos(phi)    -sin(phi);
                0       sin(phi)    cos(phi)];

        y=Rpsi*Rtheta*Rphi;
    end
    
end