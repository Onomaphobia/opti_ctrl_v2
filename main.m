clc
clear variables
close all

OPTIMAL_CONTROL = true;
OPTIMAL_CONTROL_trajectory = true;

% Links length[cm]
L1 = 30;
L2 = 30;
L3 = 30;

% initial angles[degree]
a1 = 60;
a2 = -60;
a3 = -60;

a1 = deg2rad(a1);
a2 = deg2rad(a2);
a3 = deg2rad(a3);

% a1i=0;
% a2i=0;
% a3i=0; 

% disp('Please enter the goal position and orientation of the end-effector in xyz coordinates')
xg = input('enter the x coordinate:');
yg = input('enter the y coordinate:');
% input('Please enter the desired linear velocity of the gripper in cm/s:')
v  = 10;

% initial position transformation matrices use forward kinmatics function
% the Forward function is to calculate the transformation matrix using the
% initial angles of revolute joints and initial length of prismatic joint

disp('The initial transformation matrix is ')
[Ti01,Ti02,Ti03,Ti04] = Forward(L1,L2,L3,a1,a2,a3);
Ti04
% Initial position vector of the end-effector
disp('the initial position of the end-effector')
xi = Ti04(1,4)
yi = Ti04(2,4)

% No IK, n is predefined
% n = 50;



% Calculating the linear distance from the initial position to the desired position and the linear velocity:
D=sqrt((xg-Ti04(1,4))^2 + (yg-Ti04(2,4))^2); %Distance from initial position to goal position
dt=0.05; % Time increment in seconds.
time=D/v;  % Total time of animation.
n=round(time/dt); % Number of iterations rounded up.
dt=time/n; % Adjusted time increment in seconds. 

% the constant commanded Cartesian rates
Xdot = (xg-Ti04(1,4))/n; %[cm/sec]
Ydot = (yg-Ti04(2,4))/n; %[cm/sec]
dx = [Xdot;Ydot];


%% %data initialization
% t = zeros(n+1,1);
% X = zeros(n+1,1);
% y = zeros(n+1,1);
Theta = zeros(3,n+1);
Theta_pinv = zeros(3,n+1);
Thetadot = zeros(3,n+1);
lambda = zeros(3,n+1);
U = zeros(3,n);
G = zeros(3,n);
G_pre = zeros(3,n);
Scale = zeros(3,n);
Scale_pre = zeros(3,n);
DFDQT = zeros(3,3);
DLDQ = zeros(3,1);

%qdot_noX = (I-J#J)U
Thetadot_noX = zeros(3,n+1);
cost_opti = zeros(n+1,1);
cost_goal_position = zeros(n+1,1);
count_max = 100;
sum_opti_cost = zeros(count_max,1);
sum_cost = zeros(count_max,1);
alpha = 0.1;
ALPHA = 0.002;
beta = 1;

count = 1;
t(1) = 0;
Theta_ini(1,1) = a1;
Theta_ini(2,1) = a2;    
Theta_ini(3,1) = a3;
Theta(1,1) = a1;
Theta(2,1) = a2;    
Theta(3,1) = a3;   
Theta_pinv(1,1) = a1;    
Theta_pinv(2,1) = a2;
Theta_pinv(3,1) = a3; 

cost_goal_position(1,1) = G_sum(L1,L2,L3,Theta(1,1),Theta(2,1),Theta(3,1),9.8);

% generate trajectory from pinvJ
i = 2;
while i < n + 2
    %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
    %dq = pinvJ*dx;
    
    Thetadot(:,i) = Qdot(L1,L2,L3,Theta_pinv(1,i-1),Theta_pinv(2,i-1),Theta_pinv(3,i-1),Xdot,Ydot,0,0,0);
    Theta_pinv(1,i) = Theta_pinv(1,i-1) + Thetadot(1,i);
    Theta_pinv(2,i) = Theta_pinv(2,i-1) + Thetadot(2,i);
    Theta_pinv(3,i) = Theta_pinv(3,i-1) + Thetadot(3,i);
%     cost_pinv(i,1) = cost(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i));

    i = i + 1;
end

i = 2; 
%NO IK, end-effector stay
while i < n + 2
    cost_goal_position(i,1) = G_sum(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),9.8);
    i = i + 1;
end

%% % optimal control
if OPTIMAL_CONTROL_trajectory == true
	%timer start
	tic
    
    %initialize theta_opti
    Theta_opti = Theta_pinv;
    
    %finite differentiation for dldq
    Theta_opti_dot = diff(Theta_opti,1,2);
    Theta_opti_dotdot = diff(Theta_opti_dot,1,2);
    
    Theta_opti_dot(:,size(Theta_opti,2)) = Theta_opti_dot(:,size(Theta_opti,2)-1);
    
    Theta_opti_dotdot(:,size(Theta_opti,2)) = Theta_opti_dotdot(:,size(Theta_opti,2)-2);
    Theta_opti_dotdot(:,size(Theta_opti,2)-1) = Theta_opti_dotdot(:,size(Theta_opti,2)-2);
    
	%calculation of lambda
	i = n;
	while i > 0
	    DFDQT = dfdqT(L1,L2,L3,Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),Xdot,Ydot,U(1,i),U(2,i),U(3,i));
	    DLDQ = dldq(Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),Theta_opti_dot(1,i),Theta_opti_dot(2,i),Theta_opti_dot(3,i),Theta_opti_dotdot(1,i),Theta_opti_dotdot(2,i),Theta_opti_dotdot(3,i),9.8);
	    lambda(:,i) = DLDQ + DFDQT * lambda(:,i+1) + lambda(:,i+1);
	    i = i - 1;
	end

	%calculation of G
	i = 1;
	while i < n + 1
	    DFDUT = dfduT(L1,L2,L3,Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),U(1,i),U(2,i),U(3,i));
	    DLDU = dldu(); % DLDU = 0
	    G(:,i) = DLDU + DFDUT * lambda(:,i+1);
	    i = i + 1;
	end

	Scale = -G;
    scale_max = max(max(abs(Scale),[],2));
    alpha = ALPHA/scale_max;

	U = U + alpha*Scale;
	%store G and Scale
	G_pre = G;
	Scale_pre = Scale;

	%calculation of Theta for the second iteration
	i = 2;
	while i < n + 2
	    %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
	    %dq = pinvJ*dx;
	    Thetadot(:,i) = Qdot(L1,L2,L3,Theta_opti(1,i-1),Theta_opti(2,i-1),Theta_opti(3,i-1),Xdot,Ydot,U(1,i-1),U(2,i-1),U(3,i-1));
	    
	    Theta_opti(1,i) = Theta_opti(1,i-1) + Thetadot(1,i);
	    Theta_opti(2,i) = Theta_opti(2,i-1) + Thetadot(2,i);
	    Theta_opti(3,i) = Theta_opti(3,i-1) + Thetadot(3,i);
	    
	    cost_opti(i,1) = tau_sum(Theta_opti(1,i-1),Theta_opti(2,i-1),Theta_opti(3,i-1),Theta_opti_dot(1,i-1),Theta_opti_dot(2,i-1),Theta_opti_dot(3,i-1),Theta_opti_dotdot(1,i-1),Theta_opti_dotdot(2,i-1),Theta_opti_dotdot(3,i-1),9.8);
	    
	    i = i + 1;
	end

	%count of CGM iteration
	
	sum_opti_cost(count,1) = sum(cost_opti);


    while count < count_max 
	    
	    %calculation of lambda
	    i = n;
	    while i > 0
	        DFDQT = dfdqT(L1,L2,L3,Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),Xdot,Ydot,U(1,i),U(2,i),U(3,i));
            DLDQ = dldq(Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),Theta_opti_dot(1,i),Theta_opti_dot(2,i),Theta_opti_dot(3,i),Theta_opti_dotdot(1,i),Theta_opti_dotdot(2,i),Theta_opti_dotdot(3,i),9.8);
	        lambda(:,i) = DLDQ + DFDQT * lambda(:,i+1) + lambda(:,i+1);
	        i = i - 1;
	    end

	    %calculation of G
	    i = 1;
	    while i < n + 1
	        DFDUT = dfduT(L1,L2,L3,Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),U(1,i),U(2,i),U(3,i));
            DLDU = dldu(); % DLDU = 0
            G(:,i) = DLDU + DFDUT * lambda(:,i+1);
            i = i + 1;
	    end
	    
	    beta = sum(sum(G.^2))/(sum(sum(G_pre.^2)));
	    Scale = -G + beta*Scale_pre;
%         Scale = -G;
        scale_max = max(max(abs(Scale),[],2));
        alpha = ALPHA/scale_max;
	    U = U + alpha*Scale;
	    %store G and Scale
	    G_pre = G;
	    Scale_pre = Scale;
	    
	    %count of CGM iteration
	    count = count + 1;
	    
	    %calculation of Theta for next iteration
	    i = 2;
        while i < n + 2

            %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
            %dq = pinvJ*dx;
            Thetadot(:,i) = Qdot(L1,L2,L3,Theta_opti(1,i-1),Theta_opti(2,i-1),Theta_opti(3,i-1),Xdot,Ydot,U(1,i-1),U(2,i-1),U(3,i-1));

            Theta_opti(1,i) = Theta_opti(1,i-1) + Thetadot(1,i);
            Theta_opti(2,i) = Theta_opti(2,i-1) + Thetadot(2,i);
            Theta_opti(3,i) = Theta_opti(3,i-1) + Thetadot(3,i);

            cost_opti(i,1) = tau_sum(Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i),Theta_opti_dot(1,i),Theta_opti_dot(2,i),Theta_opti_dot(3,i),Theta_opti_dotdot(1,i),Theta_opti_dotdot(2,i),Theta_opti_dotdot(3,i),9.8);
            i = i + 1;
        end

	    sum_opti_cost(count,1) = sum(cost_opti);
% 	    if (sum_cost(count-1,1) > sum_cost(count,1) && (sum_cost(count-1,1) - sum_cost(count,1))/sum_cost(count-1,1) < 8e-5) || (count > 20 &&sum_cost(count,1) > sum_cost(count-1,1))
% 	        break
% 	    end
    end

	count;
	toc
end

 
%% %data initialization
% t = zeros(n+1,1);
% X = zeros(n+1,1);
% y = zeros(n+1,1);
Theta = zeros(3,n+1);
Theta_pinv = zeros(3,n+1);
Thetadot = zeros(3,n+1);
lambda = zeros(3,n+1);
U = zeros(3,n);
G = zeros(3,n);
G_pre = zeros(3,n);
Scale = zeros(3,n);
Scale_pre = zeros(3,n);
DFDQT = zeros(3,3);
DLDQ = zeros(3,1);

%qdot_noX = (I-J#J)U
Thetadot_noX = zeros(3,n+1);
cost_goal_position = zeros(n+1,1);
count_max = 100;
sum_opti_cost = zeros(count_max,1);
sum_cost = zeros(count_max,1);
alpha = 0.1;
ALPHA = 0.002;
beta = 1;

count = 1;
t(1) = 0;
Theta_ini(1,1) = a1;
Theta_ini(2,1) = a2;    
Theta_ini(3,1) = a3;
Theta(1,1) = a1;
Theta(2,1) = a2;    
Theta(3,1) = a3;   
Theta_pinv(1,1) = a1;    
Theta_pinv(2,1) = a2;
Theta_pinv(3,1) = a3; 

cost_goal_position(1,1) = G_sum(L1,L2,L3,Theta(1,1),Theta(2,1),Theta(3,1),9.8);

% generate trajectory from pinvJ
i = 2;
while i < n + 2
    %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
    %dq = pinvJ*dx;
    
    Thetadot(:,i) = Qdot(L1,L2,L3,Theta_pinv(1,i-1),Theta_pinv(2,i-1),Theta_pinv(3,i-1),Xdot,Ydot,0,0,0);
    Theta_pinv(1,i) = Theta_pinv(1,i-1) + Thetadot(1,i);
    Theta_pinv(2,i) = Theta_pinv(2,i-1) + Thetadot(2,i);
    Theta_pinv(3,i) = Theta_pinv(3,i-1) + Thetadot(3,i);
%     cost_pinv(i,1) = cost(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i));

    i = i + 1;
end

%copy goal configuration for optimal control of goal configuration
Theta(1,:) = Theta_pinv(1,n+1);
Theta(2,:) = Theta_pinv(2,n+1);
Theta(3,:) = Theta_pinv(3,n+1);

i = 2; 
%NO IK, end-effector stay
while i < n + 2
    cost_goal_position(i,1) = G_sum(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),9.8);
    i = i + 1;
end


%% %optimal_control CGM + CHOMP
if OPTIMAL_CONTROL == true
	%timer start
	tic

	%calculation of lambda
	i = n;
	while i > 0
	    DFDQT = dQdot_noXdqT(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),U(1,i),U(2,i),U(3,i));
	    DLDQ = dG_compdq(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),9.8);
	    lambda(:,i) = DLDQ + DFDQT * lambda(:,i+1) + lambda(:,i+1);
	    i = i - 1;
	end

	%calculation of G
	i = 1;
	while i < n + 1
	    DFDUT = dQdot_noXduT(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),U(1,i),U(2,i),U(3,i));
	    DLDU = dG_compdu(); % DLDU = 0
	    G(:,i) = DLDU + DFDUT * lambda(:,i+1);
	    i = i + 1;
	end

	Scale = -G;
    scale_max = max(max(abs(Scale),[],2));
    alpha = ALPHA/scale_max;

	U = U + alpha*Scale;
	%store G and Scale
	G_pre = G;
	Scale_pre = Scale;

	%calculation of Theta for the second iteration
	i = 2;
	while i < n + 2
	    %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
	    %dq = pinvJ*dx;
	    Thetadot_noX(:,i) = Qdot_noX(L1,L2,L3,Theta(1,i-1),Theta(2,i-1),Theta(3,i-1),U(1,i-1),U(2,i-1),U(3,i-1));
	    
	    Theta(1,i) = Theta(1,i-1) + Thetadot_noX(1,i);
	    Theta(2,i) = Theta(2,i-1) + Thetadot_noX(2,i);
	    Theta(3,i) = Theta(3,i-1) + Thetadot_noX(3,i);
	    
	    cost_goal_position(i,1) = G_sum(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),9.8);
	    
	    i = i + 1;
	end

	%count of CGM iteration
	
	sum_cost(count,1) = sum(cost_goal_position);


    while count < count_max 
	    
	    %calculation of lambda
	    i = n;
	    while i > 0
	        DFDQT = dQdot_noXdqT(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),U(1,i),U(2,i),U(3,i));
	        DLDQ = dG_compdq(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),9.8);
	        lambda(:,i) = DLDQ + DFDQT * lambda(:,i+1) + lambda(:,i+1);
	        i = i - 1;
	    end

	    %calculation of G
	    i = 1;
	    while i < n + 1
	        DFDUT = dQdot_noXduT(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),U(1,i),U(2,i),U(3,i));
            DLDU = dG_compdu(); % DLDU = 0
            G(:,i) = DLDU + DFDUT * lambda(:,i+1);
            i = i + 1;
	    end
	    
	    beta = sum(sum(G.^2))/(sum(sum(G_pre.^2)));
	    Scale = -G + beta*Scale_pre;
%         Scale = -G;
        scale_max = max(max(abs(Scale),[],2));
        alpha = ALPHA/scale_max;
	    U = U + alpha*Scale;
	    %store G and Scale
	    G_pre = G;
	    Scale_pre = Scale;
	    
	    %count of CGM iteration
	    count = count + 1;
	    
	    %calculation of Theta for next iteration
	    i = 2;
        while i < n + 2
        
            %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
            %dq = pinvJ*dx;
            Thetadot_noX(:,i) = Qdot_noX(L1,L2,L3,Theta(1,i-1),Theta(2,i-1),Theta(3,i-1),U(1,i-1),U(2,i-1),U(3,i-1));

            Theta(1,i) = Theta(1,i-1) + Thetadot_noX(1,i);
            Theta(2,i) = Theta(2,i-1) + Thetadot_noX(2,i);
            Theta(3,i) = Theta(3,i-1) + Thetadot_noX(3,i);

            cost_goal_position(i,1) = G_sum(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i),9.8);
            i = i + 1;
        end

	    sum_cost(count,1) = sum(cost_goal_position);
% 	    if (sum_cost(count-1,1) > sum_cost(count,1) && (sum_cost(count-1,1) - sum_cost(count,1))/sum_cost(count-1,1) < 8e-5) || (count > 20 &&sum_cost(count,1) > sum_cost(count-1,1))
% 	        break
% 	    end
    end

	count;
	toc
    
    Theta_fin(1,1) = Theta(1,n+1);
    Theta_fin(2,1) = Theta(2,n+1);
    Theta_fin(3,1) = Theta(3,n+1);
    cost_goal_pinv = G_sum(L1,L2,L3,Theta_pinv(1,n+1),Theta_pinv(2,n+1),Theta_pinv(3,n+1),9.8);
    cost_goal_opti = G_sum(L1,L2,L3,Theta_fin(1),Theta_fin(2),Theta_fin(3),9.8);

    Theta(1,:) = linspace(Theta_ini(1),Theta_fin(1),n+1);
    Theta(2,:) = linspace(Theta_ini(2),Theta_fin(2),n+1);
    Theta(3,:) = linspace(Theta_ini(3),Theta_fin(3),n+1);
end




%CHOMP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tic
if OPTIMAL_CONTROL == true
    Theta_CHOMP = Theta;
else
    Theta_CHOMP = Theta_pinv;
end
chomp_iterator = 1;
CHOMP_iteration = 200;
tau_cost = zeros(CHOMP_iteration,1);

%finite_differencing matrix
Diff_rules(1, 1) = 0;
Diff_rules(1, 2) = -1 / 12.0;
Diff_rules(1, 3) = 16 / 12.0;
Diff_rules(1, 4) = -30 / 12.0;
Diff_rules(1, 5) = 16 / 12.0;
Diff_rules(1, 6) = -1 / 12.0;
Diff_rules(1, 7) = 0;
dimension = size(Theta_CHOMP,2);
i = 0;
finite_diff = zeros(dimension);
while i <=dimension
    i = i + 1;
    j = -3;
    while j < 3
        j = j + 1;
        index = i + j;
        if index <= 0
            continue
        end
        if index >dimension
            continue
        end
        finite_diff(i,index) = Diff_rules(1,j + 4);
    end
end
A = finite_diff'*finite_diff;
% A = eye(dimension);
%CHOMP interation
while chomp_iterator < CHOMP_iteration+1
    Thetadot = diff(Theta_CHOMP,1,2);
    Thetadotdot = diff(Thetadot,1,2);
    Thetadotdotdot = diff(Thetadotdot,1,2);
    Thetadotdotdotdot = diff(Thetadotdotdot,1,2);
    
    %fix dimension
    Thetadot(:,size(Theta_CHOMP,2)) = Thetadot(:,size(Theta_CHOMP,2)-1);
    
    Thetadotdot(:,size(Theta_CHOMP,2)) = Thetadotdot(:,size(Theta_CHOMP,2)-2);
    Thetadotdot(:,size(Theta_CHOMP,2)-1) = Thetadotdot(:,size(Theta_CHOMP,2)-2);
    
    Thetadotdotdot(:,size(Theta_CHOMP,2)) = Thetadotdotdot(:,size(Theta_CHOMP,2)-3);
    Thetadotdotdot(:,size(Theta_CHOMP,2)-1) = Thetadotdotdot(:,size(Theta_CHOMP,2)-3);
    Thetadotdotdot(:,size(Theta_CHOMP,2)-2) = Thetadotdotdot(:,size(Theta_CHOMP,2)-3);
    
    Thetadotdotdotdot(:,size(Theta_CHOMP,2)) = Thetadotdotdotdot(:,size(Theta_CHOMP,2)-4);
    Thetadotdotdotdot(:,size(Theta_CHOMP,2)-1) = Thetadotdotdotdot(:,size(Theta_CHOMP,2)-4);
    Thetadotdotdotdot(:,size(Theta_CHOMP,2)-2) = Thetadotdotdotdot(:,size(Theta_CHOMP,2)-4);
    Thetadotdotdotdot(:,size(Theta_CHOMP,2)-3) = Thetadotdotdotdot(:,size(Theta_CHOMP,2)-4);
    
    
    traj_iterator = 1;
    tau_traj = zeros(size(Theta_CHOMP,2),1);
    tau_grad = zeros(3,size(Theta_CHOMP,2));
    while traj_iterator < size(Theta_CHOMP,2) + 1
        tau_traj(traj_iterator,1) = tau_sum(Theta_CHOMP(1,traj_iterator),Theta_CHOMP(2,traj_iterator),Theta_CHOMP(3,traj_iterator),Thetadot(1,traj_iterator),Thetadot(2,traj_iterator),Thetadot(3,traj_iterator),Thetadotdot(1,traj_iterator),Thetadotdot(2,traj_iterator),Thetadotdot(3,traj_iterator),9.8);
        tau_grad(:,traj_iterator) = tau_gradient(Theta_CHOMP(1,traj_iterator),Theta_CHOMP(2,traj_iterator),Theta_CHOMP(3,traj_iterator),Thetadot(1,traj_iterator),Thetadot(2,traj_iterator),Thetadot(3,traj_iterator),Thetadotdot(1,traj_iterator),Thetadotdot(2,traj_iterator),Thetadotdot(3,traj_iterator),Thetadotdotdot(1,traj_iterator),Thetadotdotdot(2,traj_iterator),Thetadotdotdot(3,traj_iterator),Thetadotdotdotdot(1,traj_iterator),Thetadotdotdotdot(2,traj_iterator),Thetadotdotdotdot(3,traj_iterator),9.8);
    traj_iterator = traj_iterator + 1;
    end
    tau_cost(chomp_iterator,1) = sum(tau_traj);
    if chomp_iterator == 1
        tau_traj_ini = tau_traj;
    end
    if chomp_iterator == CHOMP_iteration - 1
        tau_traj_finish = tau_traj;
    end
     
    
    for joint_count = 1:3
        final_increment(joint_count,:) = transpose(A \ transpose(tau_grad(joint_count,:)));
        increment_max(joint_count) = max(final_increment(joint_count,:));
        increment_min(joint_count) = min(final_increment(joint_count,:));
        scale(joint_count) =0.005 / max(abs(increment_max(joint_count)), abs(increment_min(joint_count)));
        if scale(joint_count) > 1
            scale(joint_count) = 1;
        end
%         max_scale(joint_count) = 0.1/abs(increment_max(joint_count));
%         min_scale(joint_count) = 0.1/abs(increment_min(joint_count));
%         if max_scale(joint_count) < 1
%             scale(joint_count) = max_scale(joint_count);
%         end
%         if min_scale(joint_count) < 1
%             scale(joint_count) = min_scale(joint_count);
%         end
        final_increment(joint_count,:) = scale(joint_count) * final_increment(joint_count,:);
    end
    
    %update theta
    final_increment(1,:);
    Theta_CHOMP = Theta_CHOMP - final_increment;
    chomp_iterator = chomp_iterator + 1;
%      pause
end
toc

%% chomp without optimal control
%CHOMP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tic

Theta_CHOMP_NO_opti = Theta_pinv;

chomp_iterator = 1;
CHOMP_iteration = 200;
tau_cost = zeros(CHOMP_iteration,1);

%finite_differencing matrix
Diff_rules(1, 1) = 0;
Diff_rules(1, 2) = -1 / 12.0;
Diff_rules(1, 3) = 16 / 12.0;
Diff_rules(1, 4) = -30 / 12.0;
Diff_rules(1, 5) = 16 / 12.0;
Diff_rules(1, 6) = -1 / 12.0;
Diff_rules(1, 7) = 0;
dimension = size(Theta_CHOMP_NO_opti,2);
i = 0;
finite_diff = zeros(dimension);
while i <=dimension
    i = i + 1;
    j = -3;
    while j < 3
        j = j + 1;
        index = i + j;
        if index <= 0
            continue
        end
        if index >dimension
            continue
        end
        finite_diff(i,index) = Diff_rules(1,j + 4);
    end
end
A = finite_diff'*finite_diff;
% A = eye(dimension);
%CHOMP interation
while chomp_iterator < CHOMP_iteration+1
    Thetadot = diff(Theta_CHOMP_NO_opti,1,2);
    Thetadotdot = diff(Thetadot,1,2);
    Thetadotdotdot = diff(Thetadotdot,1,2);
    Thetadotdotdotdot = diff(Thetadotdotdot,1,2);
    
    %fix dimension
    Thetadot(:,size(Theta_CHOMP_NO_opti,2)) = Thetadot(:,size(Theta_CHOMP_NO_opti,2)-1);
    
    Thetadotdot(:,size(Theta_CHOMP_NO_opti,2)) = Thetadotdot(:,size(Theta_CHOMP_NO_opti,2)-2);
    Thetadotdot(:,size(Theta_CHOMP_NO_opti,2)-1) = Thetadotdot(:,size(Theta_CHOMP_NO_opti,2)-2);
    
    Thetadotdotdot(:,size(Theta_CHOMP_NO_opti,2)) = Thetadotdotdot(:,size(Theta_CHOMP_NO_opti,2)-3);
    Thetadotdotdot(:,size(Theta_CHOMP_NO_opti,2)-1) = Thetadotdotdot(:,size(Theta_CHOMP_NO_opti,2)-3);
    Thetadotdotdot(:,size(Theta_CHOMP_NO_opti,2)-2) = Thetadotdotdot(:,size(Theta_CHOMP_NO_opti,2)-3);
    
    Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)) = Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-4);
    Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-1) = Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-4);
    Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-2) = Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-4);
    Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-3) = Thetadotdotdotdot(:,size(Theta_CHOMP_NO_opti,2)-4);
    
    
    traj_iterator = 1;
    tau_traj_NO_opti = zeros(size(Theta_CHOMP_NO_opti,2),1);
    tau_grad = zeros(3,size(Theta_CHOMP_NO_opti,2));
    while traj_iterator < size(Theta_CHOMP_NO_opti,2) + 1
        tau_traj_NO_opti(traj_iterator,1) = tau_sum(Theta_CHOMP_NO_opti(1,traj_iterator),Theta_CHOMP_NO_opti(2,traj_iterator),Theta_CHOMP_NO_opti(3,traj_iterator),Thetadot(1,traj_iterator),Thetadot(2,traj_iterator),Thetadot(3,traj_iterator),Thetadotdot(1,traj_iterator),Thetadotdot(2,traj_iterator),Thetadotdot(3,traj_iterator),9.8);
        tau_grad(:,traj_iterator) = tau_gradient(Theta_CHOMP_NO_opti(1,traj_iterator),Theta_CHOMP_NO_opti(2,traj_iterator),Theta_CHOMP_NO_opti(3,traj_iterator),Thetadot(1,traj_iterator),Thetadot(2,traj_iterator),Thetadot(3,traj_iterator),Thetadotdot(1,traj_iterator),Thetadotdot(2,traj_iterator),Thetadotdot(3,traj_iterator),Thetadotdotdot(1,traj_iterator),Thetadotdotdot(2,traj_iterator),Thetadotdotdot(3,traj_iterator),Thetadotdotdotdot(1,traj_iterator),Thetadotdotdotdot(2,traj_iterator),Thetadotdotdotdot(3,traj_iterator),9.8);
    traj_iterator = traj_iterator + 1;
    end
    tau_cost(chomp_iterator,1) = sum(tau_traj_NO_opti);
    if chomp_iterator == 1
        tau_traj_ini = tau_traj_NO_opti;
    end
    if chomp_iterator == CHOMP_iteration - 1
        tau_traj_finish = tau_traj_NO_opti;
    end
     
    
    for joint_count = 1:3
        final_increment(joint_count,:) = transpose(A \ transpose(tau_grad(joint_count,:)));
        increment_max(joint_count) = max(final_increment(joint_count,:));
        increment_min(joint_count) = min(final_increment(joint_count,:));
        scale(joint_count) =0.005 / max(abs(increment_max(joint_count)), abs(increment_min(joint_count)));
        if scale(joint_count) > 1
            scale(joint_count) = 1;
        end
%         max_scale(joint_count) = 0.1/abs(increment_max(joint_count));
%         min_scale(joint_count) = 0.1/abs(increment_min(joint_count));
%         if max_scale(joint_count) < 1
%             scale(joint_count) = max_scale(joint_count);
%         end
%         if min_scale(joint_count) < 1
%             scale(joint_count) = min_scale(joint_count);
%         end
        final_increment(joint_count,:) = scale(joint_count) * final_increment(joint_count,:);
    end
    
    %update theta
    final_increment(1,:);
    Theta_CHOMP_NO_opti = Theta_CHOMP_NO_opti - final_increment;
    chomp_iterator = chomp_iterator + 1;
%      pause
end
toc


%% figure

line = zeros(3,n+1);
line(:,1) = Ti04(1:3,4);

figure('Name','Arm animation','NumberTitle','off')
% Plots of the Arm
hold on

% link0=plot3([0,Ti01(1,4)],[0,Ti01(2,4)],[0,Ti01(3,4)],'color',1-0.25*(1-[1,0,0]),'LineWidth',5);
% link1=plot3([Ti01(1,4), Ti02(1,4)],[Ti01(2,4), Ti02(2,4)],[Ti01(3,4), Ti02(3,4)],'color',1-0.25*(1-[0,0,1]),'LineWidth',5);
% link2=plot3([Ti02(1,4), Ti03(1,4)],[Ti02(2,4), Ti03(2,4)],[Ti02(3,4), Ti03(3,4)],'color',1-0.25*(1-[0,1,0]),'LineWidth',5);
% link3=plot3([Ti03(1,4), Ti04(1,4)],[Ti03(2,4), Ti04(2,4)],[Ti03(3,4), Ti04(3,4)],'color',1-0.25*(1-[1,0,0]),'LineWidth',5);
% joint1=plot3(Ti01(1,4),Ti01(2,4),Ti01(3,4),'-ko','MarkerEdgeColor',1-0.25*(1-[0,0,0]),'MarkerFaceColor',1-0.25*(1-[.49 1 .63]),'MarkerSize',10);
% joint2=plot3(Ti02(1,4),Ti02(2,4),Ti02(3,4),'-ko','MarkerEdgeColor',1-0.25*(1-[0,0,0]),'MarkerFaceColor',1-0.25*(1-[.49 1 .63]),'MarkerSize',10);
% joint3=plot3(Ti03(1,4),Ti03(2,4),Ti03(3,4),'-ko','MarkerEdgeColor',1-0.25*(1-[0,0,0]),'MarkerFaceColor',1-0.25*(1-[.49 1 .63]),'MarkerSize',10);
% endeffector=plot3(Ti04(1,4),Ti04(2,4),Ti04(3,4),'-o','color',1-0.25*(1-[0,0,0]),'LineWidth',5,'MarkerSize',5);


% link0_opti=plot3([0,Ti01(1,4)],[0,Ti01(2,4)],[0,Ti01(3,4)],'color',1-0.5*(1-[1,0,0]),'LineWidth',5);
% link1_opti=plot3([Ti01(1,4), Ti02(1,4)],[Ti01(2,4), Ti02(2,4)],[Ti01(3,4), Ti02(3,4)],'color',1-0.5*(1-[0,0,1]),'LineWidth',5);
% link2_opti=plot3([Ti02(1,4), Ti03(1,4)],[Ti02(2,4), Ti03(2,4)],[Ti02(3,4), Ti03(3,4)],'color',1-0.5*(1-[0,1,0]),'LineWidth',5);
% link3_opti=plot3([Ti03(1,4), Ti04(1,4)],[Ti03(2,4), Ti04(2,4)],[Ti03(3,4), Ti04(3,4)],'color',1-0.5*(1-[1,0,0]),'LineWidth',5);
% joint1_opti=plot3(Ti01(1,4),Ti01(2,4),Ti01(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
% joint2_opti=plot3(Ti02(1,4),Ti02(2,4),Ti02(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
% joint3_opti=plot3(Ti03(1,4),Ti03(2,4),Ti03(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
% endeffector_opti=plot3(Ti04(1,4),Ti04(2,4),Ti04(3,4),'-o','color',1-0.5*(1-[0,0,0]),'LineWidth',5,'MarkerSize',5);

% link0_CHOMP_NO_opti=plot3([0,Ti01(1,4)],[0,Ti01(2,4)],[0,Ti01(3,4)],'color',1-0.5*(1-[1,0,0]),'LineWidth',5);
% link1_CHOMP_NO_opti=plot3([Ti01(1,4), Ti02(1,4)],[Ti01(2,4), Ti02(2,4)],[Ti01(3,4), Ti02(3,4)],'color',1-0.75*(1-[0,1,0]),'LineWidth',5);
% link2_CHOMP_NO_opti=plot3([Ti02(1,4), Ti03(1,4)],[Ti02(2,4), Ti03(2,4)],[Ti02(3,4), Ti03(3,4)],'color',1-0.75*(1-[0,0,1]),'LineWidth',5);
% link3_CHOMP_NO_opti=plot3([Ti03(1,4), Ti04(1,4)],[Ti03(2,4), Ti04(2,4)],[Ti03(3,4), Ti04(3,4)],'color',1-0.75*(1-[1,0,0]),'LineWidth',5);
% joint1_CHOMP_NO_opti=plot3(Ti01(1,4),Ti01(2,4),Ti01(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
% joint2_CHOMP_NO_opti=plot3(Ti02(1,4),Ti02(2,4),Ti02(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
% joint3_CHOMP_NO_opti=plot3(Ti03(1,4),Ti03(2,4),Ti03(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
% endeffector_CHOMP_NO_opti=plot3(Ti04(1,4),Ti04(2,4),Ti04(3,4),'-o','color',1-0.5*(1-[0,0,0]),'LineWidth',5,'MarkerSize',5);

link0_CHOMP=plot3([0,Ti01(1,4)],[0,Ti01(2,4)],[0,Ti01(3,4)],'color',1-0.5*(1-[1,0,0]),'LineWidth',5);
link1_CHOMP=plot3([Ti01(1,4), Ti02(1,4)],[Ti01(2,4), Ti02(2,4)],[Ti01(3,4), Ti02(3,4)],'color',1-1*(1-[0,1,0]),'LineWidth',5);
link2_CHOMP=plot3([Ti02(1,4), Ti03(1,4)],[Ti02(2,4), Ti03(2,4)],[Ti02(3,4), Ti03(3,4)],'color',1-1*(1-[0,0,1]),'LineWidth',5);
link3_CHOMP=plot3([Ti03(1,4), Ti04(1,4)],[Ti03(2,4), Ti04(2,4)],[Ti03(3,4), Ti04(3,4)],'color',1-1*(1-[1,0,0]),'LineWidth',5);
joint1_CHOMP=plot3(Ti01(1,4),Ti01(2,4),Ti01(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
joint2_CHOMP=plot3(Ti02(1,4),Ti02(2,4),Ti02(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
joint3_CHOMP=plot3(Ti03(1,4),Ti03(2,4),Ti03(3,4),'-ko','MarkerEdgeColor',1-0.5*(1-[0,0,0]),'MarkerFaceColor',1-0.5*(1-[.49 1 .63]),'MarkerSize',10);
endeffector_CHOMP=plot3(Ti04(1,4),Ti04(2,4),Ti04(3,4),'-o','color',1-0.5*(1-[0,0,0]),'LineWidth',5,'MarkerSize',5);

% x=-300:100:300; % simulation plane 
% y=x;
% [x,y]=meshgrid(x,y);
% z=x*0;
% surf(x,y,z);
% % shading interp
% alpha(0.6)
grid on;
axis([-80 80 -80 80 -80 80 -80 80])
view(2);

% Plots of points of interest on the system:
initial=plot3(Ti04(1,4),Ti04(2,4),Ti04(3,4),'-ro','LineWidth',5,'MarkerSize',5);
goal=plot3(xg,yg,0,'-go','LineWidth',5,'MarkerSize',5);
title('ARM Animation'); xlabel('x, (cm)'); ylabel('y (cm)'); zlabel('z (cm)');


% pause


i = 2;
while i < n + 2
    %[pinvJ]= pinvJ_2D(L1,L2,L3,a1,a2,a3);
    %dq = pinvJ*dx;
%     Thetadot(:,i) = Qdot(L1,L2,L3,Theta(1,i-1),Theta(2,i-1),Theta(3,i-1),Xdot,Ydot,U(1,i-1),U(2,i-1),U(3,i-1));
% 
%     Theta(1,i) = Theta(1,i-1) + Thetadot(1,i);
%     Theta(2,i) = Theta(2,i-1) + Thetadot(2,i);
%     Theta(3,i) = Theta(3,i-1) + Thetadot(3,i);

%     cost_CGM(i,1) = cost(L1,L2,L3,Theta(1,i),Theta(2,i),Theta(3,i

%     [T01p,T02p,T03p,T04p] = Forward(L1,L2,L3,Theta_pinv(1,i),Theta_pinv(2,i),Theta_pinv(3,i));
%     set(link1,'XData',[T01p(1,4), T02p(1,4)],'YData',[T01p(2,4), T02p(2,4)],'ZData',[T01p(3,4), T02p(3,4)]);
%     set(link2,'XData',[T02p(1,4), T03p(1,4)],'YData',[T02p(2,4), T03p(2,4)],'ZData',[T02p(3,4), T03p(3,4)]);
%     set(link3,'XData',[T03p(1,4), T04p(1,4)],'YData',[T03p(2,4), T04p(2,4)],'ZData',[T03p(3,4), T04p(3,4)]);
%     set(joint1,'XData',T01p(1,4),'YData',T01p(2,4),'ZData',T01p(3,4));    
%     set(joint2,'XData',T02p(1,4),'YData',T02p(2,4),'ZData',T02p(3,4));
%     set(joint3,'XData',T03p(1,4),'YData',T03p(2,4),'ZData',T03p(3,4));
%     set(endeffector,'XData',T04p(1,4),'YData',T04p(2,4),'ZData',T04p(3,4)); 
%     line(:,i) = T04p(1:3,4);
%     plot3(line(1,1:i),line(2,1:i),line(3,1:i),'r','LineWidth',2);
    
%     [T01_opti,T02_opti,T03_opti,T04_opti] = Forward(L1,L2,L3,Theta_opti(1,i),Theta_opti(2,i),Theta_opti(3,i));
%     set(link1_opti,'XData',[T01_opti(1,4), T02_opti(1,4)],'YData',[T01_opti(2,4), T02_opti(2,4)],'ZData',[T01_opti(3,4), T02_opti(3,4)]);
%     set(link2_opti,'XData',[T02_opti(1,4), T03_opti(1,4)],'YData',[T02_opti(2,4), T03_opti(2,4)],'ZData',[T02_opti(3,4), T03_opti(3,4)]);
%     set(link3_opti,'XData',[T03_opti(1,4), T04_opti(1,4)],'YData',[T03_opti(2,4), T04_opti(2,4)],'ZData',[T03_opti(3,4), T04_opti(3,4)]);
%     set(joint1_opti,'XData',T01_opti(1,4),'YData',T01_opti(2,4),'ZData',T01_opti(3,4));    
%     set(joint2_opti,'XData',T02_opti(1,4),'YData',T02_opti(2,4),'ZData',T02_opti(3,4));
%     set(joint3_opti,'XData',T03_opti(1,4),'YData',T03_opti(2,4),'ZData',T03_opti(3,4));
%     set(endeffector_opti,'XData',T04_opti(1,4),'YData',T04_opti(2,4),'ZData',T04_opti(3,4)); 
%     line(:,i) = T04_opti(1:3,4);
%     plot3(line(1,1:i),line(2,1:i),line(3,1:i),'r','LineWidth',2);
    
%     [T01_CHOMP_NO_opti,T02_CHOMP_NO_opti,T03_CHOMP_NO_opti,T04_CHOMP_NO_opti] = Forward(L1,L2,L3,Theta_CHOMP_NO_opti(1,i),Theta_CHOMP_NO_opti(2,i),Theta_CHOMP_NO_opti(3,i));
%     set(link1_CHOMP_NO_opti,'XData',[T01_CHOMP_NO_opti(1,4), T02_CHOMP_NO_opti(1,4)],'YData',[T01_CHOMP_NO_opti(2,4), T02_CHOMP_NO_opti(2,4)],'ZData',[T01_CHOMP_NO_opti(3,4), T02_CHOMP_NO_opti(3,4)]);
%     set(link2_CHOMP_NO_opti,'XData',[T02_CHOMP_NO_opti(1,4), T03_CHOMP_NO_opti(1,4)],'YData',[T02_CHOMP_NO_opti(2,4), T03_CHOMP_NO_opti(2,4)],'ZData',[T02_CHOMP_NO_opti(3,4), T03_CHOMP_NO_opti(3,4)]);
%     set(link3_CHOMP_NO_opti,'XData',[T03_CHOMP_NO_opti(1,4), T04_CHOMP_NO_opti(1,4)],'YData',[T03_CHOMP_NO_opti(2,4), T04_CHOMP_NO_opti(2,4)],'ZData',[T03_CHOMP_NO_opti(3,4), T04_CHOMP_NO_opti(3,4)]);
%     set(joint1_CHOMP_NO_opti,'XData',T01_CHOMP_NO_opti(1,4),'YData',T01_CHOMP_NO_opti(2,4),'ZData',T01_CHOMP_NO_opti(3,4));    
%     set(joint2_CHOMP_NO_opti,'XData',T02_CHOMP_NO_opti(1,4),'YData',T02_CHOMP_NO_opti(2,4),'ZData',T02_CHOMP_NO_opti(3,4));
%     set(joint3_CHOMP_NO_opti,'XData',T03_CHOMP_NO_opti(1,4),'YData',T03_CHOMP_NO_opti(2,4),'ZData',T03_CHOMP_NO_opti(3,4));
%     set(endeffector_CHOMP_NO_opti,'XData',T04_CHOMP_NO_opti(1,4),'YData',T04_CHOMP_NO_opti(2,4),'ZData',T04_CHOMP_NO_opti(3,4)); 
%     line(:,i) = T04_CHOMP_NO_opti(1:3,4);
%     plot3(line(1,1:i),line(2,1:i),line(3,1:i),'r','LineWidth',2);
    
    [T01_CHOMP,T02_CHOMP,T03_CHOMP,T04_CHOMP] = Forward(L1,L2,L3,Theta_CHOMP(1,i),Theta_CHOMP(2,i),Theta_CHOMP(3,i));
    set(link1_CHOMP,'XData',[T01_CHOMP(1,4), T02_CHOMP(1,4)],'YData',[T01_CHOMP(2,4), T02_CHOMP(2,4)],'ZData',[T01_CHOMP(3,4), T02_CHOMP(3,4)]);
    set(link2_CHOMP,'XData',[T02_CHOMP(1,4), T03_CHOMP(1,4)],'YData',[T02_CHOMP(2,4), T03_CHOMP(2,4)],'ZData',[T02_CHOMP(3,4), T03_CHOMP(3,4)]);
    set(link3_CHOMP,'XData',[T03_CHOMP(1,4), T04_CHOMP(1,4)],'YData',[T03_CHOMP(2,4), T04_CHOMP(2,4)],'ZData',[T03_CHOMP(3,4), T04_CHOMP(3,4)]);
    set(joint1_CHOMP,'XData',T01_CHOMP(1,4),'YData',T01_CHOMP(2,4),'ZData',T01_CHOMP(3,4));    
    set(joint2_CHOMP,'XData',T02_CHOMP(1,4),'YData',T02_CHOMP(2,4),'ZData',T02_CHOMP(3,4));
    set(joint3_CHOMP,'XData',T03_CHOMP(1,4),'YData',T03_CHOMP(2,4),'ZData',T03_CHOMP(3,4));
    set(endeffector_CHOMP,'XData',T04_CHOMP(1,4),'YData',T04_CHOMP(2,4),'ZData',T04_CHOMP(3,4)); 
    line(:,i) = T04_CHOMP(1:3,4);
    plot3(line(1,1:i),line(2,1:i),line(3,1:i),'r','LineWidth',2);
    
    %screenshot
%     if i == n+1 || mod(i,50) == 0
%         pause;
%     end
        
    pause(0.05)
    
    i = i + 1;
    
end

Theta_pinvdot = diff(Theta_pinv,1,2);
Theta_pinvdotdot = diff(Theta_pinv,1,2);

Theta_pinvdot(:,size(Theta_pinv,2)) = Theta_pinvdot(:,size(Theta_pinv,2)-1);  

Theta_pinvdotdot(:,size(Theta_pinv,2)) = Theta_pinvdotdot(:,size(Theta_pinv,2)-2);
Theta_pinvdotdot(:,size(Theta_pinv,2)-1) = Theta_pinvdotdot(:,size(Theta_pinv,2)-2);

traj_iterator = 1;
    tau_traj_pinv = zeros(size(Theta_pinv,2),1);
    while traj_iterator < size(Theta_pinv,2) + 1
        tau_traj_pinv(traj_iterator,1) = tau_sum(Theta_pinv(1,traj_iterator),Theta_pinv(2,traj_iterator),Theta_pinv(3,traj_iterator),Theta_pinvdot(1,traj_iterator),Theta_pinvdot(2,traj_iterator),Theta_pinvdot(3,traj_iterator),Theta_pinvdotdot(1,traj_iterator),Theta_pinvdotdot(2,traj_iterator),Theta_pinvdotdot(3,traj_iterator),9.8);
    traj_iterator = traj_iterator + 1;
    end



hold off;

figure('Name','cost value vs iteration','NumberTitle','off');
plot(sum_cost(1:count,1),'LineWidth',2);

figure('Name','tau vs tau_pinv','NumberTitle','off');
plot(tau_traj,'DisplayName','chomp with opti','LineWidth',2);
hold on;
plot(tau_traj_pinv,'DisplayName','pinv','LineWidth',2);
plot(tau_traj_NO_opti,'DisplayName','chomp NO opti','LineWidth',2);
plot(cost_opti,'DisplayName','opti','LineWidth',2);
hold off;