clc, clear all

syms theta1 theta2 theta3 theta1d theta2d theta3d theta1dd theta2dd theta3dd g
% theta1 = sym('theta1(t)');
% theta2 = sym('theta2(t)');
% theta3 = sym('theta3(t)');
% theta1d = sym('theta1d(t)');
% theta2d = sym('theta2d(t)');
% theta3d = sym('theta3d(t)');
% theta1dd = sym('theta1dd(t)');
% theta2dd = sym('theta2dd(t)');
% theta3dd = sym('theta3dd(t)');



q = [theta1 theta2 theta3];
qdot = [theta1d theta2d theta3];
qDdot = [theta1dd theta2dd theta3dd];
% g = 9.8;
[Ti01,Ti12,Ti23,Ti34] = Forward_tau(30, 30, 30, theta1, theta2, theta3);
T(:,:,1) = Ti01;
T(:,:,2) = Ti12;
T(:,:,3) = Ti23;
T(:,:,4) = Ti34;

I(:,:,1) = [0 0 0; 0 0 0; 0 0 3];
I(:,:,2) = [0 0 0; 0 0 0; 0 0 2];
I(:,:,3) = [0 0 0; 0 0 0; 0 0 1];

m = [5 4 3];
Pc = [15 0 0; 15 0 0; 15 0 0];

w(:,1) = qdot(1)*[0; 0; 1];
wdot(:,1) = qDdot(1)*[0; 0; 1];
vdot(:,1) = transpose(T(1:3,1:3,1))*[0; -g; 0];
% vdot(:,1) = [0; 0; 0];
for i = 2:3
    %Angular velocity and acceleration forward propogation:
    w(:,i) = transpose(T(1:3,1:3,i))*w(:,i-1) + qdot(i)*[0; 0; 1];
    wdot(:,i) =  transpose(T(1:3,1:3,i))*wdot(:,i-1) + ...
        cross((transpose(T(1:3,1:3,i))*w(:,i-1)),(qdot(i)*[0; 0; 1])) + ...
        qDdot(i)*[0; 0; 1];
    %Linear acceleration forward propogation:
    vdot(:,i) = transpose(T(1:3,1:3,i))*(cross(wdot(:,i-1),T(1:3,4,i)) + ...
        cross(w(:,i-1),cross(w(:,i-1),T(1:3,4,i))) + vdot(:,i-1));
end
for i = 1:3
    %Acceleration of link mass centers:
    vcdot(:,i) = transpose(cross(wdot(:,i),Pc(i,:)) + ...
        cross(w(:,i),cross(w(:,i),Pc(i,:)))) + vdot(:,i);
    F(:,i) = m(i)*vcdot(:,i);
    N(:,i) = I(:,:,i)*wdot(:,i) + cross(w(:,i),(I(:,:,i)*w(:,i)));
end

% T(:,:,4) = eye(4);
f = sym('A',[3 4]);
% external force on end-effector 10 Newtons
f(:,4) = [0; -10; 0];
n = sym('A',[3 4]);
n(:,4) = [0; 0; 0];

for i = 1:3
    f(:,4-i) = T(1:3,1:3,5-i)*f(:,5-i) + F(:,4-i);
    n(:,4-i) = N(:,4-i) + T(1:3,1:3,5-i)*n(:,5-i) + ...
        transpose(cross(Pc(4-i,:),F(:,4-i))) + ...
        cross(T(1:3,4,5-i),(T(1:3,1:3,5-i)*f(:,5-i)));
    t(4-i) = transpose(n(:,4-i))*[0; 0; 1];
end

t = transpose(t);



tau_sum = simplify(transpose(t)*t);




