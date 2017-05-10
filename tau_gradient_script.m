clc
clear all
syms theta1 theta2 theta3 theta1d theta2d theta3d theta1dd theta2dd theta3dd g

%with external force 10N on end-effector
tau_sum = (676*theta1dd + 676*theta2dd + 676*theta3dd + 1350*theta1d^2*sin(theta3) + 1350*theta2d^2*sin(theta3) + 1350*theta1dd*cos(theta2 + theta3) + 1350*theta1dd*cos(theta3) + 1350*theta2dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta3) - 300)^2 + (4278*theta1dd + 4278*theta2dd + 676*theta3dd - 300*cos(theta3) - 1350*theta3^2*sin(theta3) + 4500*theta1d^2*sin(theta2) - 150*g*cos(theta1 + theta2) + 1350*theta1dd*cos(theta2 + theta3) + 4500*theta1dd*cos(theta2) + 2700*theta1dd*cos(theta3) + 2700*theta2dd*cos(theta3) + 1350*theta3dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) - 2700*theta3*theta1d*sin(theta3) - 2700*theta3*theta2d*sin(theta3) - 300)^2 + (300*cos(theta2 + theta3) - 4278*theta2dd - 676*theta3dd - 11706*theta1dd + 300*cos(theta3) + 1350*theta3^2*sin(theta3) + 4500*theta2d^2*sin(theta2) + 150*g*cos(theta1 + theta2) - 2700*theta1dd*cos(theta2 + theta3) - 1350*theta2dd*cos(theta2 + theta3) - 1350*theta3dd*cos(theta2 + theta3) + 285*g*cos(theta1) - 9000*theta1dd*cos(theta2) - 2700*theta1dd*cos(theta3) - 4500*theta2dd*cos(theta2) - 2700*theta2dd*cos(theta3) - 1350*theta3dd*cos(theta3) + 45*g*cos(theta1 + theta2 + theta3) + 1350*theta3^2*sin(theta2 + theta3) + 1350*theta2d^2*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta2 + theta3) + 2700*theta3*theta2d*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta3) + 2700*theta3*theta2d*sin(theta3) + 9000*theta1d*theta2d*sin(theta2) + 300)^2;
  

dtaudtheta1 = simplify(diff(tau_sum,theta1));
dtaudtheta2 = simplify(diff(tau_sum,theta2));
dtaudtheta3 = simplify(diff(tau_sum,theta3));

% dtaudtheta1d = simplify(diff(tau_sum,theta1d));
% dtaudtheta2d = simplify(diff(tau_sum,theta2d));
% dtaudtheta3d = simplify(diff(tau_sum,theta3d));
% 
% dtaudtheta1dd = simplify(diff(tau_sum,theta1dd));
% dtaudtheta2dd = simplify(diff(tau_sum,theta2dd));
% dtaudtheta3dd = simplify(diff(tau_sum,theta3dd));