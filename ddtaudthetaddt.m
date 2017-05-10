clc 
clear all

syms theta1 theta2 theta3 theta1d theta2d theta3d theta1dd theta2dd theta3dd g t

theta1 = sym('theta1(t)');
theta2 = sym('theta2(t)');
theta3 = sym('theta3(t)');
theta1d = sym('theta1d(t)');
theta2d = sym('theta2d(t)');
theta3d = sym('theta3d(t)');
theta1dd = sym('theta1dd(t)');
theta2dd = sym('theta2dd(t)');
theta3dd = sym('theta3dd(t)');

% % no external force
% dtaudtheta1d = 2*(2700*theta1d*sin(theta2 + theta3) - 2700*theta3*sin(theta3) + 9000*theta1d*sin(theta2))*(4278*theta1dd + 4278*theta2dd + 676*theta3dd - 1350*theta3^2*sin(theta3) + 4500*theta1d^2*sin(theta2) - 150*g*cos(theta1 + theta2) + 1350*theta1dd*cos(theta2 + theta3) + 4500*theta1dd*cos(theta2) + 2700*theta1dd*cos(theta3) + 2700*theta2dd*cos(theta3) + 1350*theta3dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) - 2700*theta3*theta1d*sin(theta3) - 2700*theta3*theta2d*sin(theta3)) + 2*(2700*theta3*sin(theta2 + theta3) + 2700*theta2d*sin(theta2 + theta3) + 2700*theta3*sin(theta3) + 9000*theta2d*sin(theta2))*(1350*theta3^2*sin(theta3) - 4278*theta2dd - 676*theta3dd - 11706*theta1dd + 4500*theta2d^2*sin(theta2) + 150*g*cos(theta1 + theta2) - 2700*theta1dd*cos(theta2 + theta3) - 1350*theta2dd*cos(theta2 + theta3) - 1350*theta3dd*cos(theta2 + theta3) + 285*g*cos(theta1) - 9000*theta1dd*cos(theta2) - 2700*theta1dd*cos(theta3) - 4500*theta2dd*cos(theta2) - 2700*theta2dd*cos(theta3) - 1350*theta3dd*cos(theta3) + 45*g*cos(theta1 + theta2 + theta3) + 1350*theta3^2*sin(theta2 + theta3) + 1350*theta2d^2*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta2 + theta3) + 2700*theta3*theta2d*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta3) + 2700*theta3*theta2d*sin(theta3) + 9000*theta1d*theta2d*sin(theta2)) + 2*(2700*theta1d*sin(theta2 + theta3) + 2700*theta1d*sin(theta3) + 2700*theta2d*sin(theta3))*(676*theta1dd + 676*theta2dd + 676*theta3dd + 1350*theta1d^2*sin(theta3) + 1350*theta2d^2*sin(theta3) + 1350*theta1dd*cos(theta2 + theta3) + 1350*theta1dd*cos(theta3) + 1350*theta2dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta3));
%  
%  
% dtaudtheta2d = 2*(2700*theta3*sin(theta2 + theta3) + 2700*theta1d*sin(theta2 + theta3) + 2700*theta2d*sin(theta2 + theta3) + 2700*theta3*sin(theta3) + 9000*theta1d*sin(theta2) + 9000*theta2d*sin(theta2))*(1350*theta3^2*sin(theta3) - 4278*theta2dd - 676*theta3dd - 11706*theta1dd + 4500*theta2d^2*sin(theta2) + 150*g*cos(theta1 + theta2) - 2700*theta1dd*cos(theta2 + theta3) - 1350*theta2dd*cos(theta2 + theta3) - 1350*theta3dd*cos(theta2 + theta3) + 285*g*cos(theta1) - 9000*theta1dd*cos(theta2) - 2700*theta1dd*cos(theta3) - 4500*theta2dd*cos(theta2) - 2700*theta2dd*cos(theta3) - 1350*theta3dd*cos(theta3) + 45*g*cos(theta1 + theta2 + theta3) + 1350*theta3^2*sin(theta2 + theta3) + 1350*theta2d^2*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta2 + theta3) + 2700*theta3*theta2d*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta3) + 2700*theta3*theta2d*sin(theta3) + 9000*theta1d*theta2d*sin(theta2)) + sin(theta3)*(5400*theta1d + 5400*theta2d)*(676*theta1dd + 676*theta2dd + 676*theta3dd + 1350*theta1d^2*sin(theta3) + 1350*theta2d^2*sin(theta3) + 1350*theta1dd*cos(theta2 + theta3) + 1350*theta1dd*cos(theta3) + 1350*theta2dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta3)) - 5400*theta3*sin(theta3)*(4278*theta1dd + 4278*theta2dd + 676*theta3dd - 1350*theta3^2*sin(theta3) + 4500*theta1d^2*sin(theta2) - 150*g*cos(theta1 + theta2) + 1350*theta1dd*cos(theta2 + theta3) + 4500*theta1dd*cos(theta2) + 2700*theta1dd*cos(theta3) + 2700*theta2dd*cos(theta3) + 1350*theta3dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) - 2700*theta3*theta1d*sin(theta3) - 2700*theta3*theta2d*sin(theta3));
%  
%  
% dtaudtheta3d = 0;




% % external force = (0,-10,0)
dtaudtheta1d = 2*(2700*theta1d*sin(theta2 + theta3) + 2700*theta1d*sin(theta3) + 2700*theta2d*sin(theta3))*(676*theta1dd + 676*theta2dd + 676*theta3dd + 1350*theta1d^2*sin(theta3) + 1350*theta2d^2*sin(theta3) + 1350*theta1dd*cos(theta2 + theta3) + 1350*theta1dd*cos(theta3) + 1350*theta2dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta3) - 300) + 2*(2700*theta1d*sin(theta2 + theta3) - 2700*theta3*sin(theta3) + 9000*theta1d*sin(theta2))*(4278*theta1dd + 4278*theta2dd + 676*theta3dd - 300*cos(theta3) - 1350*theta3^2*sin(theta3) + 4500*theta1d^2*sin(theta2) - 150*g*cos(theta1 + theta2) + 1350*theta1dd*cos(theta2 + theta3) + 4500*theta1dd*cos(theta2) + 2700*theta1dd*cos(theta3) + 2700*theta2dd*cos(theta3) + 1350*theta3dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) - 2700*theta3*theta1d*sin(theta3) - 2700*theta3*theta2d*sin(theta3) - 300) + 2*(2700*theta3*sin(theta2 + theta3) + 2700*theta2d*sin(theta2 + theta3) + 2700*theta3*sin(theta3) + 9000*theta2d*sin(theta2))*(300*cos(theta2 + theta3) - 4278*theta2dd - 676*theta3dd - 11706*theta1dd + 300*cos(theta3) + 1350*theta3^2*sin(theta3) + 4500*theta2d^2*sin(theta2) + 150*g*cos(theta1 + theta2) - 2700*theta1dd*cos(theta2 + theta3) - 1350*theta2dd*cos(theta2 + theta3) - 1350*theta3dd*cos(theta2 + theta3) + 285*g*cos(theta1) - 9000*theta1dd*cos(theta2) - 2700*theta1dd*cos(theta3) - 4500*theta2dd*cos(theta2) - 2700*theta2dd*cos(theta3) - 1350*theta3dd*cos(theta3) + 45*g*cos(theta1 + theta2 + theta3) + 1350*theta3^2*sin(theta2 + theta3) + 1350*theta2d^2*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta2 + theta3) + 2700*theta3*theta2d*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta3) + 2700*theta3*theta2d*sin(theta3) + 9000*theta1d*theta2d*sin(theta2) + 300);

dtaudtheta2d = 2*(2700*theta3*sin(theta2 + theta3) + 2700*theta1d*sin(theta2 + theta3) + 2700*theta2d*sin(theta2 + theta3) + 2700*theta3*sin(theta3) + 9000*theta1d*sin(theta2) + 9000*theta2d*sin(theta2))*(300*cos(theta2 + theta3) - 4278*theta2dd - 676*theta3dd - 11706*theta1dd + 300*cos(theta3) + 1350*theta3^2*sin(theta3) + 4500*theta2d^2*sin(theta2) + 150*g*cos(theta1 + theta2) - 2700*theta1dd*cos(theta2 + theta3) - 1350*theta2dd*cos(theta2 + theta3) - 1350*theta3dd*cos(theta2 + theta3) + 285*g*cos(theta1) - 9000*theta1dd*cos(theta2) - 2700*theta1dd*cos(theta3) - 4500*theta2dd*cos(theta2) - 2700*theta2dd*cos(theta3) - 1350*theta3dd*cos(theta3) + 45*g*cos(theta1 + theta2 + theta3) + 1350*theta3^2*sin(theta2 + theta3) + 1350*theta2d^2*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta2 + theta3) + 2700*theta3*theta2d*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta2 + theta3) + 2700*theta3*theta1d*sin(theta3) + 2700*theta3*theta2d*sin(theta3) + 9000*theta1d*theta2d*sin(theta2) + 300) + sin(theta3)*(5400*theta1d + 5400*theta2d)*(676*theta1dd + 676*theta2dd + 676*theta3dd + 1350*theta1d^2*sin(theta3) + 1350*theta2d^2*sin(theta3) + 1350*theta1dd*cos(theta2 + theta3) + 1350*theta1dd*cos(theta3) + 1350*theta2dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) + 2700*theta1d*theta2d*sin(theta3) - 300) - 5400*theta3*sin(theta3)*(4278*theta1dd + 4278*theta2dd + 676*theta3dd - 300*cos(theta3) - 1350*theta3^2*sin(theta3) + 4500*theta1d^2*sin(theta2) - 150*g*cos(theta1 + theta2) + 1350*theta1dd*cos(theta2 + theta3) + 4500*theta1dd*cos(theta2) + 2700*theta1dd*cos(theta3) + 2700*theta2dd*cos(theta3) + 1350*theta3dd*cos(theta3) - 45*g*cos(theta1 + theta2 + theta3) + 1350*theta1d^2*sin(theta2 + theta3) - 2700*theta3*theta1d*sin(theta3) - 2700*theta3*theta2d*sin(theta3) - 300);
 
dtaudtheta3d = 0;
 
ddtaudtheta1ddt = simplify(diff(dtaudtheta1d,t))
ddtaudtheta2ddt = simplify(diff(dtaudtheta2d,t))
ddtaudtheta3ddt = simplify(diff(dtaudtheta3d,t))




 