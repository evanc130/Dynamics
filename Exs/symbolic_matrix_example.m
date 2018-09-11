% symboic_matrix_example.m
%
% last modified 9/6/18 CLee
%
syms theta phi  % define variables as 'symbolic,' can be placed anywhere in code
%
% define matrices 
r1 = [ cos(theta)  -sin(theta);  sin(theta) cos(theta)]
r2 = [ cos(phi)      -sin(phi);  sin(phi)  cos(phi)]
%
r1r2 = r1*r2  % use '*' for matrix mult. use '.*' to mult indvid. elements of matrix
%
% substitue a number in symbolic form for an angle (in symbolic form)
r1sub = subs(r1, theta, 0.4)
% evalutate symbolic values numerically 
n = 6   % number of sig figs for vpa
r1num = vpa(r1sub, n)