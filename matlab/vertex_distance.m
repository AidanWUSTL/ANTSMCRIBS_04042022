function [D] = vertex_distance(V1, V2)

% for each vertex in V1, 
% D(I, J) is the distance from V1(I, :) to V2(J, :)

XC = V1(:, 1) - V2(:, 1)';
YC = V1(:, 2) - V2(:, 2)';
ZC = V1(:, 3) - V2(:, 3)';

D = sqrt(XC .* XC + YC .* YC + ZC .* ZC);