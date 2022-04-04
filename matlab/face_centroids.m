function [C] = face_centroids(V, F)

%^C = mean(cat(3, V(F(:, 1), :), V(F(:, 2), :), V(F(:, 3), :)), 3); 

% the method below is much faster
VA = V(F(:, 1), :);
VB = V(F(:, 2), :);
VC = V(F(:, 3), :);

C = (VA + VB + VC) / 3;

