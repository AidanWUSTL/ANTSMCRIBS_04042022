clear;

P = pwd;

G = regexp(P, filesep, 'split');

SubjID = G{end};
%SubjID = 'V2_010_2011.08.11';
clear G P;

SubjTempDir = fullfile('SurfReconDeformable', SubjID, 'temp');
SubjMeshesDir = fullfile('SurfReconDeformable', SubjID, 'meshes');

[WhiteV, WhiteF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'white_world.surf'));
[WhiteInternalV, WhiteInternalF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'white+internal_world.surf'));
[WhiteRHV, WhiteRHF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'white-rh_world.surf'));
[WhiteLHV, WhiteLHF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'white-lh_world.surf'));

[PialV, PialF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial_world.surf'));
[PialInternalV, PialInternalF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial+internal_world.surf'));
[PialRHV, PialRHF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial-rh_world.surf'));
[PialLHV, PialLHF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial-lh_world.surf'));

% 
% ILH = find_equal_vertices(WhiteLHV, WhiteInternalV);
% IRH = find_equal_vertices(WhiteRHV, WhiteInternalV);
% 
% IPLH = find_equal_vertices(PialLHV, PialInternalV);
% IPRH = find_equal_vertices(PialRHV, PialInternalV);
% 
% IP65 = find_equal_vertices(Pial6V, Pial5V);
[Pial1V, Pial1F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-1_world.surf'));
[Pial2V, Pial2F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-2_world.surf'));
% [Pial4V, Pial4F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-4_world.surf'));
% [Pial5V, Pial5F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-5_world.surf'));
% [Pial8V, Pial8F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-8_world.surf'));
% Pial8Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'pial-8.vtp'));
% 
% clf;
%^patch('Vertices', Pial8V, 'Faces', Pial8F, 'FaceColor', 'flat', 'FaceVertexCData', Pial8Curv.CollisionType, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
return;

VA = Pial8V(Pial8F(:, 1), :);
VB = Pial8V(Pial8F(:, 2), :);
VC = Pial8V(Pial8F(:, 3), :);
FaceCentroids = (VA + VB + VC) / 3;
R = 5;
%FaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;

BadFaces = Pial8Curv.CollisionType > 0;

VX = FaceCentroids(BadFaces, 1) - FaceCentroids(:, 1)';
VY = FaceCentroids(BadFaces, 2) - FaceCentroids(:, 2)';
VZ = FaceCentroids(BadFaces, 3) - FaceCentroids(:, 3)';

VMAG = sqrt(VX .* VX + VY .* VY + VZ .* VZ);
[~, J] = find(VMAG < 1);
J = unique(J);

patch('Vertices', Pial8V, 'Faces', Pial8F(J, :), 'FaceColor', 'flat', 'FaceVertexCData', Pial8Curv.CollisionType(J), 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);

return;
%%

F = any(Pial5V ~= Pial4V, 2);

BadFaces = any(ismember(Pial5F, find(F)), 2);

AX = zeros(1, 2);
clf;
patch('Vertices', Pial5V, 'Faces', Pial5F(BadFaces, :), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
hold on;
patch('Vertices', Pial4V, 'Faces', Pial4F(BadFaces, :), 'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
axis equal;

%%

%IP54 = find_equal_vertices(Pial5V, Pial4V);
%

% clf;
% 
% patch('Vertices', WhiteLHV, 'Faces', WhiteLHF, 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', double(I), 'EdgeAlpha', 0.4);
% axis equal;

%keyboard;
[PialV, PialF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial_world.surf'));
[PialRHV, PialRHF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial-rh_world.surf'));
[PialLHV, PialLHF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial-lh_world.surf'));
[PialInternalV, PialInternalF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'pial+internal_world.surf'));

WhiteRHCurv = deformable_extract_curvs(fullfile(SubjMeshesDir, 'white-rh.vtp'));
WhiteLHCurv = deformable_extract_curvs(fullfile(SubjMeshesDir, 'white-lh.vtp'));
PialRHCurv = deformable_extract_curvs(fullfile(SubjMeshesDir, 'pial-rh.vtp'));
PialLHCurv = deformable_extract_curvs(fullfile(SubjMeshesDir, 'pial-lh.vtp'));

WhiteCurv = deformable_extract_curvs(fullfile(SubjMeshesDir, 'white.vtp'));
[Pial1V, Pial1F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-1_world.surf'));
[Pial2V, Pial2F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-2_world.surf'));
[Pial3V, Pial3F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-3_world.surf'));
[Pial4V, Pial4F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-4_world.surf'));
[Pial5V, Pial5F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-5_world.surf'));
[Pial6V, Pial6F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-6_world.surf'));

return;

%%
AX = zeros(1, 2);
clf;
AX(1) = subplot(1, 2, 1);
patch('Vertices', WhiteLHV, 'Faces', WhiteLHF(1:100, :), 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', (1:size(WhiteLHV, 1))', 'EdgeAlpha', 0.4);
axis equal;

AX(2) = subplot(1, 2, 2);
patch('Vertices', PialLHV, 'Faces', PialLHF(1:100, :), 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', (1:size(PialLHV, 1))', 'EdgeAlpha', 0.4);
axis equal;
linkprop(AX(:), 'View');

%%
clf;
patch('Vertices', Pial1V, 'Faces', Pial1F, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial1Curv.CortexMask, 'EdgeAlpha', 0.4);
axis equal;

%%

RAS = [-33, 18.38, 21.49];
VA = Pial4V(Pial4F(:, 1), :);
VB = Pial4V(Pial4F(:, 2), :);
VC = Pial4V(Pial4F(:, 3), :);
FaceCentroids = (VA + VB + VC) / 3;
R = 5;
FaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;


clf;
patch('Vertices', Pial4V, 'Faces', Pial4F(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial4Curv.IsOuterFace(FaceI), 'EdgeAlpha', 0.4);
axis equal;

%%
RAS = [13, -19.16, 0.90];
VA = Pial3V(Pial3F(:, 1), :);
VB = Pial3V(Pial3F(:, 2), :);
VC = Pial3V(Pial3F(:, 3), :);
FaceCentroids = (VA + VB + VC) / 3;
R = 5;
FaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;
%FaceI = FaceI & ismember(Pial2Curvs{I}.CollisionType, [5, 6]);
clf;
subplot(1, 2, 1);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType(FaceI), 'EdgeAlpha', 0.4);
patch('Vertices', Pial3V, 'Faces', Pial3F(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial3Curv.IsOuterFace(FaceI), 'EdgeAlpha', 0.4);


VA = Pial4V(Pial4F(:, 1), :);
VB = Pial4V(Pial4F(:, 2), :);
VC = Pial4V(Pial4F(:, 3), :);
FaceCentroids = (VA + VB + VC) / 3;
R = 5;
FaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;
%FaceI = FaceI & ismember(Pial2Curvs{I}.CollisionType, [5, 6]);

subplot(1, 2, 2);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType(FaceI), 'EdgeAlpha', 0.4);
patch('Vertices', Pial4V, 'Faces', Pial3F(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial4Curv.IsOuterFace(FaceI), 'EdgeAlpha', 0.4);

%%
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.Repulsion_magnitude, 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType(FaceI), 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType, 'EdgeAlpha', 0.4);

%[Pial2ZeroV, Pial2ZeroF] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-2-output_000_world.surf'));
%return;
%clf;
%patch('Vertices', Pial1V, 'Faces', Pial1F, 'FaceColor', 'flat', 'FaceVertexCData', Pial1Curv.RegionId);
%#patch('Vertices', Pial1V, 'Faces', Pial1F, 'FaceColor', 'interp', 'FaceVertexCData', Pial1Curv.Status);
%axis equal;
% 
% AX = zeros(1, 2);
% AX(1) = subplot(1, 2, 1);
% 
% patch('Vertices', WhiteV, 'Faces', WhiteF, 'FaceColor', 'flat', 'EdgeAlpha', 0.4, 'FaceVertexCData', WhiteCurv.RegionId);
% axis equal;
% AX(2) = subplot(1, 2, 2);
% 
% patch('Vertices', WhiteV, 'Faces', WhiteF, 'FaceColor', 'flat', 'EdgeAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask);
% axis equal;
% linkprop(AX(:), 'View');

%[Pial5V, Pial5F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-5_world.surf'));

%Pial5Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'pial-5.vtp'));

D = dir(fullfile(SubjTempDir, 'pial-3-output_*.vtp'));
[DN{1:length(D)}] = deal(D.name);

Pial2Curvs = cell(1, length(DN));
Pial2V = cell(1, length(DN));
Pial2F = cell(1, length(DN));

for z = 1:length(DN)
    Pial2Curvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
    T = strrep(DN{z}, '.vtp', '_world.surf');
    [Pial2V{z}, Pial2F{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
end
Pial2Names = DN;

% D = dir(fullfile(SubjTempDir, 'pial-3-output_*.vtp'));
% [DN{1:length(D)}] = deal(D.name);
% 
% Pial3Curvs = cell(1, length(DN));
% Pial3V = cell(1, length(DN));
% Pial3F = cell(1, length(DN));
% 
% for z = 1:length(DN)
%     Pial3Curvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
%     T = strrep(DN{z}, '.vtp', '_world.surf');
%     [Pial3V{z}, Pial3F{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
% end
% 
% Pial2Curvs = [Pial2Curvs, Pial3Curvs];
% Pial2V = [Pial2V, Pial3V];
% Pial3V = [Pial2F, Pial3F];
%Pial2Names = [Pial2Names, DN];


return;
%%
I = 10;
VA = Pial2V{I}(Pial2F{I}(:, 1), :);
VB = Pial2V{I}(Pial2F{I}(:, 2), :);
VC = Pial2V{I}(Pial2F{I}(:, 3), :);
FaceCentroids = (VA + VB + VC) / 3;
R = 5;
FaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;
%FaceI = FaceI & ismember(Pial2Curvs{I}.CollisionType, [5, 6]);
clf;
subplot(1, 1, 1);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType(FaceI), 'EdgeAlpha', 0.4);
patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.Status, 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.Repulsion_magnitude, 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType(FaceI), 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType, 'EdgeAlpha', 0.4);

xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
hold on;
plot3(RAS(1), RAS(2), RAS(3), '*', 'MarkerSize', 50);
view(-25, 14);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}, 'FaceColor', 'flat', 'FaceVertexCData', Pial2Curvs{I}.CollisionType, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.4);
%%
I = 11;
VA = Pial3V{I}(Pial3F{I}(:, 1), :);
VB = Pial3V{I}(Pial3F{I}(:, 2), :);
VC = Pial3V{I}(Pial3F{I}(:, 3), :);
FaceCentroids = (VA + VB + VC) / 3;
R = 15;
FaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;

clf;
subplot(1, 1, 1);
patch('Vertices', Pial3V{I}, 'Faces', Pial3F{I}(FaceI, :), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}(FaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType(FaceI), 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2V{I}, 'Faces', Pial2F{I}, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2Curvs{I}.CollisionType, 'EdgeAlpha', 0.4);

xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
hold on;
plot3(RAS(1), RAS(2), RAS(3), '*', 'MarkerSize', 50);
view(4, 66);

%subplot(1, 2, 2);
%FF = Pial2Curvs{I}.CollisionType ~= 0;
%plot3(FaceCentroids(FF, 1), FaceCentroids(FF, 2), FaceCentroids(FF, 3), '*');
%axis equal;
% in world coordinates +ve is right hemi, -ve is left hemi
% clf;
% patch('Vertices', Pial5V, 'Faces', Pial5F, 'FaceColor', 'flat', 'FaceVertexCData', Pial5Curv.CollisionMask, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1);
% axis equal;
%%
% SubjID = 'MOD2281_V1_a';
% SubjTempDir = fullfile('SurfReconDeformable', SubjID, 'temp');
% 
% [CorpusCallosumMaskNII, CorpusCallosumMaskIMG] = load_nii(fullfile(SubjTempDir, 'corpus-callosum-mask.nii.gz'));
% 
% InvCorpusCallosumMaskQForm = inv(CorpusCallosumMaskNII.hdr.qform);
% 
% [CerebrumRHHullV, CerebrumRHHullF] = freesurfer_read_surf(fullfile(SubjTempDir, 'cerebrum-rh-hull_world.surf'));
% 
% Cerebrum1ImplicitSurfaceFillMask = freesurfer_read_curv(fullfile(SubjTempDir, 'cerebrum-rh-1.ImplicitSurfaceFillMask.curv'));
