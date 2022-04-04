clear;

P = pwd;

G = regexp(P, filesep, 'split');

SubjID = G{end};
%SubjID = 'MOD1517_V1_a';
clear G P;

SubjTempDir = fullfile('SurfReconDeformable', SubjID, 'temp');
SubjMeshesDir = fullfile('SurfReconDeformable', SubjID, 'meshes');

FullD = dir(fullfile(SubjTempDir, '*.vtp'));
[FullDN{1:length(FullD)}] = deal(FullD.name);
clear FullD;

[WhiteV, WhiteF] = freesurfer_read_surf(fullfile(SubjMeshesDir, 'white_world.surf'));
WhiteCurv = deformable_extract_curvs(fullfile(SubjMeshesDir, 'white.vtp'));

[White2V, White2F] = freesurfer_read_surf(fullfile(SubjTempDir, 'white-2_world.surf'));
White2Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'white-2.vtp'));

[Pial1V, Pial1F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-1_world.surf'));
[Pial2V, Pial2F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-2_world.surf'));
%[Pial3V, Pial3F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-3_world.surf'));
%[Pial4V, Pial4F] = freesurfer_read_surf(fullfile(SubjTempDir, 'pial-4_world.surf'));

Pial1Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'pial-1.vtp'));
Pial2Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'pial-2.vtp'));
%Pial3Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'pial-3.vtp'));
%Pial4Curv = deformable_extract_curvs(fullfile(SubjTempDir, 'pial-4.vtp'));

D = dir(fullfile(SubjTempDir, 'pial-2-output_*_eval.vtp'));
[DN{1:length(D)}] = deal(D.name);

Pial2OutputN = DN;
Pial2OutputCurvs = cell(1, length(DN));
Pial2OutputV = cell(1, length(DN));
Pial2OutputF = cell(1, length(DN));

for z = 1:length(DN)
    Pial2OutputCurvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
    T = strrep(DN{z}, '.vtp', '_world.surf');
    [Pial2OutputV{z}, Pial2OutputF{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
end
clear DN;

D = dir(fullfile(SubjTempDir, 'pial-2-curvature_surface_*_world.vtp'));
[DN{1:length(D)}] = deal(D.name);

%Pial2CurvSurfs = cell(1, length(DN));
Pial2CurvSurfsV = cell(1, length(DN));
Pial2CurvSurfsF = cell(1, length(DN));

for z = 1:length(DN)
    
    %#Pial2OutputCurvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
    T = strrep(DN{z}, '.vtp', '_world.surf');
    [Pial2CurvSurfsV{z}, Pial2CurvSurfsF{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
end
clear DN;
D = dir(fullfile(SubjTempDir, 'pial-2-curvature_surface_*.vtp'));
[DN{1:length(D)}] = deal(D.name);

%Pial2CurvSurfs = cell(1, length(DN));
Pial2CurvSurfsV = cell(1, length(DN));
Pial2CurvSurfsF = cell(1, length(DN));

for z = 1:length(DN)
    
    %#Pial2OutputCurvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
    T = strrep(DN{z}, '.vtp', '_world.surf');
    [Pial2CurvSurfsV{z}, Pial2CurvSurfsF{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
end
clear DN;
D = dir(fullfile(SubjTempDir, 'pial-2-gradient*.vtp'));
[DN{1:length(D)}] = deal(D.name);

%Pial2CurvSurfs = cell(1, length(DN));
Pial2GradientV = cell(1, length(DN));
Pial2GradientF = cell(1, length(DN));

for z = 1:length(DN)
    
    %#Pial2OutputCurvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
    T = strrep(DN{z}, '.vtp', '_world.surf');
    [Pial2GradientV{z}, Pial2GradientF{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
end
clear DN;
D = dir(fullfile(SubjTempDir, 'pial-2-curvature_gradient*.vtp'));
[DN{1:length(D)}] = deal(D.name);

%Pial2CurvSurfs = cell(1, length(DN));
Pial2CurvGradientV = cell(1, length(DN));
Pial2CurvGradientF = cell(1, length(DN));

for z = 1:length(DN)
    
    %#Pial2OutputCurvs{z} = deformable_extract_curvs(fullfile(SubjTempDir, DN{z}));
    T = strrep(DN{z}, '.vtp', '.gradient.surf');
    [Pial2CurvGradientV{z}, Pial2CurvGradientF{z}] = freesurfer_read_surf(fullfile(SubjTempDir, T));
end

FaceI = 201597;

return;
%%
clf;
patch('Vertices', Pial1V, 'Faces', Pial1F, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial1Curv.CortexMask, 'EdgeAlpha', 0.4);
axis equal;

%%
clf;
patch('Vertices', WhiteV, 'Faces', WhiteF, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial1Curv.CortexMask, 'EdgeAlpha', 0.4);
axis equal;

%%
I = 1;

%RAS = [-33, 18.38, 21.49];

%Pial2OutputV{1} = Pial1V

FaceCentroids = face_centroids(White2V, White2F);
WhiteFaceCentroids = face_centroids(WhiteV, WhiteF);

FaceI = 201597;
RAS = FaceCentroids(FaceI, :);
R = 10;
%ValidFaceI = sqrt(sum((VA - RAS) .* (VA - RAS), 2)) < R & sqrt(sum((VB - RAS) .* (VB - RAS), 2)) < R & sqrt(sum((VC - RAS) .* (VC - RAS), 2)) < R;

XC = FaceCentroids - RAS;

ValidFaceI = sqrt(sum(XC .* XC, 2)) < R;
clear XC;

AX = zeros(2, 2);
SR = 2;
SC = 2;

clf;
AX(1) = subplot(SR, SC, 1);

%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2OutputCurvs{I}.CollisionType(ValidFaceI), 'EdgeAlpha', 0.4);
patch('Vertices', WhiteV, 'Faces', WhiteF(ValidFaceI, :), 'FaceColor', 'interp', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.InitialStatus, 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);%patch('Vertices', Pial1V, 'Faces', Pial1F(ValidFaceI, :), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
axis equal;
view(-42, -10);
AX(2) = subplot(SR, SC, 2);

%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2OutputCurvs{I}.CollisionType(ValidFaceI), 'EdgeAlpha', 0.4);
patch('Vertices', WhiteV, 'Faces', WhiteF(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);%patch('Vertices', Pial1V, 'Faces', Pial1F(ValidFaceI, :), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
axis equal;
view(-42, -10);

%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2OutputCurvs{I}.CollisionType(ValidFaceI), 'EdgeAlpha', 0.4);
patch('Vertices', WhiteV, 'Faces', WhiteF(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);%patch('Vertices', Pial1V, 'Faces', Pial1F(ValidFaceI, :), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
axis equal;
view(-42, -10);
AX(3) = subplot(SR, SC, 3);
%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', Pial2OutputCurvs{I}.CollisionType(ValidFaceI), 'EdgeAlpha', 0.4);
patch('Vertices', Pial2V, 'Faces', Pial2F(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);
%patch('Vertices', Pial2OutputV{I}, 'Faces', Pial2OutputF{I}(ValidFaceI, :), 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', WhiteCurv.CortexMask(ValidFaceI), 'EdgeAlpha', 0.4);%patch('Vertices', Pial1V, 'Faces', Pial1F(ValidFaceI, :), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
axis equal;
view(-42, -10);
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
