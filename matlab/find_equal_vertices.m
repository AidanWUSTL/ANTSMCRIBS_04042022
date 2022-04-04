function [IDX] = find_equal_vertices(V1, V2)

VV = [V1; V2];

[VVSorted, VVSortedI] = sortrows(VV);

%V1IInSorted = VVSortedI(1:size(V1, 1));
%[~, V1IInSorted] = ismember((1:size(V1, 1))', VVSortedI);

V1IInSorted = VVSortedI;
V1IInSorted(VVSortedI) = 1:length(VVSortedI);
V1IInSorted = V1IInSorted(1:size(V1, 1));

%keyboard;
%I = all(VVSorted(V1IInSorted, :) == VVSorted(mod(V1IInSorted, size(V1, 1)) + 1, :), 2);
I = all(VVSorted(V1IInSorted, :) == VVSorted(mod(V1IInSorted, size(VVSorted, 1)) + 1, :), 2);
IDX = repmat(-1, size(V1, 1), 1);
IDX(I) = VVSortedI(V1IInSorted(I) + 1) - size(V1, 1);

%keyboard;