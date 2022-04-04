function [C] = deformable_extract_curvs(VTPFileName)

if exist(VTPFileName, 'file') ~= 2
    error('VTP not found: ' + VTPFileName)
end

if ~endsWith(VTPFileName, '.vtp')
    error('file must have extension .vtp')
end

[FILEPATH, VTPFilePrefix, ~] = fileparts(VTPFileName);

D = dir(fullfile(FILEPATH, [VTPFilePrefix '.*.curv']));
if isempty(D)
    C = [];
else
    [DN{1:length(D)}] = deal(D.name);

    for z = 1:length(DN)
        
        T = regexp(DN{z}, ['^' VTPFilePrefix '.([^\.]+).curv$'], 'tokens');
        if ~isempty(T)
            CurvName = T{1}{1};
            %disp(CurvName);
            C.(CurvName) = freesurfer_read_curv(fullfile(FILEPATH, DN{z}));
        end
    end
end