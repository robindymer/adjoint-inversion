function removeSubpaths()
% Removes this project (and its local dependencies) from the MATLAB path.

baseDir = fileparts(mfilename('fullpath'));
rmpath(baseDir);
rmpath(fullfile(baseDir, 'sbplib'));
end