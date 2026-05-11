function addSubpaths()
% Adds this project (and its local dependencies) to the MATLAB path.

baseDir = fileparts(mfilename('fullpath'));
addpath(baseDir);
addpath(fullfile(baseDir, 'sbplib'));
end