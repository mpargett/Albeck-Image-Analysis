%IMAN_PATHDEFS
%   Define relevant paths for this installation.
%
%   EDIT this function to reflect the installation paths on your computer
%   or network.

%   This version for: Albeck Lab, with network storage in MCB, UC DAVIS

function p = iman_pathdefs()
%Declare proper path for image processing code
% p.code = '';  %Default sets no new path
p.code = ['\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis\;',...
    '\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis\uTrack Mods\'];
%     ['\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis\Test Scripts\'];
        
%Declare proper path for Bio-Formats package
% p.bioformats = '';    %Default sets no new path, but must be set unless
%                       bioformats.jar is in the Working Directory when
%                       processing is run
p.bioformats = '\\albecklab.mcb.ucdavis.edu\data\bfmatlab\';

%Declare proper path for uTrack
% p.utrack = '';    %Default sets no new path
p.utrack = '\\albecklab.mcb.ucdavis.edu\data\u-track_2.1.0\software\';
        
%Declare common problematic paths
%   If desired can list undesirable paths (warning will show when they are
%   not present), or can search the current path for them via regular
%   expressions
pth = path;
p.bad = regexpi(pth, ...
    '[^;]*\\albecklab.mcb.ucdavis.edu\\data\\Code[^;]*;', 'match');
p.bad = [p.bad,regexpi(pth, ...
    '[^;]*\\albecklab.mcb.ucdavis.edu\\data\\u-track[^;]*;', 'match')];
end

