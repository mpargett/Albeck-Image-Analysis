%IMAN_ASSERTJAVAPATHS
%   Ensures that Java path definitions are uniformly defined. The
%   iman_pathdefs function must be customized to any new computing or file
%   storage system. 
%

function [pth] = iman_assertjavapaths
%Define desired paths
%--------------------------------------------------------------------------
%   Collect pre-defined paths
p = iman_pathdefs();
p = structfun(@(x)regexprep(x, '([^\\|/])$', '$1\\'), p, 'Un', 0);

%Remove typically problematic MATLAB paths
rmpath(['', p.bad{:}]);

%IF bioformats provided, ensure proper path is used
if isempty(p.bioformats);    pth.jc = [];
else
    pth.jc = {[p.bioformats,'bioformats_package.jar']};
    %NOTE: To be rid of the "log4j:WARN No appenders found for logger" warning,
    %   edit classpath.txt, and add the following path to the end of it
    %   \\mcb.ucdavis.edu\Shared\Data\mcb_labs\albeck\bfmatlab
    
    %Identify sources of the Bioformats java package
    jcp_d = javaclasspath('-dynamic');  jcp_s = javaclasspath('-static');
    %Check for conflicts in static java class path, and report error as needed
    hasbf = ~cellfun(@isempty, regexpi(jcp_s, 'bioformats', 'start'));
    if any(hasbf)   %IF found, check if is the desired package, if not, Error.
        if nnz(hasbf) > 1 || ~strcmpi(jcp_s{hasbf}, pth.jc)
            error(['Conflicting Bioformats Java package found in static ',...
                'Java Class Path.  Edit classpath.txt to remove.']);
        end
    end
    %Now check dynamic path
    hasbf = ~cellfun(@isempty, regexpi(jcp_d, 'bioformats', 'start'));
    if any(hasbf)   %IF found, check if is the desired package
        goodbf = strcmpi(pth.jc, jcp_d(hasbf));
        if any(hasbf & ~goodbf)
            %Undesired Bioformats packages found. Remove from class path
            javarmpath(jcp_d{hasbf & ~goodbf});
        end
    end
    
    %Check for paths with addtitional pre-existing bioformats packages
    pbf_l = which('loci_tools.jar', '-ALL');
    pbf_b = which('bioformats_package.jar', '-ALL');
    %Get paths containing bioformats packages and remove them
    pbf = regexpi([pbf_l; pbf_b], '(.*)\\\w*\.jar', 'tokens');
    if ~isempty(pbf); pbf = [pbf{:}]; pbf = [pbf{:}]; rmpath(pbf{:}); end
end

%Add proper paths
%--------------------------------------------------------------------------
%Code for Image Processing, Bio-Formats, and u-Track
addpath(p.code, p.bioformats, ...
            p.utrack, [p.utrack,'mex'], [p.utrack,'kdtree']);

%Add desired bioformats package to dynamic java class path
if ~isempty(pth.jc) %Only if a new path is available
    jcp_d = javaclasspath('-dynamic');  %Get existing path to append to
    javaclasspath(jcp_d, pth.jc);   %Must include here to have available
end    

end