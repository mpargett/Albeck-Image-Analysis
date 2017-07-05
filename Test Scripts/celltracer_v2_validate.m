%CELLTRACER_V2_VALIDATE
%   Validation run for celltracer version 2
%   


function [] = celltracer_v2_validate()

%Run test for basic operations, ND2 file
fprintf('\nRunning basic ND2 tests.\n');
try 
    celltracer_v2_test_basic_nd2;  fprintf('Basic test OK.\n\n');
catch me; fprintf('Failure.\n\n'); disp(me);
end

%Run test for frame shift corrections, tiff files (1 per image)
fprintf('\nRunning tiff and frame shift tests.\n');
try 
    celltracer_v2_test_tiff_shift;  fprintf('Tiff and shift test OK.\n\n');
catch me; fprintf('Failure.\n\n'); disp(me);
end

%Run test for multipage tiff files
fprintf('\nRunning MultiPage tiff tests.\n');
try 
    celltracer_v2_test_mptiff;  fprintf('MultiPage tiff test OK.\n\n');
catch me; fprintf('Failure.\n\n'); disp(me);
end

%Run test for multiline light source (SPECTRAX), ND2 file
fprintf('\nRunning SPECTRAX light source test.\n');
try 
    celltracer_v2_test_spectrax;  fprintf('SPECTRAX test OK.\n\n');
catch me; fprintf('Failure.\n\n'); disp(me);
end

%Run test for spectral unmixing, ND2 file
fprintf('\nRunning spectral unmixing test.\n');
try 
    celltracer_v2_test_unmix;  fprintf('Unmixing test OK.\n\n');
catch me; fprintf('Failure.\n\n'); disp(me);
end

%Remove output files
try   delete('celltracer_v2_test_*.mat'); end %#ok<TRYNC>


%% Run image analysis procedure
% [dao, GMD, dmx] = iman_celltracer(ip, op);



end

