% This is a small script to compile the necessary MEX files.
% This file is modified from the following file in varbvs package.
% https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/install.m

% These are the files containing the main functions implemented in C. Note
% that not all these files are needed to compile each of the MEX files.
corefiles = {'C/doublevectormatlab.c '
	     'C/vectorops.c '
	     'C/sigmoid.c '
	     'C/rssvarbvsr.c '
	     'C/sparsematrixmatlab.c '};
	    
% These are the commands to build the build the MEX shared library files.
options = '-O -largeArrayDims -IC COPTIMFLAGS="-std=gnu99"';

eval(['mex ',options,' C/rss_varbvsr_update_matlab.c ',corefiles{1:5}]);

fprintf('Compilation of MEX files is complete.\n');
