% This is a small script to compile the necessary MEX files.

% Here is the opportunity to set some additional flags here that may be
% passed to the compiler. These flags tell the GCC compiler to use the ISO
% C99 standard, and to optimize the code as much as possible. Depending
% on the compiler you use to build the MEX shared library files, you may
% want to change these variables, or set them to the empty string ('').
cflags  = '-std=gnu99 -O3 -Os -Dchar16_t=UINT16_T';
ldflags = '-O3 -Os';

% These are the files containing the main functions implemented in C. Note
% that not all these files are needed to compile each of the MEX files.
corefiles = {'C/doublevectormatlab.c '
	     'C/vectorops.c '
	     'C/sigmoid.c '
	     'C/rssvarbvsr.c '
	     'C/sparsematrixmatlab.c '};
	    
% These are the commands to build the build the MEX shared library files.
options = sprintf(['-O -largeArrayDims -IC -I%s ' ...
		   'COPTIMFLAGS="%s" LDOPTIMFLAGS="%s" '],...
		   'C/', cflags,ldflags);

eval(['mex ',options,'C/rss_varbvsr_update_matlab.c ',corefiles{1:5}]);

fprintf('Compilation of MEX files is complete.\n');
