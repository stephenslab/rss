This page describes how to install the
[MATLAB scripts](https://github.com/stephenslab/rss/tree/master/src_vb)
that implement the VB algorithms for RSS.

Step 1. Download the RSS repository: <https://github.com/stephenslab/rss>.

```bash
wget https://github.com/stephenslab/rss/archive/master.zip
unzip master.zip
rm master.zip
mv rss-master rss
```

Step 2. Go to `rss/src_vb` and run `install.m` in MATLAB.

```bash
cd rss/src_vb
matlab -nodesktop
```

The output may look like this:

```matlab
>> run install.m
Building with 'gcc'.
MEX completed successfully.
Compilation of MEX files is complete.
```

After successfully running `install.m`,
you will find a file `rss_varbvsr_update_matlab.mexa64` (in Linux)
in the current directory (`rss/src_vb`), which is the workhorse for VB.

Please note that these scripts have only been tested in Linux.
If you have trouble installing these scripts,
please open an [issue](https://github.com/stephenslab/rss/issues) and
tell us the details of your compiling environment.

### MEX files in MATLAB

The tricky part of this installation is to compile MEX files in MATLAB.

Before running `install.m`, please make sure that
you have the compiler compatible with your version of MATLAB,
and that you can compile the MEX files in the tutorials given on the MathWorks website.

In addition, please ensure that you compile and run these
[VB scripts](https://github.com/stephenslab/rss/tree/master/src_vb)
in the same version of MATLAB.
For example, if you have compiled MEX files using MATLAB Release R2014b,
please use the compiled output `rss_varbvsr_update_matlab.mexa64` only in MATLAB R2014b.

For more information on MEX, please see:

- <http://www.mathworks.com/help/matlab/write-cc-mex-files.html>
- <http://www.mathworks.com/help/matlab/matlab_external/introducing-mex-files.html>
- <http://www.mathworks.com/support/sysreq/previous_releases.html>    
