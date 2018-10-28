This page shows how to install the Monte Carlo Markov chain (MCMC) scripts
in [rss/src/](https://github.com/stephenslab/rss/tree/master/src).

Step 1. Download the RSS repository: <https://github.com/stephenslab/rss>.
```
wget https://github.com/stephenslab/rss/archive/master.zip
unzip master.zip
rm master.zip
mv rss-master rss
```
Step 2. Download and install the [lightspeed](http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/) MATLAB toolbox (author: Tom Minka).
```
cd rss/src
wget http://ftp.research.microsoft.com/downloads/db1653f0-1308-4b45-b358-d8e1011385a0/lightspeed.zip
unzip lightspeed.zip
rm lightspeed.zip
cd lightspeed/
matlab -nodisplay < install_lightspeed.m
```

Step 3. Download and install the [lapack](http://www.mathworks.com/matlabcentral/fileexchange/16777-lapack) MATLAB package (author: Tim Toolan).
```
cd rss/src
unzip lapack.zip
rm lapackhelp.m lapack.zip license.txt
```
Notice that if an appropriate compiled version of `lapack.c` does not exist, the package will ask whether to build one.  

Step 4. Download the [progress](http://www.mathworks.com/matlabcentral/fileexchange/8564-progress) MATLAB package (author: Martinho Marta-Almeida).
```
unzip progress.zip
rm progress.zip license.txt 
```
Files downloaded from Step 2-4 must be placed under `rss/src`.
```
-rw-rw-r-- 1 xiangzhu xiangzhu   5261 2015-11-11 08:46 calc_posterior_bvsr.m
-rw-rw-r-- 1 xiangzhu xiangzhu    972 2015-11-11 08:46 compute_pve.m
-rw-rw-r-- 1 xiangzhu xiangzhu  86772 2015-11-11 12:01 lapack.c
-rw-rw-r-- 1 xiangzhu xiangzhu   5927 2015-11-11 12:01 lapack.m
-rwxrwxr-x 1 xiangzhu xiangzhu 120123 2015-11-11 12:12 lapack.mexa64
drwx------ 6 xiangzhu xiangzhu  32768 2015-11-11 11:54 lightspeed
-rw-rw-r-- 1 xiangzhu xiangzhu   2573 2015-11-11 12:01 progress.m
-rw-rw-r-- 1 xiangzhu xiangzhu   4402 2015-11-11 08:46 propose_gamma.m
-rw-rw-r-- 1 xiangzhu xiangzhu  11629 2015-11-11 08:46 rss_ash.m
-rw-rw-r-- 1 xiangzhu xiangzhu  22182 2015-11-11 08:46 rss_bslmm.m
-rw-rw-r-- 1 xiangzhu xiangzhu   8363 2015-11-11 08:46 rss_bvsr.m
-rw-rw-r-- 1 xiangzhu xiangzhu  11348 2015-11-11 08:46 update_betatilde.m
-rw-rw-r-- 1 xiangzhu xiangzhu   4796 2015-11-11 08:46 update_bz.m
-rw-rw-r-- 1 xiangzhu xiangzhu  17205 2015-11-11 08:46 update_zlabel.m
```
