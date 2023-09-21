This external folder contains some binary programs from Advanced Normalization Tools (ANTs).

There are only a few pre-compiled programs for MacOS and Linux.
We donwloaded a package of programs for each platform on website: https://github.com/stnava/ANTs/releases
Then, we extracted only the required programs, since the entire package size is large.

It was created a function in iftUtil.h called iftSetANTsEnvironment which configures the environment to run ANTs programs inside libIFT.

If you want more ANTs programs, download the entire package of programs for MacOS (Yosemite) and Linux, and extracted the desired packages into subfolders MacOS and Linux.
