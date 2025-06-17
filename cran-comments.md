## R CMD check results

0 errors | 0 warnings | 0 notes


This is a patch to resolve a problem notified to me by the maintainer of ggplot2: my code was making an invalid call to ggplot2 which
would raise an error in the forthcoming version.

I fixed the problem and have tested against the forthcoming version. Have tested on the package build with R CMD check on windows, 
ubuntu and macos - latest, ubuntu-latest devel.
