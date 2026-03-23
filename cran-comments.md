## R CMD check results

0 errors | 0 warnings | 0 notes


This is a patch to resolve a problem notified to me by CRAN maintainers: changes in Bioconductor packages have caused a warning about dependencies. The fix is to update the minimum required versions of R and Bioconductor in the installation instructions, and to update the NEWS file to reflect this change.

I fixed the problem and have tested against the forthcoming version. Have tested on the package build with R CMD check on windows, 
ubuntu and macos - latest, ubuntu-latest devel.
