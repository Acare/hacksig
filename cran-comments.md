## Resubmission
This is a resubmission. In this version I have:

* Fixed a bug in the `hack_sig()` function related to the `"singscore"` method.

## Test environments

- Ubuntu 20.04.3 LTS (local): R 4.1,2;
- Win-builder: R-devel, R-release and R-oldrelease;
- R-hub builder:
  * Fedora Linux, R-devel, clang, gfortran;
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC;
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit;
- R-CMD-check GitHub Actions:
  * macOS-latest (R-release);
  * Windows-latest (R-release);
  * Ubuntu-latest (R-devel);
  * Ubuntu-latest (R-release);
  * Ubuntu-latest (R-oldrelease);

## R CMD check results

0 errors | 0 warnings | 1 note

**the NOTE appears only when using R-hub builder on Windows Server 2022, R-devel, 64 bit**

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
