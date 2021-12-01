## Resubmission
This is a resubmission. In this version I have:

* Replaced `\dontrun{}` with `\donttest{}` because examples in 'hack_sig.R' take more than 5 seconds to run;

* Added references for the three available score computation methods in the Description section of the DESCRIPTION file. Also, I have added a sentence for using the `get_sig_info()` function to retrieve information (including references) about the 19 gene signatures implemented so far. Listing all the 19 references would have unnecessarily extended the Description section.

## Test environments

- Ubuntu 20.04.3 LTS (local): R 4.0.5;
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

There was a message about possibly misspelled words in DESCRIPTION (**Foroutan** at 24:24, **GSEA** at 23:24, **al** at 22:38, 23:43, 24:36, **et** at 22:35, 23:40, 24:33, **singscore** at 24:10, **transcriptomics** at 19:37).

- All these words are actually spelled correctly: "single sample GSEA" and "singscore" are two computation methods, whereas "Foroutan" is the first author of the singscore publication;
- "et al." indicates other authors for listed references in the Description section of DESCRIPTION.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
