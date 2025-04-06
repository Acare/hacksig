## Resubmission
This is a resubmission. In this version I have:

* Added new signatures to the package;
* Fixed vignette-building bugs related to new versions of suggested packages

## Test environments

- Fedora Linux 41 (local): R 4.4.3;

## R CMD check results

#### Local

0 errors | 0 warnings | 3 notes

1. `Package suggested but not available for checking: ‘covr’`: The package is actually listed in DESCRIPTION;
2. `unable to verify current time`: The *worldclockapi* website is unavailable and this could cause the note;
3. `Packages unavailable to check Rd xrefs: ‘GSVA’, ‘singscore’`: These are Bioconductor packages, they're just referenced in the help pages.
