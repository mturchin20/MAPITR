# Version News and Updates

Version updates will be tracked and explained here. Major updates & releases will be particularly highlighted.

## MAPITR v1.1.2

###### Summary
* Fixed pathway naming issue in results output

###### Since previous version (v1.1.1)
* Better hard-coded results output so pathway names stay intact
* Also switched 'Results' output to be a data.frame and not a matrix 

###### Notes

###### Next steps (if applicable)


## MAPITR v1.1.1

###### Summary
* Further tweaks for no-OpenMP version to match CRAN specifications

###### Since previous version (v1.1.0)
* Included proper 'if/else' statement in .cpp to get no-OpenMP versions working

###### Notes
* Followed some examples from https://github.com/stephenslab/mashr to fix this on recommendation from Peter
* Used R-hub for additional test environment checks

###### Next steps (if applicable)


## MAPITR v1.1.0

###### Summary
* Included no-OpenMP versions

###### Since previous version (v1.0.6)
* Included no-OpenMP versions of the main 'MAPITR' code
* Updated Roxygen comments, vignette, and unit tests accordingly

###### Notes

###### Next steps (if applicable)


## MAPITR v1.0.6

###### Summary
* Edits in prep for manuscript submission

###### Since previous version (v1.0.5)
* Edits to DESCRIPTION, README, and vignette in prep for manuscript submission

###### Notes

###### Next steps (if applicable)


## MAPITR v1.0.5

###### Summary
* Edits for CRAN resubmission

###### Since previous version (v1.0.4)
* Removed use of '.GlobalEnv' in unit tests as requested
* Made sure 'cores' is always set to 1 in vignettes and unit tests as requested
* Made a formatting edit to 'Description' field in DESCRIPTION as requested

###### Notes

###### Next steps (if applicable)


## MAPITR v1.0.4

###### Summary
* Tweaks for CRAN resubmission

###### Since previous version (v1.0.3)
* Another tweak to the CRAN-specific MIT license format

###### Notes

###### Next steps (if applicable)


## MAPITR v1.0.3

###### Summary
* Tweaks for CRAN resubmission

###### Since previous version (v1.0.2)
* Update MIT license format for CRAN template version

###### Notes

###### Next steps (if applicable)


## MAPITR v1.0.2

###### Summary
* Tweaks for CRAN resubmission

###### Since previous version (v1.0.1)
* Made leaner version of simulated dataset for unit tests and vignette

###### Notes
* Updated unit tests and vignette to match new, leaner simulated dataset

###### Next steps (if applicable)


## MAPITR v1.0.1

###### Summary
* Reopened functionality for covariate usage and included example in current vignette

###### Since previous version (v1.0.0)
* Added the `Covariate` option back into the `MAPITR` function
* Included example of `Covariate` option usage into vignette

###### Notes

###### Next steps (if applicable)


## MAPITR v1.0.0

###### Summary
* First CRAN-ready release and submission

###### Since previous version (v0.1.5)
* Cleaned up files as needed for CRAN submission
* Cleaned up README.md

###### Notes

###### Next steps (if applicable)


## MAPITR v0.1.5

###### Summary
* Added unit tests, Roxygen2 comments, a vignette, and a pkgdown Github page
* Got Travis-CI page going

###### Since previous version (v0.1.0)
* Pausing within-code data check lines for now (only covering obvious cases, and most important check -- the genotype matrix -- too time intensive as a specific, internal step for the main function)
* Finished unit tests
* Finished vignette
* Completed pkgdown + Github page setup
* Added repo to Travis-CI page
* Added 'wCovs' version

###### Notes
* Added Roxygen2 comments for .rda files in '/data'
* Worked through 'devtools::check()' notes, warnings, and errors
* Peter Carbonetto provided some edits that account for v0.1.2, v0.1.3, & v0.1.4 (all edits for getting last parts of Travis-CI setup working)

###### Next steps (if applicable)
* Add in 'with covariate' version and possible 'no OpenMP' versions


## MAPITR v0.1.0

###### Summary
* Minimal working example now running

###### Since previous version (v#.#.#)

###### Notes

###### Next steps (if applicable)
* Finish off rest of wrapper functions
* Implement 'with covariate' and 'without OpenMP' options
* Build up unit tests
* Build up vignette


<!---
## MAPITR v#.#.#

###### Summary

###### Since previous version (v#.#.#)

###### Notes

###### Next steps (if applicable)
-->

