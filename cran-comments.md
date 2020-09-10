## Resubmission
This is a resubmission. In this version I have:

* Made sure all package names, software names and 
  API (application programming interface) names are 
  in single quotes in title and description.
* Added "Gregory Darnell" to the Authors@R field
  as an author and contributor.
* Removed the modification of '.GlobalEnv' in the
  unit tests
* Made sure 'cores' is always set to 1 in vignette
  and all unit tests
* In the comments from the earlier submission it 
  was requested to include a reference for this
  method if possible in DESCRIPTION. Currently a 
  reference is not publicly available online, 
  however we will update the package accordingly 
  once there is. We are aiming to have a preprint 
  up soon, but would like to have the package set 
  up with CRAN before submission to a journal.

## Test environments
* local OS X install, R 4.0.2
* ubuntu 16.04.6 (on travis-ci), R 4.0.0
* win-builder (release & devel), R 4.0.2 & R 2020-08-25 r79074

## R CMD check results
There were no ERRORs or WARNINGs.

There were 1 NOTEs (win-builder):

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michael Turchin <michael_turchin@brown.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  MAPITR (15:33, 19:51)
  MArginal (2:8)
  ePIstasis (2:17)
  epistasis (14:36, 17:54)
  epistatic (22:5)
  iteratively (17:19)
  phenotypes (16:31)

These words are spelled correctly.
