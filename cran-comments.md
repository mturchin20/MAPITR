## Test environments
* local OS X install, R 4.0.2
* ubuntu 16.04.6 (on travis-ci), R 4.0.0
* win-builder (release & devel), R 4.0.2 & R 2020-08-25 r79074

## R CMD check results
There were no ERRORs or WARNINGs.

There were 3 NOTEs (win-builder):

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michael Turchin <michael_turchin@brown.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  MAPITR (13:14, 15:25, 19:42)
  MArginal (2:8)
  ePIstasis (2:17)
  epistasis (14:26, 17:44)
  epistatic (22:34)
  iteratively (17:9)
  phenotypes (16:28)

These words are spelled correctly.


** running examples for arch 'i386' ... [13s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
        user system elapsed
MAPITR 18.79   0.55   12.12

** running examples for arch 'x64' ... [12s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
        user system elapsed
MAPITR 17.77   0.53   10.67

Both of these examples are close to 10s.
