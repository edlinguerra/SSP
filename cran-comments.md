## Test environments
* local Windows 10 Pro, R 3.6.1
* win-builder (devel, release, older)
* ubuntu 18.04.3 LTS, R 3.6.2
* Ubuntu 16.04.6 LTS (on travis-ci: release, older and devel)
* Mac OS X (on travis-ci: release, older and devel)

## R CMD check results
Copy/Paste of the R CMD check results

-- R CMD check results ------------------------------------------ SSP 1.0.1 ----
Duration: 1m 25.7s

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

## Downstream dependencies
I have also run R CMD check on downstream dependencies of SSP without any errors, warnings or notes


## Resubmission
This is a resubmission. In this version I have:

* Added key references and respective doi, as suggested, in the Description field of the DESCRIPTION file.
* Incorporated the article doi that describes the SSP methods in all the help files.
