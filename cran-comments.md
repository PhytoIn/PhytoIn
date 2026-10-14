## Resubmission of an archived package

PhytoIn was archived on 2026-09-15 ("policy violation"). The maintainer did not
receive an e-mail about the problem before the archival (the maintainer address is
active: it received the CRAN maintainer check of 2026-07-14). From the CRAN check
records, the only problem reported was in the "donttest" additional check
(r-devel, 2026-08-24): the AGB() example used the TNRS web service through
BIOMASS::correctTaxo(); when the service could not be reached, AGB() stopped with
an error instead of failing gracefully, contrary to the CRAN policy on Internet
resources. If the archival had any other reason, we would be grateful to be told,
and we will fix it promptly.

Changes in version 0.2.1 (the archived version was 0.2.0):

* AGB() no longer accesses the Internet by default (correct.taxon = FALSE is now
  the default).
* When requested (correct.taxon = TRUE), AGB() handles an unavailable TNRS service: it gives an informative warning
  and uses the taxon names as provided (no error).
* Examples run by R CMD check (including --run-donttest) no longer use the
  Internet (correct.taxon = FALSE); the online example is in \dontrun{}.
* The optional cache directory is only created when the user requests it
  (cache = TRUE or NULL, as agreed during the first submission in 2025) and never
  falls back to the user's home directory.

## Test environments

* GitHub Actions, R CMD check --as-cran --run-donttest:
  - windows-latest: R 4.6.1 and R-devel (2026-10-13 r90668)
  - macos-latest: R 4.6.1
  - ubuntu-latest: R-devel (2026-10-12 r90666), R 4.6.1 and R 4.5.3
  All: Status OK.
* win-builder: R 4.6.1 (1 NOTE, see below).

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  - New submission / Package was archived on CRAN: this is the resubmission
    of the archived package, as explained above.
  - Possibly misspelled words in DESCRIPTION (DBH, FAO, Heringer, IBGE,
    Phytosociology, al, et, brasileira, da, phytosociological, quadrat, ...):
    these are technical terms, acronyms, author names and words of the
    Portuguese title of a cited reference; they are correct.
