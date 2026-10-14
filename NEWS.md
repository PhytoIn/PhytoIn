# PhytoIn 0.2.1

- AGB(): the default of `correct.taxon` is now `FALSE`, so AGB() does not access the
  internet unless requested. Use `correct.taxon = TRUE` to standardize taxon names
  with the online TNRS service (previous default).
- AGB(): taxon standardization via the online TNRS service (`correct.taxon = TRUE`)
  now fails gracefully. If the service cannot be reached, a warning is issued and the
  taxon names are used as provided, instead of stopping with an error.
- AGB(): the optional cache directory no longer falls back to the user's home directory
  (R < 4.0.0 now uses the session `tempdir()`); the cache is created only on request
  (`cache = TRUE` or `NULL`).
- AGB(): examples that run during checks no longer require an internet connection
  (`correct.taxon = FALSE`); the online example is wrapped in `\dontrun{}`.
- Documentation of `correct.taxon` and `cache` clarified.

# PhytoIn 0.2.0

- New function summary.param() (S3) with concise and robust output.
- collector.curve() and rarefaction() now include confidence intervals (ribbons) and configurable themes.
- AGB() no longer depends on rappdirs; optional cache via tools::R_user_dir() and BIOMASS::createCache().
- BAplot() fixed (global bindings) and compatible with rectangular plots and individual coordinates.
- New datasets: quadrat2_plot.df, quadrat2_tree.df, quadrat3_rect.df.
- Documentation, examples, and NAMESPACE/DESCRIPTION revised.

## First public release on CRAN

- phytoparam(), summary.param(), plot.param() for phytosociological parameters and diversity indices.
- AGB() for aboveground biomass, carbon, and CO₂ equivalent (wrapper of BIOMASS).
- stratvol() for stratified volume by DBH classes.
- collector.curve() for species accumulation curves.
- rarefaction() for individual-based rarefaction with confidence intervals.
- BAplot() for visualization of basal areas (rectangular plots and individual coordinates).
- Included datasets: quadrat.df, point.df, quadrat2_plot.df, quadrat2_tree.df, quadrat3_rect.df.

# PhytoIn 0.1.0 (development version)

- Initial internal version.
- Prototypes of the core functions for phytosociological analysis.
- Not submitted to CRAN.





