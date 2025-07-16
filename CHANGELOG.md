* 1.6.1 (2025.07.16) - Adding zenodo as a backup data source.
* 1.6.0 (2024.12.21) - Moved data to datacentral.org.au, now downloads maps only as needed.
                       BaseObserver.generate() now allow for the horizon to be set (thanks D. McKenna).
                       Removed case statements to support older Python versions (thanks @ sjoerd-bouma)
* 1.5.4 (2024.04.15) - Added observed_gsm property to `BaseObserver`,
                       Added `hpix2sky` and `sky2hpix` helpers,
                       Now using `query_disc` to find pixels below horizon (neater code).
* 1.5.3 (2024.04.15) - Added `init_gsm()` and `init_observer()` helper functions
* 1.5.2 (2024.03.29) - Fix PyPi distribution files (thanks to D. McKenna)

* 1.5.0 (2024.03.29) - Fix to `GlobalSkyModel16`, correct T_CMB monopole addition (thanks to R. Francis)
* 1.4.1 (2023.04.28) - Added `include_cmb` option to sky models.
* 1.4.0 (2023.04.27) - Added new interpolation method (courtesy R. Braun)
* 1.3.x (2023.03.17) - Reworked horizon masking for Observer classes (thanks to C. Dilullo)
