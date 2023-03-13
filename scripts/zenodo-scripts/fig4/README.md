To reconstruct figure 4:

- `pisces_extract_equatorial_values_full_series.py`: extract NEMO/Pisces data along the equator. Need to be applied to all the variables (`thetao`, `uo`, `PLK`)
- `pisces_compute_hovmoller_composites_strong_ninos.py`: computes the strong Nino composites. Need to be applied to all the variables (`thetao`, `uo`, `PLK`)
- `apecosm_extract_equatorial_values_full_series.py`: extract Apecosm data along the equator. Need to be applied to all the variables (`OOPE`, `gamma1`, ...)
- `apecosm_compute_hovmoller_composites_strong_ninos.py`: computes the strong Nino composites. Need to be applied to all the variables (`OOPE`, `gamma1`, ...)
- `fig4.py`: Draws the figure


