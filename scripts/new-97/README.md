# Description of the processing

## Processing

### Step 0

- Extract the Pacific mesh mask using `pacific_grid_extraction.py`

### Figure 3

**Mean and OND-97 OOPE anomalies**

- Extract OOPE variables on the Pacific using `pacific_apecosm_exraction.py`
- Extract the 97 period using `extract_pacific_apecosm_97.py`
- Compute the climatology using `pacific_apecosm_climatology.py`
- Plot the figure using `plot_maps_mean_anom_ond.py`

### Figure 4

**Mean and OND-97 anomalies for NEMO/PISCES and Apecosm profiles**

- Extract FORAGE variables on the Pacific using `pacific_apecosm_exraction.py`
- Extract the 97 period using `extract_pacific_apecosm_97.py`
- Compute the climatology using `pacific_apecosm_climatology.py`
- Compute the mean and OND-97 forage anomalies using `compute_forage_mean_ond_97.py`

- Extract NEMO/Pisces variables on the Pacific using `pacific_pisces_extraction.py`
- Compute NEMO/Pisces climatologies using `pacific_pisces_climatology.py`
- Extract data on the 97-period using `extract_pacific_pisces_97.py`
- Compute the mean and OND-97 profile anomalies using `compute_nemo_mean_ond_97.py`

- Plot the figure using `plot_profiles_mean_ond97.py`

### Figure 5

**Hovmoller diagrams for NEMO/PISCES and Apecosm**

- Compute the vertically averaged Pisces fields using `compute_pisces_0-200_mean.py`
- Compute the mean and anomalies for the vertically averaged Pisces fields using `compute_pisces_0-200_anoms.py`

- Compute the surface OOPE anomalies using `compute_apecosm_surface_anoms.py`

- Plot the figure using `plot_all_hovmoller_phys_oope.py`

### Figure 6

**Hovmoller diagrams for Apecosm variables**

- Compute the surface Apecosm (for all the variables) anomalies using `compute_apecosm_surface_anoms.py`
- Extract Apecosm variables on the Pacific using `pacific_apecosm_exraction.py`
- Extract the 97 period using `extract_pacific_apecosm_97.py`
- Compute the climatology using `pacific_apecosm_climatology.py`
- Compute the surface anomalies using `compute_apecosm_surface_anoms.py`
- Plot the figure using `plot_all_hovmoller_apecosm.py`

## Files

### Grid

- **`pacific_grid_extraction.py`: extraction of the Pacific grid**

### Pisces 

- **`compute_nemo_mean_ond_97.py`: compute the OND97 anomalies of NEMO profiles**
- **`compute_pisces_0-200_anoms.py`: compute the anomalies based on the vertically averaged fields**
- **`compute_pisces_0-200_mean.py`: compute the vertical mean for the clim and the raw fields for the 97 Nino**
- `compute_pisces_profile_anomalies.py`: compute Pisces profile anomalies
- **`extract_pacific_pisces_97.py`: extraction of the data for the 97 Nino**
- `compute_pisces_profile_anomalies.py`: compute the profile anomalies for the 97 Nino.
- `plot_pisces_profile_anomalies.py`: plot the profile anomalies for the 97 Nino
- `plot_pisces_surface_hovmoller.py`: plot the Hovmoller diagram for the 97 El Nino
- `extract_pacific_pisces_97.py`: extract NEMO/Pisces for the 97 subperiod
- **`pacific_pisces_climatology.py`: Computation of the climatology for the NEMO/Pisces data**
- **`pacific_pisces_extraction.py`: extraction of NEMO/PISCES data on the Pacifc.**
- `plot_pisces_0-200m_hovmoller.py`: plot Pisces Hovmoller diagram for a single variable averaged between 0 and zmax
- `plot_pisces_profile_anomalies.py`: plot monthly profile anomalies
- `plot_pisces_surface_anomalies.py`: plot monthly surface maps
- `plot_pisces_surface_hovmoller.py`: plot surface hovmoller diagrams

### Apecosm

- **`compute_apecosm_surface_anoms.py`: compute Apecosm anomalies for the 97 Nino**
- `compute_covariance_full_domain_oni.py`: compute lead-lag covariance between detrended OOPE and ONI index
- `compute_covariance_full_domain_pcs.py`: compute 0-lag covariance between detrended OOPE and EOF PCs.
- `compute_covariance_full_domain_tpi.py`: compute lead-lag covariance between detrended OOPE and ONI index
- `compute_eof_pacific_apecosm.py`: compute the EOF for all the Epi. size classes
- `compute_forage_mean_ond_97.py`: compute the mean forage anomalies for the 97 El Nino (averaged over OND months)
- `compute_sardara_time_series.py`: compute the time-series over the saradara domain
- `detrend_pacific_apecosm_anomalies.py`: detrend the full Apecosm time-series of Epi. OOPE
- `extract_pacific_apecosm_97.py`: extract Pacific time-series for the the 97 subperiod
- `pacific_apecosm_anomalies_full_series.py`: computation of the full time-series of Epi. OOPE anomalies
- `pacific_apecosm_climatology.py`: computation of climatology
- `pacific_apecosm_extraction.py`: extraction of Apecosm variables on Pacific
- `plot_all_hovmoller_apecosm.py`: plot Hovmoller diagrams for all the Apecosm variables
- `plot_apecosm_yearly_surface_anomalies.py`: plot yearly OOPE anomalies for 97 Nino
- `plot_apecosm_surface_anomalies.py`: plot monthly anomalies maps
- `plot_apecosm_surface_hovmoller.py`: plot Hovmoller anomalies for Apecosm variables
- `plot_eof_covariances.py`: plot EOF covariance maps (whole Pacific domain)
- `plot_eof_maps.py`: plot EOF maps
- `plot_eof_pcs.py`: plot EOF pcs
- `plot_integrated_adv_hovmoller_apecosm.py`: plot the Hov diagrams for integrated adv. trends
- `plot_maps_mean_anom_ond.py`: plot the maps of OOPE mean and OND anomalies
- 

### Apecosm + Pisces  

- **`plot_all_hovmoller_phys_oope.py`: plot Hovmoller diagrams for NEMO and Apecosm (biomass) variables**
- **`plot_profiles_mean_ond97.py`: plot mean and profile anomalies for NEMO/Pisces variables**

## Misc

- `plot_oni_index_subperiod.py`: plot the ONI uindex for the 97 subperiod


