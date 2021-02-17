.. _description:

#################
ERSEM area models
#################

ERSEM is used in a variety of different models, the table bellow outlines several
of these models.


+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| Area                          | Includes       | Spatial Scale                    | Quality (data used)               |
+                               +----------+-----+-----------------------+----------+                                   +
|                               | Hydro    | BGC | Domain                | Res.(km) |                                   |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| Water column (GOTM-ERSEM)     | X        | X   | Any                   | N/A      | 4,6,8,14,16,22                    |
|                               |          |     | location              |          |                                   |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| Global (NEMO)                 | X        | X   | Global                | 111      | 8,11,19,20,21                     |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| North Atlantic (NEMO)         | X        | X   | 205-BON               | 28       | 8,11,19,20,21                     |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| North Atlantic (NEMO shelf)   | X        | X   | 20N-80N               | 7        | 8,11,19,20,21                     |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| Irish and Celtic Seas (FVCOM) | X        | X   | 41.8-56.8N9, .6-Z.SW  | 3.5      | 7,23,24,25                        |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| NW European shelf (NEMO)      |          | X   | 40N-6SN, ZOW-13E      | 7        | 3,8,11,12,13,14,19                |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| NW European shelf (POLCOMS)   | X        | X   | 40N-6SN, 20W-13E      | 10       | 3,8,9,10,11,12,13,14,15,16,17,19  |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| NW European shelf (POLCOMS    | X        | X   | 40N-6S                | 10       | 2,3,8,9,10,12,13,14,15,17         |
| + data assimilation)          |          |     | N,ZOW-13E             |          |                                   |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| AMM7NW European-Baltric-GCOMS | X        | X   | 46.4-63N, 17.SW-13E   | 10       | 3,8,11,14,17,18                   |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| WEC                           | X        | X   | 41.5-50.9N, 7.6-1.2SW | 1.9      | 2,4,8,12,13,14                    |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+
| WEC (data assimilation)       | X        | X   | 41.5-50.9N, 7.6-1.2SW | 7        | 1,2,4,8,12,13,14                  |
+-------------------------------+----------+-----+-----------------------+----------+-----------------------------------+


Key
~~~

.. hlist::
    :columns: 3

    - 1. Glocolour ocean colour
    - 2. ESA CCI ocean colour
    - 3. North Sea Project cru1ses
    - 4. Western Channel Observatory
    - 5. Tide sauses 
    - 6. Smart Buoy data 
    - 7. NCEP reanalysis
    - 8. World Ocean Atlas
    - 9. ICES temperature andsalinity
    - 10. ICES nutrients and chlorophyll
    - 11. IPCC climate forcing (HADGEMI,PSL, ECHAM)
    - 12. EA riverine nutrients data 
    - 13. Europeanriver data
    - 14. ECMWF reanalysis meteorolosv
    - 15. DMIreanalysis meteorology
    - 16. ESSC ocean reanalysis
    - 17. GLORYS reanalysis
    - 18. GLOBAL NEWS river data
    - 19. GLODAP
    - 20. DFS atmospheric reanalysis
    - 21. GEM/GLORI river data
    - 22. BATS data
    - 23. HF Radar and CTD Data.
    - 24. NTSLF data
    - 25. CEH river discharge and temperature.


.. note::
    The table of area models is not yet complete, please ERSEM's issue
    tracking for updates.

Existing uses:
~~~~~~~~~~~~~~

- **Natural resources:** understanding biogeochemical cycles, biodiversity 
  and valuation of ecosystem services.
- **Resilience to environmental hazards:** eutrophication, microplastics, 
  fishing pressure, invasive species and harmful algal blooms. Impacts 
  of offshore renewable energy and ecological risks of carbon capture 
  and storage.
- **Environmental change:** climate change impacts, ocean acidification, 
  multiple stressors, and sustainable fisheries.
- Run operationally by the UK Met Office to predict water quality.
- Estimation of the carbon budget of the UK shelf.


Potential new uses:
~~~~~~~~~~~~~~~~~~~

- Understanding shelf seas carbon ("blue carbon") and nutrient budgets 
  {past and future climate).
- Expansion to represent biodiversity-relevant processes over a range 
  of spatial and temporal scales, and simulate changes in function in 
  the context of ecosystem services.
- Implementation and testing scalable models of differing complexity.

Key modelling issues:
~~~~~~~~~~~~~~~~~~~~~

- Setting inputs (parameterisation} and testing outputs against real 
  data (calibration) is an essential, but resource-intensive and 
  on-going process to ensure quality and improve predictions. 
  Understanding the impact of changing inputs on the outputs from the 
  models (sensitivity) and the effect of uncertainty in model parameters 
  on robustness of model predictions.
- Challenge to assess model capability with respect to seasonal variability, 
  long-term changes, regime shift and tipping points due to limitations of 
  the data available.
- Complexity of model leads to a need for significant interpretation and 
  explanation for stakeholders.
- Potential mismatch between scales of model output and data sets. 
- Significant expertise needed to operate system and high performance parallel 
  computing facility required for three-dimensional full scale simulations. 
