.. _description:

#################
ERSEM area models
#################

ERSEM is a planktonic ecosystem model which has been coupled to a number 
of different hydrodynamic models. It describes the biogeochemical cycling 
of carbon and the nutrients nitrogen, phosphorous, silicon, oxygen and iron. 
The ecosystem is subdivided into three functional types: producers 
(phytoplankton), decomposers (bacteria} and consumers (zooplankton), and 
then further subdivided by trait (size, silica uptake) to create a foodweb. 
Physiological and population processes are included in the descriptions of 
functional group dynamics. Four phytoplankton, three zooplankton and one 
bacteria are represented, along with the cycling of carbon, nitrogen, 
phosphorous, silicon and oxygen through pelagic and benthic ecosystems.


+-----------+----------+-----+-----------+----------+-----------+
| Area      | Includes       | Spatial              | Quality   |
| Modelled  |                | Scale                | (data     |
|           |                |                      | used)     |
+           +----------+-----+-----------+----------+           +
|           | Hydro    | BGC | Domain    | Res.(km) |           |
+-----------+----------+-----+-----------+----------+-----------+
| Water     | X        | X   | Any       | N/A      | 4,6,8     |
| column    |          |     | location  |          | ,14,16,22 |
| (GO       |          |     |           |          |           |
| TM-ERSEM) |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| Glo       | X        | X   | Global    | 111      | 8,11,     |
| bal(NEMO) |          |     |           |          | 19,20,21, |
+-----------+----------+-----+-----------+----------+-----------+
| North     | X        | X   | 205-BON,  | 28       | 8,11      |
| Atlantic  |          |     |           |          | ,19,20,21 |
| (NEMO)    |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| North     | X        | X   | 20N-80N,  | 7        | 8,11      |
| Atlantic  |          |     |           |          | ,19,20,21 |
| (NE       |          |     |           |          |           |
| Mo-shelf} |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| Irish     | X        | X   | 41.       | 3.5      | 7         |
| andCeltic |          |     | 8-56.8N9, |          | ,23,24,25 |
| Seas      |          |     | .6-Z.SW   |          |           |
| (FVCOM)   |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| NW        |          | X   | 40N-6S    | 7        | 3,8,11,   |
| European  |          |     | N,ZOW-13E |          | 12,       |
| shelf     |          |     |           |          | 13,14,19  |
| (NEMO)    |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| NW        | X        | X   | 40N-6S    | 10       | 3,8,9     |
| European  |          |     | N,20W-13E |          | ,10,11,12 |
| shelf     |          |     |           |          | ,13,14,15 |
| (POLCOMS) |          |     |           |          | ,16,17,19 |
+-----------+----------+-----+-----------+----------+-----------+
| NW        | X        | X   | 40N-6S    | 10       | 2,        |
| European  |          |     | N,ZOW-13E |          | 3,8,9     |
| shelf     |          |     |           |          | ,10,12,13 |
| (POLCOMS  |          |     |           |          | ,14,15,17 |
| + data    |          |     |           |          |           |
| assi      |          |     |           |          |           |
| milation) |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| AMM7NW    | X        | X   | 46.4-63N, | 10       | 3,8,11,   |
| Eur       |          |     | 17.SW-13E |          | 14,17,18  |
| opean-Bal |          |     |           |          |           |
| tic-GCOMS |          |     |           |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| WEC       | X        | X   | 41        | 1.9      | 2.4,8     |
|           |          |     | .5-50.9N, |          | ,12.13,14 |
|           |          |     | 7.6-1.2SW |          |           |
+-----------+----------+-----+-----------+----------+-----------+
| WEC (data | X        | X   | 41        | 7        | 1,2,4,8   |
| assi      |          |     | .5-50.9N, |          | ,12,13,14 |
| milation) |          |     | 7.6-1.2SW |          |           |
+-----------+----------+-----+-----------+----------+-----------+




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
