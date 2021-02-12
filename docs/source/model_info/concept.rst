.. _concept:

##############
What is ERSEM?
##############


ERSEM (European Regional Seas Ecosystem Model) :cite:`ERSEM,code` is a model for the 
pelagic and benthic
dynamics of the lower trophic level of the marine ecosystem. From its very beginning 
in the early nineties :cite:`Baretta1995` it has been based on the functional group 
approach, concentrating in its food-web structure on the function of a group of 
organisms in the ecosystem often associated with size rather than taxa. A 
particularity of the model with respect to the greater part of other biogeochemical 
models and ecosystem models for the lower trophic level of the marine food-web
is the full resolution of the major biogeochemical cycles of carbon, nitrogen, 
phosphorus, silicon and optionally iron with fully dynamic stochiometrioc relationships 
in its states. Additionally, with its main components phytoplankton, zooplankton, 
bacteria, dissolved organic matter and particulate matter it fully resolves the 
microbial part of the marine food-web. Contrary to its name which derives from its
original purpose of an ecosystem model specifically developed for the shelf seas of 
the North-West European shelf it has been extended over the recent years in applications 
to a variety of shelf seas and coastal ocean regions as well entire ocean basins 
and the global ocean.


The ecosystem in ERSEM  is divided into functional types, which are further 
subdivided by traits such as size. In the pelagic, ERSEM by default 
distinguishes:

* 4 types of phytoplankton: diatoms, picophytoplankton, nanophytoplankton,
  microphytoplankton
* 3 types of zooplankton: nanoflagellates, microzooplankton, mesozooplankon
* bacteria

The benthic system includes:

* 3 types of infauna: meiofauna, suspension feeders, deposit feeders
* 2 types of bacteria: aerobic and anaerobic

In addition, ERSEM tracks the concentrations of phosphate, nitrate, ammonium,
silicate, iron, oxygen, dissolved inorganic carbon and alkalinity in the
pelagic and in sediment porewaters. It includes several classes of particulate
and dissolved organic matter in pelagic and sediment. A carbonate system module
calculates pH and calcium carbonate saturation.

Through `FABM <http://fabm.net>`__, ERSEM can be coupled to a wide range of
hydrodynamic models including NEMO, FVCOM, ROMS and GOTM. FABM also makes it
easy to customise the default set of functional types described above, and to
combine ERSEM with modules representing other parts of the ecosystem, including
fish communities, shellfish, seagrass meadows and spectraly resolved irradiance.

.. figure:: ../../images/ERSEM.png
   :alt: ERSEM diagram
   :width: 100.0%

