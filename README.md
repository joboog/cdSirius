# cdSirius
**cdSirius** is an implementation of the [Sirius suite of computational mass spectrometry tools](https://bio.informatik.uni-jena.de/software/sirius/) maintained by the Böcker lab at University of Jena within [Compound Discoverer v3.3](https://www.thermofisher.com/ch/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/compound-discoverer-software.html) from ThermoFisher Scientific.  This software is intended to allow annotation of metabolites and other small molecules with molecular formulas, 2D structures, and compound classes based on acquisition of high-resolution, accurate mass MS/MS through data-dependent LC-MS/MS analysis of complex samples.

Implementation of Sirius within Compound Discoverer uses the scripting node functionality within CD as an extensible interface between independent software programs.  Within this implementation, Sirius runs on the host PC and is called as a background service during post-processing through the Sirius API.  Mass spectral data for discrete compounds detected by Compound Discoverer are passed to Sirius for processing, and results are reported back to CD and persisted to the result file for interpretation.  Linkages among CD compounds and Sirius results are maintained within the resulting tables.
## Dependencies
Installation of cdSirius requires a number of dependencies to be installed on the host PC.  _Important_: It is recommended to install all of these dependencies as administrator so that all users will have access.  This is particularly important for the python packages.  Also, the packages should be installed in the base python environment.
* Fully licensed installation of Compound Discoverer v3.3, SP 3
* [Sirius](https://v6.docs.sirius-ms.io/install/) minimum program version 6.0.1
* [Python v.3.11](https://www.python.org/downloads/release/python-3110/).  Note that more recent versions of Python may work but have not been tested.  _Note:_ it is recommended to install Python at `C:/Program Files/Python311/python.exe`.
* [pySirius](https://github.com/sirius-ms/sirius-client-openAPI/tree/master/client-api_python) python package implementing the Sirius REST API for interfacing
* [pyEDS](https://github.com/thermofisherlsms/pyeds/tree/master) python package for programmatic access to mass spectra within Compound Discoverer result files
* [RDKit](https://pypi.org/project/rdkit-pypi/) cheminformatics python package (for implementation of future functionality in cdSirius)
* [pandas](https://pypi.org/project/pandas/) data science python package
## Installation
1. Download the source code and unpack it to a location accessible by all users.  An example might be `C:/python/cdSirius`.
2. Create a new folder at `C:/Program Files/Thermo/Compound Discoverer 3.3/Tools/Scripts/cdSirius` and copy the following files from the source code root directory to the newly created folder.  You will need administrator privileges for this:
   - `node.json`
   - `IMG_16x16.png`
   - `IMG_32x32.png`
3. Edit the node.json file you just copied to correct the paths in `"ScriptProcessorArguments"` for `"ExecutableCommandLineArguments"` and `"ExecutablePath"` and in 
<img width="674" alt="image" src="https://github.com/user-attachments/assets/db45bb6c-96a9-40a2-ac88-aeb97f988c5b" />
