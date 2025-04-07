# cdSirius
**cdSirius** is an implementation of the [Sirius suite of computational mass spectrometry tools](https://bio.informatik.uni-jena.de/software/sirius/) maintained by the Böcker lab at University of Jena within [Compound Discoverer v3.3](https://www.thermofisher.com/ch/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/compound-discoverer-software.html) from ThermoFisher Scientific.  This software is intended to allow annotation of metabolites and other small molecules with molecular formulas, 2D structures, and compound classes based on acquisition of high-resolution, accurate mass MS/MS through data-dependent LC-MS/MS analysis of complex samples.  Note that cdSirius will only work with high-resolution MS and MS/MS data, and has only been tested with data-dependent analysis (DDA) results.  It may work with DIA data, but no guarantees are given.

Implementation of Sirius within Compound Discoverer uses the scripting node functionality within CD as an extensible interface between independent software programs.  Within this implementation, Sirius runs on the host PC and is called as a background service during post-processing through the Sirius API.  Mass spectral data for discrete compounds detected by Compound Discoverer are passed to Sirius for processing, and results are reported back to CD and persisted to the result file for viewing and interpretation.  Linkages among CD compounds and Sirius results are maintained within the resulting tables.
## Dependencies
Installation of cdSirius requires a number of dependencies to be installed on the host PC.  _Important_: It is recommended to install all of these dependencies as administrator so that all users will have access.  This is particularly important for the python packages.  Also, the packages should be installed in the base python environment.
* Fully licensed installation of Compound Discoverer v3.3, SP 3
* [Sirius](https://v6.docs.sirius-ms.io/install/) minimum program version 6.0.1
* [Python v.3.11](https://www.python.org/downloads/release/python-3110/).  Note that more recent versions of Python may work but have not been tested.  _Note:_ it is recommended to install Python at `C:/Program Files/Python311/python.exe`.
* [PySirius](https://github.com/sirius-ms/sirius-client-openAPI/tree/master/client-api_python) python package implementing the Sirius REST API for interfacing
* [PyEDS](https://github.com/thermofisherlsms/pyeds/tree/master) python package for programmatic access to mass spectra within Compound Discoverer result files
* [RDKit](https://pypi.org/project/rdkit-pypi/) cheminformatics python package (for implementation of future functionality in cdSirius)
* [pandas](https://pypi.org/project/pandas/) data science python package
## Sirius user account
To run Sirius with CSI:FingerID and CANOPUS functionality, you will need a Sirius user account.  You can create a user account as described in the [Sirius Wiki](https://v6.docs.sirius-ms.io/account-and-license/).  The username and password will be used in the cdSirius node for authentication.  Do not re-use a sensitive password for this user account, as the password will not be encrypted and will be visible in plain text within the CD method editor. _Note_: If you are an academic user, you will qualify for a free account, but you must use your institutional email address when you register your account.  
## Installation
1. After dependencies above are fulfilled, download the source code and unpack it to a location accessible by all users.  An example might be `C:/python/cdSirius`.
2. Create a new folder at `C:/Program Files/Thermo/Compound Discoverer 3.3/Tools/Scripts/cdSirius` and copy the following files from the source code root directory to the newly created folder.  You will need administrator privileges for this:
   - `node.json`
   - `IMG_16x16.png`
   - `IMG_32x32.png`
3. Edit the node.json file you just copied to correct the paths in lines 19, 20, and 30 according to your local installation.

   <img width="500" alt="image" src="https://github.com/user-attachments/assets/9c965c4d-73cd-4ccb-9d09-5d451a725f1f" />

   **Figure 1.** node.json file section with paths to relevant locations
5. Launch Compound Discoverer 3.3 SP3 and navigate to the Help -> License Manager dialogue.  Run "Scan for Missing Features":

   <img width="693" alt="image" src="https://github.com/user-attachments/assets/1b5c8aa4-cf06-4425-9251-429dfb610424" />

   **Figure 2.** Scanning for missing features within the CD license manager dialogue

6.  Close and re-start Compound Discoverer to complete installation and allow new nodes to be registered.
## Using cdSirius within a Compound Discoverer workflow
The cdSirius node is a post-processing node that can be appended to an existing full processing node, or it can be included in a "reprocessing" workflow to retrospectively add Sirius results to the cdResult file.  Either way, you will find the new Sirius node within the Workflow Editor Node menu, in the _10. Post-Processing_ sub-menu:
   
      <img width="194" alt="image" src="https://github.com/user-attachments/assets/c5a6a15b-ec71-4182-91bf-ca05b04d28fd" />

After adding the node to the workflow, the processing configuration dialogue is available for editing.  Default parameters are provided for most settings:

      <img width="567" alt="image" src="https://github.com/user-attachments/assets/76b6c9a8-f21d-40dd-8cb0-7f5e6265e309" />

      **Figure 3.** Sirius node parameter configuration.

### cdSirius parameter settings
1.  **Sirius Program Settings:**  These are global settings for the Sirius program service.
   - <ins>Sirius Program Path</ins>: Here you can set the program path for the Sirius executable.  The default should be suitable for most installations, but can be changed for non-standard installations.
   - <ins>Save Sirius Result</ins>: Setting this to `True` will enable the .sirius workspace to be persisted as a permanent file, saved to the same directory with the cdResult file.  This is useful if e.g. you plan to re-open and analyze data using the Sirius GUI at a later time.
   - <ins>Sirius Username</ins>: Your username as configured on the BrightGiant web portal as described above
   - <ins>Sirius Password</ins>: Your password as configured on the BrightGiant web portal as described above (**Note**: Do not re-use sensitive passwords here as this information is not encrypted or hidden)
2.  **Compound Selection Settings:**  These settings modify which Compounds within the CD processing environment are selected for submission to Sirius
   - <ins>Checked Feature Status Handling</ins>: This is a switch that allows for down-selection of only compounds of interest when cdSirius is used in reprocessing of existing cdResult files.  Selecting "Checked" will pass only "checked" compounds to Sirius for calculation. The default of "All" will pass all non-background compounds (which have not been pre-filtered within the workflow by peak quality thresholds) to the Sirius service.
   - <ins>Peak Quality Threshold</ins>: This threshold allows a more stringent selection of chromatographic peaks for Sirius processing.  See the CD documentation for an explanation of peak quality metrics.
   - <ins>Maximum MW</ins>: Sirius computation becomes extremely slow for large molecules, so it is best to limit the upper MW range to only those of interest for a particular analysis.
3.  **Molecular Formula Prediction Settings:**  This set of parameters controls the "base" molecular formula calculation functions within Sirius and is necessary for all further processing
   - <ins>Predict Formulas</ins>: Must be set to "True" for any processing to occur.  Set to False only for program debugging purposes.
   - 


