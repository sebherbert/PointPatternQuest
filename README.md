
# GUIDELINES for using the spatial dispersion analysis codes

Sebastien Herbert herbert.sebastien@gmail.com

***

Required Matlab libraries
- https://fr.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids



***
## For cell type dispersion in a dynamic patern:

### Preprocess the data files =>  `extractDuplicateAndFormatDyn.m`

Reformat the cell types exported in .csv format by Imaris into a single MATLAB file using the `extractDuplicateAndFormatDyn.m` function. The processing also checks that there are no duplicates in the files.

- This function will ask for the 3 expected cell types. Look at the csv files for an example. (Attention: there can be more than 1 input file for the type 1 population)
- Each file will be converted into a MATLAB table concatenated with the others and tested for duplicates. Duplicates are of 3 kinds:
  - duplicated cells are type2 and mother
  - duplicated cells are of the same cell type
  - duplicated cells are of different cell types (not mother) => only this situation will bring the function to throw an error message prompting the user to check which format should used (usually by going back to Imaris).
- If the files were properly set and provided, the output should be an 8 variable table by N lines (being the number of cells). Variables are :
  - PositionX
  - PositionY
  - PositionZ
  - Unit
  - time
  - ID
  - cellType
  - oldID (original Imaris cell ID, not always unique since the cells were concatenated).
- The table is saved under the name **liveDataCurated.mat** and contains the fullDataLive table

### Launch the analysis of the dynamic pattern => `pointDynPatternMain.m`

This function will both handle the nearest neighbour analysis with simulation of the possible spatial distribution effects and their display. This function will automatically launch a parallel pool to accelerate calculations. The resulting analysis is stored under the **completeAnalysis.mat** file. By modifying the parameters in the code you can specify:
- which populations should be tested against which (For example mother against mother by permuting into the type 1+2 population).
- the number of random permutations for each model (10000 by default)
- the ranges of R and S to model
- that analyses were already calculated (if you just want to re-plot the displays), the algorithm will then ask you to select the analysis file (**completeAnalysis.mat**).

Once you are satisfied with the parameters to use, simply call the function and select the **livedataCurated.mat** created at the previous step.
By default, this function will run output 4 images in .fig format and 4  data files in .mat format.

Images: *Name is automatically appended with the permutation population*
- **RMSEisomap_permutIn_type1type2.fig**: RMSE map with isolines separation
- **RMSEmap_permutIn_type1type2.fig**: RMSE map with interpolation between sampled positions
- **allDistFreq_type1type2.fig**: All distances frequency between the source and target populations
- **averageDistFreq_type1type2.fig**: All distances frequency between the source and target populations averaged per individual simulation

DataFiles:
- **completeAnalysis.mat**: Contains the next three files but can be quite big.
- **NNanalysis.mat**: Contains only the processed and analysed dispersion results
- **NNdispersion.mat**: Contains only the dispersion results prior to analysis
- **NNdistances.mat**: Contains only the distances