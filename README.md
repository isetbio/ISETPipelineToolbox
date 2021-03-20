## Image Reconstruction Toolbox
This repository contains the Bayesian image reconstruction algorithm we used in our project: Bayesian Image Reconstruction from Cone Mosaic Signal. For details, visit the main [Project Page](https://github.com/isetbio/ISETImagePipeline).

## Structure of the Repository
Please also see the Methods section of our paper. For actual usage, the best way is to start with the instructions in the [Project Page](https://github.com/isetbio/ISETImagePipeline). 
```
...
├── coneMosaic                      
    ├── ConeResponse.m              # Base Class, Wrapper of ISETBio, compute cone responses to natural images
    ├── ConeResponsePeripheral.m    # Extend the base class to use peripheral cone mosaic and optics
    ├── PeripheralModel.m           # Static function that returns optical model (PSF) at different visual eccentricity
    ├── ConeResponseCmosaic.m       # Extend the base class to use the newly updated cone mosaic 
├── imageEstimator                  # Deprecated; Regression-based method for image reconstruction    
    ├── ...
├── patchEstimator                  # The main algorithm used in our paper
    ├── PatchEstimator.m            # Base class for the estimator
    ├── GaussianPatch...            # Implement the Gaussian image prior
    ├── SparsePatch...              # Implement the sparse coding image prior
    ├── Poisson...                  # Implement the Poisson likelihood function
├── imageHelper                     # Some helper function for image processing
    ├── MarkovPrior.m               # Gaussian prior for which we have control over its spatial and chromatic correlation
    ├── computeBasisPCA.m           # PCA on dataset of natural images
    ├── whitening.m                 # Apply the whitening transformation to images
    ├── sampleImage.m               # Crop image to a specific size
    ├── ... 
├── visualHelper                    # Some helper functions for visualization
    ├── ... 
├── fminlbfgs.m                     # MATLAB community implementation of limited memory BFGS.
...
```
Please see here for the original source of [fminlbfgs](https://www.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer).

## Contact
Feel free to contact me if you any questions or comments at   
**lingqiz at sas dot upenn dot edu**