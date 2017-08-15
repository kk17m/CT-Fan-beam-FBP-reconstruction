# Computed-tomography Fan-beam FBP reconstruction

This repository contains CT image reconstruction using Fan-beam Filtered back-projection. The reconstruction algorithm is applicable to short scan protocol as well. Appropriate sinogram weighting, like differential weighting and parker weighting, can be applied depending on the requirements.   


### FILE CONTENTS

* **FFBP_Weighted.m :** *Fan-beam Filtered back-projection function for two dimensional (slice) reconstruction.*

* **SAMPLE_FFBP.m :** *Sample script to execute the Fan-beam FBP function with appropriate reconstruction parameters. This example script can be   used to reconstruct the sample sinogram file " Sample_sinogram.sino ".*

* **designFilter.m :** *Use it to design the filter with window of choice. Certain general windows are already incorporated like the shepp-logan, cosine, hann and hamming filter.*

* **fbp2_window.m :** *Filter windows as used by Jeffrey A Fessler, "Michigan Image Reconstruction Toolbox,"*

* **fbp_sino_filter.m :** *Filtering (Curved/ Arc detector) as used by Jeffrey A Fessler, "Michigan Image Reconstruction Toolbox,"*
                  
* **filterProjections.m :** *Filtering as in iradon function (Copyright 1993-2013 The MathWorks, Inc.).*             
