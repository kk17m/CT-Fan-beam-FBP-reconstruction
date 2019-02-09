# Computed-tomography Fan-beam FBP reconstruction

This repository contains CT image reconstruction using Fan-beam Filtered back-projection. Appropriate weighting measures like Differential and Parker weighting can be applied. The reconstruction algorithm is applicable to short scan protocol as well.  


## FILE CONTENTS

* **FFBP_Weighted.m :** *Fan-beam Filtered back-projection function for two dimensional reconstruction.*

* **SAMPLE_FFBP.m :** *Sample script to execute the Fan-beam FBP function with appropriate reconstruction parameters. This example script can be used to reconstruct the sample sinogram file " Sample_sinogram.sino ".*

* **designFilter.m :** *Use it to design the filter with window of choice. Certain general windows are already incorporated like the shepp-logan, cosine, hann and hamming filter.*

* **fbp2_window.m :** *Filter windows as used by Jeffrey A Fessler, "Michigan Image Reconstruction Toolbox,"*

* **fbp_sino_filter.m :** *Filtering (Curved/ Arc detector) as used by Jeffrey A Fessler, "Michigan Image Reconstruction Toolbox,"*
                  
* **filterProjections.m :** *Filtering as in iradon function (Copyright 1993-2013 The MathWorks, Inc.).*             

### FILE CONTENTS

Taguchi, K. Image Reconstruction Algorithms for X-Ray CT, 189–216. Medical Imaging Technology and Applications Troy Farncombe and Kris Iniewski CRC Press 2013.
