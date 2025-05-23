# 3DBoneMapper

![3DBoneMapper_Logo](https://github.com/Domineuq/3DBoneMapper/assets/114567200/29f67a39-015b-4435-b771-09a44eb44b6b)


### The Identification Tool Using 3D Bone Segmentation and Registration

3DBoneMapper is an automated radiologic identification tool to identify unknown deceased via comparison of postmortem with antemortem computed tomography data. 
The anatomical structures used for identification are the sternal bone and T5 vertebra. However, potentially any bone structure may be used.

Created by the [Forensic Medicine and Imaging Research Group](https://dbe.unibas.ch/en/research/imaging-modelling-diagnosis/forensic-medicine-imaging-research-group/).
If you use it, please cite our publication: 
Neuhaus D, Wittig H, Scheurer E, Lenz C. Computational forensic identification of deceased using 3D bone segmentation and registration. Forensic Science International. 2025;367:112380. doi:10.1016/j.forsciint.2025.112380

# Requirements
+ nibabel
+ numpy
+ SimpleITK
+ os
+ csv

When using FSL FLIRT, FSL has to be installed and the following package is required within the python script:
+ subprocess


# Pipeline
Preparatory steps:
+ Convert the CT DICOM files to NIfTI using dcm2niix
+ Segment the bones (here sternal bone and T5 vertebra) using TotalSegmentator

Script:
+ Import the files
+ Crop the files (improving registration quality and efficiency)
+ Register the files
+ Compute Dice score between postmortem case in question and all registered antemortem cases
+ Highest Dice score indicates a match. If the Dice score is below a certain threshold, a message shows that the identification is not necessarily correct; other identification methods might have to be used.
  

# Usage
+ Download the python script.
+ Define the paths where you have located your antemortem and postmortem data, respectively. 
+ Define the path where you want to save the cropped and registered antemortem data. (The cropping is only performed in the first comparison to improve efficiency. The registered data is deleted after the first postmortem case is processed to save disk space.)

# MIT License
