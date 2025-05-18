# ClockGUV-Model

## Requirements
Code requires MATLAB version R2021a or greater with packages:
(1 of 5) Image Processing Toolbox, version 11.3 or greater
(2 of 5) Curve Fitting Toolbox, version 3.5.13 or greater
(3 of 5) Signal Processing Toolbox, version 8.6 or greater
(4 of 5) Statistics and Machine Learning Toolbox, version 12.1 or greater.
(5 of 5) Computer Vision Toolbox, version 10.0 or greater.

(Note: not all dependencies may be strictly required but were used in the testing environment)

## Description
Models the behavior of the post-translation oscillator (PTO) from the circadian clock of the cyanobacteria, Synechococcus elongatus, inside a giant unilamellar vesicle (GUV). The model simulates the encapsulation of the PTO (comprised of three Kai clock proteins) inside various GUVs sizes.
There are two versions of this model, one simulating the environment within a GUV used in our experiments (Clock Model - GUV Model) and one simulating the behavior of a transcription-translation feedback loop (TTFL) system (Clock Model - TTFL Model). Clock behavior is determined based on experimental data with lookup tables incorporated.  

## Instructions
1. Run either "ClockModel_GUVModel.m" or "ClockModel_TTFLModel.m" from their relevant folders. Requires the subfolder "Amplitude Simulation" and "Period Simulation" as this contains data necessary for the code to run.
2. Model results will be saved to a "Model_Results.mat" file in the same directory. A sample result file is included.

## Additional Notes
The model .m files contain modifiable input parameters in the beginning.
