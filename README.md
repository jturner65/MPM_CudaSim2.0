# MPM_CudaSim2.0
Update of MPM Cuda project with recompiled kernel, and updated to Jcuda for CUDA 9.2

I've merged a few of the kernels and fixed a few bugs in the calculations, restructured how the kernel calls are built and called, and generally cleaned up and streamlined the java code.  I've also changed some of the default values so that the sim isn't so brittle, which enabled me to increase the timestep by a factor of 10.