# MPM_CudaSim2.0
4/2025 Update : I've been expanding the simulator a lot, along with adding (and rewriting) the old multi-threaded cpu implementation. 

12/2022 Update : I've expanded what the sim supports - now any number of snowballs, up to 20, can be made, and they will all target a random point in a fairly tight zone around their mutual global COM, so they make delightful splashes around each other. I've also made rebuilding the sim more efficient, only rebuilding the objects or the cuda kernel when necessary due to modifications of dependent values.  Currently working on using a vertex buffer to track particle positions, and having the cuda kernel directly modify the values within, so that no copying is necessary. Should substantially increase framerate.

Update of MPM Cuda project with recompiled kernel, and updated to Jcuda for CUDA 9.2

I've merged a few of the kernels and fixed a few bugs in the calculations, restructured how the kernel calls are built and called, and generally cleaned up and streamlined the java code.  I've also changed some of the default values so that the sim isn't so brittle, which enabled me to increase the timestep by a factor of 10.

Here's a link to a video : 

https://www.dropbox.com/s/8idue9mlvg6lee0/MPM%20Snow%20JCuda.mp4?dl=0
