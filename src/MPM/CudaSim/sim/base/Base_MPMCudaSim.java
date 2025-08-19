package MPM.CudaSim.sim.base;

import static jcuda.driver.JCudaDriver.cuCtxCreate;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import static jcuda.driver.JCudaDriver.cuMemAlloc;
import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;
import static jcuda.driver.JCudaDriver.cuMemcpyHtoD;
import static jcuda.driver.JCudaDriver.cuMemsetD32;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleLoad;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import MPM.BaseSim.sim.Base_MPMSim;
import MPM.BaseSim.sim.SimResetProcess;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Render_Interface.IGraphicsAppInterface;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import jcuda.runtime.JCuda;

/**
 * abstract class describing a world for a CUDA-driven MPM simulation.
 * Called by sim window to execute simulation and render results.
 * Instancing classes can hold different configurations/initialization setups
 */
public abstract class Base_MPMCudaSim extends Base_MPMSim{    
    
    /**
     * CUDA kernel file name
     */
    
    private String ptxFileName = "data/MPM_CUDA_Sim_New.ptx";    
    //private String ptxFileName = "data/MPM_ABS_Sim_With_Init.ptx";        
    
    ////////////////////////////////////////////////////
    // Maps holding CUDA function pointers and parameter pointers
    // This facilitates CUDA calcs and access to appropriately configured CUDA args
    private TreeMap<String, Pointer> kernelParams;
    private TreeMap<String, CUfunction> cuFuncs;
    /**
     * Lists of names of kernel functions to perform for MPM algorithm.   
     */
    private String[] CUFileFuncNames = new String[]{
            "projectToGridandComputeForces",
            "projectToGridInit",
            "computeVol",
            "updPartVelocities",
            "compGridVelocities",
            "partCollAndUpdPos",
            "gridCollisions",
            "clearGrid",
            "updDeformationGradient"};    
    /**
     * These are the individual kernel functions to execute in order for initialization of simulation
     */
    private String[] initStepFuncKeys = new String[]{
            "clearGrid",
            "projectToGridInit",
            "computeVol",
            "compGridVelocities",
            "gridCollisions",
            "updDeformationGradient",
            "updPartVelocities",
            "partCollAndUpdPos"}; 
    /**
     * These are the individual kernel functions to execute in order for each sim step
     */
    private String[] simStepFuncKeys = new String[]{
            "clearGrid",
            "projectToGridandComputeForces",
            "compGridVelocities",
            "gridCollisions",
            "updDeformationGradient",
            "updPartVelocities",
            "partCollAndUpdPos"}; 
   
    private HashMap<String, int[][]> funcGridDimAndMemSize;
    
    ////////////////////////////////////////////////////
    // CUDA references
    private CUdevice dev;
    private CUcontext context;
    private CUmodule module;
    //private CUgraphicsResource pCudaResource;
    // CUDA Device ptr constructions
    private CUdeviceptr devPtrPartMass, devPtrPartVolume;
    private CUdeviceptr[] devPtrPartPos, devPtrPartVel; 
    private CUdeviceptr[] devPtrGridVel, devPtrGridNewVel, devPtrGridForce;
    private CUdeviceptr[][] devPtrPartElasticF, devPtrPartPlasticF;    
    private CUdeviceptr devPtrGridMass;
    // CUDA calc helper variables

    /**
     * # of cuda blocks to use for particles
     */
    private int numBlocksParticles;
    /**
     * # of cuda blocks to use for grid
     */
    private int numBlocksGrid;
    /**
     * # parts * size of float & num grid cells * size float 
     */
    protected long numPartsFloatSz, numGridFloatSz;
    /**
     * # of cuda threads per block
     */
    private final int numCUDAThreads = 1024;
    private final int[] blkThdDims = new int[] {numCUDAThreads, 1, 1};
    private int[] partGridDims;
    private int[] part4GridDims;
    private int[] gridGridDims;
    private int[] shrdMemSize = new int[] {0};
    private final int[] partVelShrdMemSize = new int[] {numCUDAThreads*6*Sizeof.FLOAT};
    
    /**
     * Raw initial particle values - position and velocity
     */
    private TreeMap<String, ArrayList<float[]>> partVals;
    
    ////////////////////////////////////////////////////
    //representations for rendering
    /**
     * local representation of particle vector quantities for rendering
     */
    private float[][] hostPartPos, hostPartVel;
    /**
     * local representation of grid vector quantities for rendering
     */
    private float [][] hostGridPos, hostGridVel, hostGridAccel;
    /**
     * local rep of grid scalars for rendering
     */
    private float[] hostGridMass;     
    /**
     * particle colors based on initial location
     */
    private int[][] hostPartClrAra, hostPartGreyAra;
        
    /**
     * 
     * @param _pa
     * @param _win
     * @param _simName
     * @param _currUIVals
     */
    public Base_MPMCudaSim(IGraphicsAppInterface _pa, Base_MPMSimWindow _win, String _simName, MPM_SimUpdateFromUIData _currUIVals) {
        super(_pa, _win, _simName, new float[] {0, 0, -9.8f}, _currUIVals);
        //redundant, but placed to specify that cuda kernel needs to be loaded
        ((MPM_CudaSimFlags) simFlags).setCudaDevInit(false);
        //Hold sim setup particle values
        partVals = new TreeMap<String, ArrayList<float[]>>();
        initPartArrays();
        //initialize cuda device pointers
        buildCudaDeviceConstructs();
        //set up grid and initialize sim with UI values and reset sim
        updateSimVals_FromUI(_currUIVals);
    }//Base_MPMCudaSim ctor
    
    /**
     * release all the device pointers
     */
    private final void _freeCUDADevPtrs() {
        //reset active ptrs
        //NOTE these destructors will not be called when this Java object goes out of scope or is otherwise destroyed/GC'ed.
        // It is the owning object's responsibility to clean these up when this object is about to be destroyed.
        JCuda.cudaFree(devPtrPartMass);          
        JCuda.cudaFree(devPtrPartVolume); 
        JCuda.cudaFree(devPtrGridMass); 
        for(int i=0;i<devPtrPartPos.length;++i) {
            JCuda.cudaFree(devPtrPartPos[i]);
            JCuda.cudaFree(devPtrPartVel[i]);
            JCuda.cudaFree(devPtrGridVel[i]); 
            JCuda.cudaFree(devPtrGridNewVel[i]);
            JCuda.cudaFree(devPtrGridForce[i]);
            
            for(int j=0;j<devPtrPartElasticF[0].length;++j) {
                JCuda.cudaFree(devPtrPartElasticF[i][j]);                    
                JCuda.cudaFree(devPtrPartPlasticF[i][j]);                    
            }
        }        
    }//_freeCUDADevPtrs()
    
    /**
     * Instance-specific reset code
     * @param rebuildSim
     */
    @Override
    protected final void resetSim_Indiv(SimResetProcess rebuildSim) {    
        _freeCUDADevPtrs();
        //sim start time - time from when sim object was first instanced
        //simStartTime = getCurTime();    
        
        msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "resetSim_Indiv","Start resetting sim");
        
        if (!((MPM_CudaSimFlags)simFlags).getCudaDevInit()) {
            //init cuda device and kernel file if not done already - only do 1 time
            msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "resetSim_Indiv","CUDA Module load/init");
            initCUDAModuleSetup();
        }

        //Set context
        JCudaDriver.cuCtxSetCurrent(context);    
        //initialize cude device and local/host constructs
        initHostAndDevPtrs();
        //rebuild cuda kernel configurations
        cudaSetup();          
        msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "resetSim_Indiv","Finished resetting sim");
   }//resetSim_Indiv  
    
    protected final void initHostAndDevPtrs() {
        //init ptrs to particle-based arrays - numparts and numPartsFloatSz need to be initialized by here
        initValues_Parts();        
        //init local reps of grid values and pointers to grid-based arrays
        initValues_Grid();
    }

    /**
     * Initialize environmental layout particle holders/arrays -  Reinitialize partVals map, called before sim environment is built
     */
    @Override
    protected final void initPartArrays() {
          partVals.clear();
        partVals.put("pos",new ArrayList<float[]>());
        partVals.put("vel",new ArrayList<float[]>());
        //Used to map colors of particles
        partVals.put("minMaxVals",new ArrayList<float[]>());
        //initialize min and max values at idx's 0 and 1
        partVals.get("minMaxVals").add(new float[] {100000.0f,100000.0f,100000.0f});
        partVals.get("minMaxVals").add(new float[] {-100000.0f,-100000.0f,-100000.0f}); 
    }//initPartArrays()
    
    /**
     * build initial layout for particles for this simulation
     * @param partVals [OUT] map of particle locs, initial velocities and min/max vals being constructed
     */
    @Override
    protected final void buildPartLayouts() {buildPartLayoutMap(partVals);}
    
    /**
     * Reinitialize existing sim - will resynthesize sample points but will not rederive locations of sim objs.
     * @param partVals
     */
    @Override
    protected final void reinitSimObjects() {reinitSimObjects(partVals);}

    /**
     * build initial layout for particles for this simulation
     * @param partVals [OUT] map of particle locs, initial velocities and min/max vals being constructed
     */
    protected abstract void buildPartLayoutMap(TreeMap<String, ArrayList<float[]>> partVals);    
    
    /**
     * Reinitialize existing sim - will resynthesize sample points but will not rederive locations of sim objs.
     * @param partVals
     */
    protected abstract void reinitSimObjects(TreeMap<String, ArrayList<float[]>> partVals);    
           
    /**
     * Update the CUDA-Specific simulator values whenever UI values change - called from base class, followed by resetSim call
     * @param upd
     */
    @Override
    protected final void updateSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd) {
        //float size of particle arrays, for malloc
        numPartsFloatSz = numParts * Sizeof.FLOAT; 
        //# cuda block/grid dims for particle functions
        numBlocksParticles = (numParts + numCUDAThreads -1)/numCUDAThreads;
        partGridDims = new int[] {numBlocksParticles, 1, 1};
        part4GridDims = new int[] {numBlocksParticles*4, 1, 1};

        //float size of grid arrays, for malloc       
        numGridFloatSz = ttlGridCount * Sizeof.FLOAT;
        //# cuda block/grid dims for grid based functions            
        numBlocksGrid = (ttlGridCount + numCUDAThreads -1)/numCUDAThreads;    
        gridGridDims = new int[] {numBlocksGrid, 1, 1};
        
        //update instancing sim values
        updateCudaSimVals_FromUI_Indiv(upd);        
    }//updateSimVals_FromUI_Indiv
    
    /**
     * Update instancing class variables based on potential UI changes
     * @param upd
     */
    protected abstract void updateCudaSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd);
    
    /**
     * run 1 time to load kernel and assign function pointers to functions
     */
    private void initCUDAModuleSetup() {
        // Enable exceptions and omit all subsequent error checks
        JCudaDriver.setExceptionsEnabled(true); 
        
        // Initialize the driver and create a context for the first device. (device 0)
        // build maps that hold values for cuda grid dim and shared mem size        
        cuInit(0);
        dev = new CUdevice();
        context = new CUcontext();
        cuDeviceGet(dev, 0);
        cuCtxCreate(context, 0, dev);
        
        // Load the ptx file.
        module = new CUmodule();
//        try {
//            compilePtxFile("src\\MPM.CudaKernels\\MPM_ABS_Sim.cu","MPM_ABS_Sim.ptx");
//        } catch (Exception e) {
//            System.out.println(e.getMessage());
//        }
        
        //Load MPM.CudaKernels module from ptxFileName
        cuModuleLoad(module, ptxFileName);        
        
        //Obtain a function pointer to each function in cuda kernel file and save it in cuFuncs
        cuFuncs = new TreeMap<String, CUfunction>();
        for (int i =0;i<CUFileFuncNames.length; ++i) {
            String key = CUFileFuncNames[i];            
            msgObj.dispInfoMessage("Base_MPMCudaSim("+simName+")", "initOnceCUDASetup","\tRegistering Kernel Function Key : " + key);
            CUfunction c = new CUfunction();            
            cuModuleGetFunction(c, module, key);
            cuFuncs.put(key,  c);
        }
    
        ((MPM_CudaSimFlags) simFlags).setCudaDevInit(true);
    }//initCUDAModuleSetup
    
    /**
     * allocate dev mem for all objects based on number of particles
     */
    @Override
    protected final void initValues_Parts() {
        float h_part_mass[] = new float[numParts];
        float h_part_eye[] = new float[numParts];
        //making class variables so can be rendered
        // for x,y,z values
        hostPartPos = new float[3][];
        hostPartVel = new float[3][];
        for(int i=0;i<hostPartPos.length;++i) {
            hostPartPos[i] = new float[numParts];
            hostPartVel[i] = new float[numParts];
        }       
        
        float[] minVals = partVals.get("minMaxVals").get(0);
        float[] maxVals = partVals.get("minMaxVals").get(1);       
        float[] posAra, velAra;
        
        hostPartClrAra = new int[numParts][3]; 
        hostPartGreyAra = new int[numParts][3]; 
        ArrayList<float[]> allPosAra = partVals.get("pos");
        ArrayList<float[]> allVelAra = partVals.get("vel");
        
        for(int i = 0; i < numParts; ++i){
            h_part_mass[i] = particleMass;
            posAra = allPosAra.get(i);
            velAra = allVelAra.get(i);            
            for(int j=0;j<hostPartPos.length;++j) {
                hostPartPos[j][i] = posAra[j];
                hostPartVel[j][i] = velAra[j];
            }            
            hostPartClrAra[i] = getClrValInt(posAra,minVals,maxVals);
            hostPartGreyAra[i] = getGreyValInt(posAra,minVals,maxVals);
            h_part_eye[i]=1.0f;
        }
        
        //Allocate memory for particle-related cuda pointers
        cuMemAlloc(devPtrPartVolume, numPartsFloatSz); 
        cuMemsetD32(devPtrPartVolume, 0, numParts);    //part_vol is a calculated quantity
        cuMemAlloc(devPtrPartMass, numPartsFloatSz);  
        cuMemcpyHtoD(devPtrPartMass, Pointer.to(h_part_mass), numPartsFloatSz);
        //allocate for all particle constructs
        for(int i=0;i<devPtrPartPos.length;++i) {
               cuMemAlloc(devPtrPartPos[i], numPartsFloatSz); 
               cuMemAlloc(devPtrPartVel[i], numPartsFloatSz);
            cuMemcpyHtoD(devPtrPartPos[i], Pointer.to(hostPartPos[i]), numPartsFloatSz);
            cuMemcpyHtoD(devPtrPartVel[i], Pointer.to(hostPartVel[i]), numPartsFloatSz);
            for(int j=0;j<devPtrPartElasticF[0].length;++j) {        
                //build identity matrices for this
                cuMemAlloc(devPtrPartElasticF[i][j], numPartsFloatSz);
                cuMemAlloc(devPtrPartPlasticF[i][j], numPartsFloatSz); 
                if(i==j) {
                    cuMemcpyHtoD(devPtrPartElasticF[i][j], Pointer.to(h_part_eye), numPartsFloatSz);                    
                    cuMemcpyHtoD(devPtrPartPlasticF[i][j], Pointer.to(h_part_eye), numPartsFloatSz);                            
                } else {
                    //set nonDiagonal values to zero
                    cuMemsetD32(devPtrPartElasticF[i][j], 0, numParts);                    
                    cuMemsetD32(devPtrPartPlasticF[i][j], 0, numParts);                 
                }            
            }                    
        }
        
    }//initValues_Parts
    
    /**
     * allocate dev mem for all objects based on number of grid cells
     */
    @Override
    protected final void initValues_Grid() {
        // for x,y,z values
        hostGridPos = new float[3][];
        hostGridVel = new float[3][];
        hostGridAccel = new float[3][];
        for(int i=0;i<hostGridPos.length;++i) {
            hostGridPos[i] = new float[ttlGridCount];
            hostGridVel[i] = new float[ttlGridCount];
               hostGridAccel[i] = new float[ttlGridCount];
        }
        
        hostGridMass = new float[ttlGridCount];        
        //build grid locations
        int gridDim=0;
        for(int i=0;i<gridSideCount;++i) {
            float xPos = (i+.5f)*cellSize;
            for(int j=0;j<gridSideCount;++j) {        
                float yPos = (j+.5f)*cellSize;
                for(int k=0;k<gridSideCount;++k) {
                    float zPos = (k+.5f)*cellSize;
                    //weird orientation stuffs
                    //TODO fix grid orientation
                    hostGridPos[0][gridDim] = zPos;
                    hostGridPos[1][gridDim] = yPos;
                    hostGridPos[2][gridDim] = xPos;                    
                    ++gridDim;
                }
            }
        }            
        
        //Allocate memory and initialize values for cuda grid pointers
        cuMemAlloc(devPtrGridMass, numGridFloatSz); 
        cuMemsetD32(devPtrGridMass, 0, ttlGridCount);
        for(int i=0;i<devPtrGridVel.length;++i) {
            cuMemAlloc(devPtrGridVel[i], numGridFloatSz);  
            cuMemAlloc(devPtrGridNewVel[i], numGridFloatSz);
            cuMemAlloc(devPtrGridForce[i], numGridFloatSz);
            
            cuMemsetD32(devPtrGridVel[i], 0, ttlGridCount);  
            cuMemsetD32(devPtrGridNewVel[i], 0, ttlGridCount);
            cuMemsetD32(devPtrGridForce[i], 0, ttlGridCount);               
        }
    }//initValues_Grids
    
    /**
     * Only performed from ctor
     */
    private final void buildCudaDeviceConstructs() {
        // particle scalars
        devPtrPartMass = new CUdeviceptr();           devPtrPartVolume = new CUdeviceptr();
        // grid scalar
        devPtrGridMass = new CUdeviceptr();
        // particle vectors - each idx is x,y,z
        devPtrPartPos = new CUdeviceptr[3];            devPtrPartVel = new CUdeviceptr[3];
        // grid vectors - each idx is x,y,z
        devPtrGridVel = new CUdeviceptr[3];            devPtrGridNewVel = new CUdeviceptr[3];            devPtrGridForce = new CUdeviceptr[3];
        // deformation matrices 3x3
        devPtrPartElasticF = new CUdeviceptr[3][];
        devPtrPartPlasticF = new CUdeviceptr[3][];
        for(int i=0;i<3;++i) {            
            devPtrPartPos[i] = new CUdeviceptr(); 
            devPtrPartVel[i] = new CUdeviceptr();        
        
            devPtrGridVel[i] = new CUdeviceptr();            
            devPtrGridNewVel[i] = new CUdeviceptr();         
            devPtrGridForce[i] = new CUdeviceptr();
            devPtrPartElasticF[i] = new CUdeviceptr[3];
            devPtrPartPlasticF[i] = new CUdeviceptr[3];            
            for(int j=0;j<3;++j) {
                devPtrPartElasticF[i][j] = new CUdeviceptr();
                devPtrPartPlasticF[i][j] = new CUdeviceptr();
            }
        }        
    }//buildCudaDeviceConstructs    
    
    /**
     * create a sphere with given center, with passed # of particles -0 returns array of [start IDX,end IDX] within posVals array for this sphere
     * @param partVals map holding all appropriate particle values
     * @param ballRad desired radius of sphere
     * @param numParts
     * @param initVel
     * @param ctr
     * @return array with initial, and final, size of position values == array of [start IDX, end IDX] within posVals array for the points making up this sphere.
     */           
    protected final int[] createSphere(TreeMap<String, ArrayList<float[]>> partVals, float ballRad, int numParts, myPointf ctr) {            
        float[] minVals = partVals.get("minMaxVals").get(0);
        float[] maxVals = partVals.get("minMaxVals").get(1); 
        ArrayList<float[]> posAra = partVals.get("pos");
        int[] returnIdxs = new int[2];
        //start at beginning of current posMap
        returnIdxs[0] = posAra.size();
        for (int i=0;i<numParts;++i) {
            //float[] posVals = getRandPosInSphereAra(ballRad, ctr); 
            float[] posVals = MyMathUtils.getRandPosInSphere(ballRad, ctr).asArray(); 
            //find min/max values for all sphere particles
            for (int v = 0; v < 3; ++v) {
                if (posVals[v] < minVals[v]) {                    minVals[v] = posVals[v];                } 
                else if (posVals[v] > maxVals[v]) {                maxVals[v] = posVals[v];                }    
            }
            posAra.add(posVals);
        }
        msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "createSphere",
                "Created a sphere of radius " + ballRad + " with "+numParts+" particles, centered at [" +ctr.toStrBrf() + "].");
        //Ending at final size of posMap
        returnIdxs[1] = posAra.size();
        return returnIdxs;
    }//createSphere
    
    protected final void setPartInitVelocities(TreeMap<String, ArrayList<float[]>> partVals, int stIdx, int endIdx, float[] initVel) {
        ArrayList<float[]> velMap = partVals.get("vel");
        for (int i=stIdx; i<endIdx;++i) {            velMap.add(initVel);        }
    }//setPartInitVelocities

    
    /**
     * launch the kernel specified by the string key
     * @param key string key of kernel to launch and kernel function grid, block and mem dims to access
     */
    protected void launchKernel(String key) {
        //Set context properly before launching kernels
        JCudaDriver.cuCtxSetCurrent(context);
        int[][] kernelDims = funcGridDimAndMemSize.get(key);
        cuLaunchKernel(cuFuncs.get(key),                                         // Kernel function
                kernelDims[0][0],kernelDims[0][1],kernelDims[0][2],               // Grid XYZ dimensions
                kernelDims[1][0],kernelDims[1][1],kernelDims[1][2],                // Block XYZ dimensions
                kernelDims[2][0],                                                // Shared memory size
                null,                                                             // Stream 
                kernelParams.get(key), null);                                    // Kernel- and extra parameters
        //Allow kernel to complete
        cuCtxSynchronize();
    }
    /**
     * Build grid, block and shared mem dims into map holding these for kernel launch
     * @param key
     * @param gridDims
     * @param blockThdDims
     * @param sharedMemSize
     */
    protected void putFuncGridMemSize(String key, int[] gridDims, int[] blockThdDims, int[] sharedMemSize) {
        funcGridDimAndMemSize.put(key, new int[][] {gridDims, blockThdDims, sharedMemSize});
    }
    
    /**
     * Set up all essential cuda kernels and launch them for initial pass
     */
    private void cudaSetup() {    
        msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Start CUDA Init.");
         //Re initialize maps of parameters and functions
        kernelParams = new TreeMap<String, Pointer>();
        funcGridDimAndMemSize = new HashMap<String, int[][]>();
        Pointer numPartsPtr = Pointer.to(new int[] {numParts});
        Pointer partMassPtr = Pointer.to(devPtrPartMass);
        Pointer partVolPtr = Pointer.to(devPtrPartVolume);
        Pointer gridMassPtr = Pointer.to(devPtrGridMass);
        Pointer[] partPosPtrAra = new Pointer[]{Pointer.to(devPtrPartPos[0]),Pointer.to(devPtrPartPos[1]),Pointer.to(devPtrPartPos[2])};
        Pointer[] partVelPtrAra = new Pointer[]{Pointer.to(devPtrPartVel[0]),Pointer.to(devPtrPartVel[1]),Pointer.to(devPtrPartVel[2])};
        
        Pointer cellSizePtr = Pointer.to(new float[] {cellSize});
        Pointer numGridSidePtr = Pointer.to(new int[]{gridSideCount});
        Pointer ttlGridCountPtr = Pointer.to(new int[] {ttlGridCount});
        Pointer[] gridVelPtrAra = new Pointer[]{Pointer.to(devPtrGridVel[0]),Pointer.to(devPtrGridVel[1]),Pointer.to(devPtrGridVel[2])};
        Pointer[] gridNewVelPtrAra = new Pointer[]{Pointer.to(devPtrGridNewVel[0]),Pointer.to(devPtrGridNewVel[1]),Pointer.to(devPtrGridNewVel[2])};
        Pointer[] gridFrcPtrAra = new Pointer[]{Pointer.to(devPtrGridForce[0]),Pointer.to(devPtrGridForce[1]),Pointer.to(devPtrGridForce[2])};
        float[] gravity = getGravity().asArray();
        Pointer[] gravPtrAra = new Pointer[] {Pointer.to(new float[] {gravity[0]}),Pointer.to(new float[] {gravity[1]}),Pointer.to(new float[] {gravity[2]})};
        
        // Scalar global quantities
        Pointer minSimBndsPtr = Pointer.to(new float[] {minSimBnds});
        Pointer maxSimBndsPtr = Pointer.to(new float[] {maxSimBnds});
        Pointer deltaTPtr = Pointer.to(new float[] {deltaT});
        Pointer wallFricPtr = Pointer.to(new float[] {getWallFric()});        
        // Compute the grid forces by projecting the particle quantities to the grid using appropriate weighting,
        // 
        kernelParams.put("projectToGridandComputeForces",Pointer.to(
                numPartsPtr, numGridSidePtr, cellSizePtr, minSimBndsPtr,
                Pointer.to(mat.getLambda0Ptr()), Pointer.to(mat.getMu0Ptr()), Pointer.to(mat.getHardeningCoeffPtr()),
                partMassPtr, partVolPtr,
                partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
                partVelPtrAra[0], partVelPtrAra[1], partVelPtrAra[2],
                //elastic matrix
                Pointer.to(devPtrPartElasticF[0][0]), Pointer.to(devPtrPartElasticF[0][1]), Pointer.to(devPtrPartElasticF[0][2]),
                Pointer.to(devPtrPartElasticF[1][0]), Pointer.to(devPtrPartElasticF[1][1]), Pointer.to(devPtrPartElasticF[1][2]),
                Pointer.to(devPtrPartElasticF[2][0]), Pointer.to(devPtrPartElasticF[2][1]), Pointer.to(devPtrPartElasticF[2][2]),
                //plastic matrix
                Pointer.to(devPtrPartPlasticF[0][0]), Pointer.to(devPtrPartPlasticF[0][1]), Pointer.to(devPtrPartPlasticF[0][2]),
                Pointer.to(devPtrPartPlasticF[1][0]), Pointer.to(devPtrPartPlasticF[1][1]), Pointer.to(devPtrPartPlasticF[1][2]),
                Pointer.to(devPtrPartPlasticF[2][0]), Pointer.to(devPtrPartPlasticF[2][1]), Pointer.to(devPtrPartPlasticF[2][2]),

                gridMassPtr,
                gridVelPtrAra[0], gridVelPtrAra[1], gridVelPtrAra[2],
                gridFrcPtrAra[0], gridFrcPtrAra[1], gridFrcPtrAra[2]));
        putFuncGridMemSize("projectToGridandComputeForces", part4GridDims, blkThdDims, shrdMemSize);
        // initial grid projection of particles' mass using position only, for calculation of particle volume 
        kernelParams.put("projectToGridInit", Pointer.to(
                numPartsPtr, numGridSidePtr, cellSizePtr, minSimBndsPtr,
                partMassPtr,
                partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
                //Output
                gridMassPtr));
        putFuncGridMemSize("projectToGridInit", partGridDims, blkThdDims, shrdMemSize);
        // compute the approx volume of each particle - first step only
        // Determine the grid cells' density rho using each's mass over its volume, then project
        // back to each particle, weighting appropriately. The particle's volume will then be the mass*rho  
        // 
        kernelParams.put("computeVol", Pointer.to(
                numPartsPtr, numGridSidePtr, cellSizePtr, minSimBndsPtr,
                partMassPtr, partVolPtr,
                partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
                gridMassPtr));
        putFuncGridMemSize("computeVol", partGridDims, blkThdDims, shrdMemSize);
        
        kernelParams.put("updPartVelocities",Pointer.to(
                numPartsPtr, numGridSidePtr, cellSizePtr, minSimBndsPtr,
                Pointer.to(mat.getAlphaPicFlipPtr()),
                partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
                partVelPtrAra[0], partVelPtrAra[1], partVelPtrAra[2],
                gridVelPtrAra[0], gridVelPtrAra[1], gridVelPtrAra[2],        
                gridNewVelPtrAra[0], gridNewVelPtrAra[1], gridNewVelPtrAra[2]));  
        putFuncGridMemSize("updPartVelocities", part4GridDims, blkThdDims, partVelShrdMemSize);
        
        kernelParams.put("compGridVelocities", Pointer.to(
                ttlGridCountPtr, gravPtrAra[0], gravPtrAra[1], gravPtrAra[2],deltaTPtr,
                gridMassPtr,        
                gridVelPtrAra[0], gridVelPtrAra[1], gridVelPtrAra[2],            
                gridNewVelPtrAra[0], gridNewVelPtrAra[1], gridNewVelPtrAra[2],
                gridFrcPtrAra[0], gridFrcPtrAra[1], gridFrcPtrAra[2]));
        putFuncGridMemSize("compGridVelocities", gridGridDims, blkThdDims, shrdMemSize);
        
        kernelParams.put("clearGrid", Pointer.to(
                ttlGridCountPtr, gridMassPtr,        
                gridVelPtrAra[0], gridVelPtrAra[1], gridVelPtrAra[2],            
                gridNewVelPtrAra[0], gridNewVelPtrAra[1], gridNewVelPtrAra[2],
                gridFrcPtrAra[0], gridFrcPtrAra[1], gridFrcPtrAra[2]));
        putFuncGridMemSize("clearGrid", gridGridDims, blkThdDims, shrdMemSize);
        
        //Only currently supports wall collisions
        kernelParams.put("gridCollisions", Pointer.to(
                ttlGridCountPtr, numGridSidePtr, cellSizePtr, 
                minSimBndsPtr, maxSimBndsPtr,
                wallFricPtr, deltaTPtr, gridMassPtr,
                gridNewVelPtrAra[0], gridNewVelPtrAra[1], gridNewVelPtrAra[2]));
        putFuncGridMemSize("gridCollisions", gridGridDims, blkThdDims, shrdMemSize);
        // Update the elastic and plastic deformation gradients for each particle by 
        // finding the next timestep's full deformation gradient.
        // F_t+1 = Fehat_t+1 * Fphat_t+1 = (I + delT*Grad(vp_t+1))Fe_t * Fp_t
        //  - first attribute all the deformation changes to the last timestep's elastic defgrad
        //  - next, derive Fe_t+1 from Fehat_t+1 by taking SVD of Fehat_t+1 = UZhatV', clamping the singular 
        // values in Zhat to [1-critComp, 1+critStretch] -> Z, and then reconstructing Fe_t+1 = UZV'
        // lastly, derive Fp_t+1 = inv(Fe_t+1)*F_t+1 = V*inv(Z)*U'*F_t+1
        kernelParams.put("updDeformationGradient", Pointer.to(
                numPartsPtr, numGridSidePtr,deltaTPtr, cellSizePtr, minSimBndsPtr,
                // these are used to clamp the singular values of the updated elastic deformation gradient
                Pointer.to(mat.getCriticalCompressionPtr()), Pointer.to(mat.getCriticalStretchPtr()),
                partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
                //elastic matrix
                Pointer.to(devPtrPartElasticF[0][0]), Pointer.to(devPtrPartElasticF[0][1]), Pointer.to(devPtrPartElasticF[0][2]),
                Pointer.to(devPtrPartElasticF[1][0]), Pointer.to(devPtrPartElasticF[1][1]), Pointer.to(devPtrPartElasticF[1][2]),
                Pointer.to(devPtrPartElasticF[2][0]), Pointer.to(devPtrPartElasticF[2][1]), Pointer.to(devPtrPartElasticF[2][2]),
                //plastic matrix
                Pointer.to(devPtrPartPlasticF[0][0]), Pointer.to(devPtrPartPlasticF[0][1]), Pointer.to(devPtrPartPlasticF[0][2]),
                Pointer.to(devPtrPartPlasticF[1][0]), Pointer.to(devPtrPartPlasticF[1][1]), Pointer.to(devPtrPartPlasticF[1][2]),
                Pointer.to(devPtrPartPlasticF[2][0]), Pointer.to(devPtrPartPlasticF[2][1]), Pointer.to(devPtrPartPlasticF[2][2]),
                //new grid velocities for velocity gradient
                gridNewVelPtrAra[0], gridNewVelPtrAra[1], gridNewVelPtrAra[2]));     
        putFuncGridMemSize("updDeformationGradient", partGridDims, blkThdDims, shrdMemSize);
        
        //Only addresses wall collisions
        kernelParams.put("partCollAndUpdPos", Pointer.to(
                numPartsPtr, minSimBndsPtr, maxSimBndsPtr,
                wallFricPtr, deltaTPtr,
                partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
                partVelPtrAra[0], partVelPtrAra[1], partVelPtrAra[2]            
                ));   
        putFuncGridMemSize("partCollAndUpdPos", partGridDims, blkThdDims, shrdMemSize);
   
        msgObj.dispInfoMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Finished CUDA Init | Launch first MPM Pass.");
         //launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {
            msgObj.dispInfoMessage("Base_MPMCudaSim("+simName+")", "cudaSetup", "Launching kernel "+j);
            launchKernel(initStepFuncKeys[j]);
            
        }

        simFlags.setSimIsBuilt(true);
        msgObj.dispInfoMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Finished first MPM Pass.");
    }//cudaSetup

    /**
     * Fully reads the given InputStream and returns it as a byte array
     *  
     * @param inputStream The input stream to read
     * @return The byte array containing the data from the input stream
     * @throws IOException If an I/O error occurs
     */
    private static byte[] toByteArray(InputStream inputStream) throws IOException{
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        byte buffer[] = new byte[8192];
        while (true){
            int read = inputStream.read(buffer);
            if (read == -1){       break;     }
            bytes.write(buffer, 0, read);
        }
        return bytes.toByteArray();
    }    
    
    /**
     * Compiles Ptx file from file in passed file name -> cuFileName needs to have format "xxxxx.cu"
     * @param krnFileName MPM.CudaKernels source code kernel file name
     * @param ptxFileName File name for output MPM.CudaKernels ptx file.
     * @throws IOException
     */
    public final void compilePtxFile(String krnFileName, String ptxFileName) throws IOException {
        File cuFile = new File(krnFileName);
        if (!cuFile.exists()) {
            throw new IOException("Kernel file not found: " + krnFileName);
        }
        String modelString = "-m" + System.getProperty("sun.arch.data.model");
        //build compilation command
        String command = "nvcc " + modelString + " -ptx " + cuFile.getPath() + " -o " + ptxFileName;
        //execute compilation
        msgObj.dispInfoMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","Executing\n" + command);
        Process process = Runtime.getRuntime().exec(command);

        String errorMessage = new String(toByteArray(process.getErrorStream())), 
                outputMessage = new String(toByteArray(process.getInputStream()));
        int exitValue = 0;
        try {exitValue = process.waitFor();} 
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while waiting for nvcc output", e);
        }

        if (exitValue != 0) {
            msgObj.dispErrorMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","nvcc process error : exitValue : " + exitValue);
            msgObj.dispErrorMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","errorMessage :\n" + errorMessage);
            msgObj.dispErrorMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","outputMessage :\n" + outputMessage);
            throw new IOException("Could not create .ptx file: " + errorMessage);
        }
        msgObj.dispInfoMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","Finished compiling PTX file : "+ ptxFileName);
    }//compilePtxFile
    
    /**
     * Instance-specific per-sim cycle simulation execution code
     * @param modAmtMillis
     * @return
     */
    @Override
    protected final boolean simMe_Indiv(float modAmtMillis) {
        //for every sim step, launch each kernel by key specified in simStepFuncKeys
        for (int j=0;j<simStepFuncKeys.length;++j) {                launchKernel(simStepFuncKeys[j]);        }
        return false;
    }
    
    /**
     * Instance-specific post-sim cyle code
     * @param modAmtMillis
     * @return
     */
    @Override
    protected final boolean simMePost_Indiv(float modAmtMillis) {
        msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "simMePost_Indiv","Start copy relevant device buffers to host.");
        //copy from device data to host particle position or velocity arrays
         if(simFlags.getShowParticles() || simFlags.getShowPartVels()) {
             for(int i=0;i<hostPartPos.length;++i) {            cuMemcpyDtoH(Pointer.to(hostPartPos[i]), devPtrPartPos[i], numPartsFloatSz);}
         }             
         if(simFlags.getShowPartVels()) {
             for(int i=0;i<hostPartVel.length;++i) {         cuMemcpyDtoH(Pointer.to(hostPartVel[i]), devPtrPartVel[i], numPartsFloatSz);}
         }
         //copy from device data to host grid velocity, accel or mass arrays
        if(simFlags.getShowGridVel()) {
             for(int i=0;i<hostGridVel.length;++i) {         cuMemcpyDtoH(Pointer.to(hostGridVel[i]), devPtrGridNewVel[i], numGridFloatSz);}
        }
        if(simFlags.getShowGridAccel()) {
             for(int i=0;i<hostGridAccel.length;++i) {         cuMemcpyDtoH(Pointer.to(hostGridAccel[i]), devPtrGridForce[i], numGridFloatSz);}        
        }        
        if(simFlags.getShowGridMass()) {                    cuMemcpyDtoH(Pointer.to(hostGridMass), devPtrGridMass, numGridFloatSz);}
        msgObj.dispDebugMessage("Base_MPMCudaSim("+simName+")", "simMePost_Indiv","End copy relevant device buffers to host.");
        return false;
    }

    private final int[] gridVelClr = new int[] {250,0,0}, gridAccelClr = new int[] {0,0,250}, gridMassClr = new int[] {0,150,50};
    //private final int[] whitePoints = new int[] {255,255,255}; 

    /**
     * Draw instance class particles
     * @param animTimeMod
     * @param showLocColors
     */
    @Override
    protected final void _drawParts(float animTimeMod, boolean showLocColors, boolean isGlblAppDebug) {
        ri.pushMatState();
        ri.drawPointCloudWithColors(hostPartPos[0].length, drawPointIncr, (showLocColors ?  hostPartClrAra : hostPartGreyAra), hostPartPos[0], hostPartPos[1], hostPartPos[2]);
        ri.popMatState();
    }//_drawParts
    
    /**
     * Draw instance class particle velocities
     * @param animTimeMod
     * @param minMag minimum magnitude per axis to draw vector
     * @param pincr
     */
    @Override
    protected final void _drawPartVel(float animTimeMod, float vecScale, int pincr, boolean isGlblAppDebug) {
        ri.pushMatState();        
        for(int i=0;i<=hostPartVel[0].length-pincr;i+=pincr) {                    
            if(        (Math.abs(hostPartVel[0][i]) > vecScale) || 
                    (Math.abs(hostPartVel[1][i]) > vecScale) || 
                    (Math.abs(hostPartVel[2][i]) > vecScale)) {
                ri.pushMatState();
                ri.setStroke(hostPartClrAra[i], 255);
                ri.translate(hostPartPos[0][i], hostPartPos[1][i], hostPartPos[2][i]);
                ri.drawLine(0,0,0, vecLengthScale*hostPartVel[0][i],vecLengthScale*hostPartVel[1][i],vecLengthScale*hostPartVel[2][i]);
                ri.popMatState();
            }
        }            
        ri.popMatState();        
    }//_drawPartVel
    
    /**
     * Draw instance class grid velocities - use _drawGridVec method
     * @param animTimeMod
     * @param minMag minimum magnitude per axis to draw vector
     */
    @Override
    protected final void _drawGridVel(float animTimeMod, float minMag, boolean isGlblAppDebug) {        
        _drawGridVec(gridVelClr, hostGridVel, hostGridPos, minMag, isGlblAppDebug);
    }

    /**
     * Draw instance class grid accelerations - use _drawGridVec method
     * @param animTimeMod
     * @param minMag minimum magnitude per axis to draw vector
     */
    @Override
    protected final void _drawGridAccel(float animTimeMod, float minMag, boolean isGlblAppDebug) {    
        _drawGridVec(gridAccelClr, hostGridAccel, hostGridPos, minMag, isGlblAppDebug);
    }
    
    /**
     * Draw instance class grid masses - use _drawGridScalar method
     * @param animTimeMod
     * @param minMag minimum magnitude to draw scalar mass
     */
    @Override
    protected final void _drawGridMass(float animTimeMod, float minMag, boolean isGlblAppDebug) {        
        _drawGridScalar(gridMassClr, hostGridMass, hostGridPos, minMag, isGlblAppDebug);
    }    
    
    /**
     * Draw any colliders if they exist
     * @param animTimeMod
     */
    @Override
    protected final void _drawColliders(float animTimeMod, boolean isGlblAppDebug) {
        ri.pushMatState();    
        drawColliders_Indiv(animTimeMod);
        ri.popMatState();
    }
    /**
     * draw internal-to-sim colliders, if they exist
     * @param animTimeMod
     */
    protected abstract void drawColliders_Indiv(float animTimeMod);
    
}//class Base_MPMCudaSim 




