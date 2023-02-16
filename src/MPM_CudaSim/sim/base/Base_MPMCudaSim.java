package MPM_CudaSim.sim.base;

import static jcuda.driver.JCudaDriver.*;

import java.io.*;
import java.util.*;

import MPM_SimMain.sim.SimResetProcess;
import MPM_SimMain.sim.Base_MPMSim;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.windowUI.base.Base_DispWindow;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import jcuda.*;
import jcuda.driver.*;
import jcuda.runtime.JCuda;

/**
 * abstract class describing a simulation world. Called by sim window to 
 * execute simulation and render results. Instancing classes can hold 
 * different configurations/initialization setups
 */
public abstract class Base_MPMCudaSim extends Base_MPMSim{	
	
	/**
	 * CUDA kernel file name
	 */
	private String ptxFileName = "MPM_CUDA_Sim_New.ptx";	
	////////////////////////////////////////////////////
    // Maps holding CUDA function pointers and parameter pointers
	// This facilitates CUDA calcs and access to appropriately configured CUDA args
	protected TreeMap<String, Pointer> kernelParams;
	protected TreeMap<String, CUfunction> cuFuncs;
	/**
	 * Lists of names of kernel functions to perform for MPM algorithm.   
	 */
	protected String[] CUFileFuncNames = new String[] {
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
	protected String[] initStepFuncKeys = new String[] {
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
	protected String[] simStepFuncKeys = new String[] {
			"clearGrid",
			"projectToGridandComputeForces",
			"compGridVelocities",
			"gridCollisions",
			"updDeformationGradient",
			"updPartVelocities",
			"partCollAndUpdPos"}; 
   
	protected HashMap<String, int[][]> funcGridDimAndMemSize;
	
	////////////////////////////////////////////////////
    // CUDA references
    protected CUdevice dev;
    protected CUcontext context;
    protected CUmodule module;
    protected CUgraphicsResource pCudaResource;
	// CUDA Device ptr constructions
	protected CUdeviceptr partMass, partVolume;
	protected CUdeviceptr[] partPos, partVel; 
    protected CUdeviceptr[] gridVel, gridNewVel, gridForce;
	protected CUdeviceptr[][] partElasticF, partPlasticF;    
    protected CUdeviceptr gridMass;
    // CUDA calc helper variables

    /**
     * # of cuda blocks to use for particles
     */
    protected int numBlocksParticles;
    /**
     * # of cuda blocks to use for grid
     */
    protected int numBlocksGrid;
    /**
     * # parts * size of float & num grid cells * size float 
     */
	protected long numPartsFloatSz, numGridFloatSz;
	/**
	 * # of cuda threads
	 */
    protected final int numCUDAThreads = 512;
    protected final int[] blkThdDims = new int[] {numCUDAThreads, 1, 1};
    protected int[] partGridDims;
    protected int[] part4GridDims;
    protected int[] gridGridDims;
    protected int[] shrdMemSize = new int[] {0};
    protected int[] partVelShrdMemSize = new int[] {numCUDAThreads*6*Sizeof.FLOAT};
    
    /**
     * Raw initial particle values
     */
    TreeMap<String, ArrayList<float[]>> partVals;
    
    ////////////////////////////////////////////////////
    //representations for rendering
    /**
     * local representation of particle vector quantities for rendering
     */
    protected float[][] hostPartPos, hostPartVel;
    /**
     * local representation of grid vector quantities for rendering
     */
    protected float [][] hostGridPos, hostGridVel, hostGridAccel;
    /**
     * local rep of grid scalars for rendering
     */
	protected float[] hostGridMass;     
	/**
	 * particle colors based on initial location
	 */
    protected int[][] hostPartClrAra, hostPartGreyAra;
        
	/**
	 * 
	 * @param _pa
	 * @param _win
	 * @param _simName
	 * @param _currUIVals
	 */
	public Base_MPMCudaSim(IRenderInterface _pa, Base_MPMSimWindow _win, String _simName, MPM_SimUpdateFromUIData _currUIVals) {
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
	}//MPM_ABS_Sim
	
	/**
	 * Instance-specific reset code
	 * @param rebuildSim
	 */
    @Override
    protected final void resetSim_Indiv(SimResetProcess rebuildSim) {	
		//reset active ptrs
		JCuda.cudaFree(partMass);		  
        JCuda.cudaFree(partVolume); 
        JCuda.cudaFree(gridMass); 
        for(int i=0;i<partPos.length;++i) {
        	JCuda.cudaFree(partPos[i]);
        	JCuda.cudaFree(partVel[i]);
	        JCuda.cudaFree(gridVel[i]); 
	        JCuda.cudaFree(gridNewVel[i]);
	        JCuda.cudaFree(gridForce[i]);
	        
	        for(int j=0;j<partElasticF[0].length;++j) {
	        	JCuda.cudaFree(partElasticF[i][j]);                    
	        	JCuda.cudaFree(partPlasticF[i][j]);		        	
	        }
        }
		//sim start time - time from when sim object was first instanced
		//simStartTime = getCurTime();	
		
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "resetSim_Indiv","Start resetting sim");
		
		if (!((MPM_CudaSimFlags)simFlags).getCudaDevInit()) {
            //init cuda device and kernel file if not done already - only do 1 time
			win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "resetSim_Indiv","CUDA Module load/init");
            initCUDAModuleSetup();
        }

		//Set context
		JCudaDriver.cuCtxSetCurrent(context);    
        //init ptrs to particle-based arrays - numparts and numPartsFloatSz need to be initialized by here
		initValues_Parts();		
        //init local reps of grid values and pointers to grid-based arrays
		initValues_Grids();
		
		//rebuild cuda kernel configurations
		cudaSetup();  	    
    }//resetSim_Indiv    

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
        //initialize min and max values
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
		//# cuda grid dims for particle functions
		numBlocksParticles = (numParts + numCUDAThreads -1)/numCUDAThreads;
	    partGridDims = new int[] {numBlocksParticles, 1, 1};
	    part4GridDims = new int[] {numBlocksParticles*4, 1, 1};

	    //float size of grid arrays, for malloc       
        numGridFloatSz = ttlGridCount * Sizeof.FLOAT;
        //# cuda grid dims for grid based functions			
		numBlocksGrid = (ttlGridCount + numCUDAThreads -1)/numCUDAThreads;	
	    gridGridDims = new int[] {numBlocksGrid, 1, 1};
		
		//update instancing sim values
		updateCudaSimVals_FromUI_Indiv(upd);		
	}//updateSimVals_FromUI
	
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
//		try {
//			compilePtxFile("src\\Cuda\\MPM_ABS_Sim.cu","MPM_ABS_Sim.ptx");
//			ptxFileName = "MPM_ABS_Sim.ptx";
//			//ptxFileName = MPM_ABS_Sim.preparePtxFile("src\\Cuda\\MPM_ABS_Sim.cu");
//		} catch (Exception e) {
//			System.out.println(e.getMessage());
//		}
		
		//Load Cuda module from ptxFileName
		cuModuleLoad(module, ptxFileName);    	
		
        //Obtain a function pointer to each function in cuda kernel file and save it in cuFuncs
		cuFuncs = new TreeMap<String, CUfunction>();
		for (int i =0;i<CUFileFuncNames.length; ++i) {
			String key = CUFileFuncNames[i];			
			win.getMsgObj().dispInfoMessage("Base_MPMCudaSim("+simName+")", "initOnceCUDASetup","\tRegistering Kernel Function Key : " + key);
			CUfunction c = new CUfunction();			
			cuModuleGetFunction(c, module, key);
			cuFuncs.put(key,  c);
		}
    
        ((MPM_CudaSimFlags) simFlags).setCudaDevInit(true);
	}//loadModuleAndSetFuncPtrs
	
	/**
	 * allocate dev mem for all objects based on number of particles
	 */
	@Override
	protected final void initValues_Parts() {
        float h_part_mass[] = new float[numParts];
        float h_part_eye[] = new float[numParts];
        //making class variables so can be rendered
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
        cuMemAlloc(partVolume, numPartsFloatSz); 
        cuMemsetD32(partVolume, 0, numParts);	//part_vol is a calculated quantity
        cuMemAlloc(partMass, numPartsFloatSz);  
        cuMemcpyHtoD(partMass, Pointer.to(h_part_mass), numPartsFloatSz);
        
        for(int i=0;i<partPos.length;++i) {
           	cuMemAlloc(partPos[i], numPartsFloatSz); 
           	cuMemAlloc(partVel[i], numPartsFloatSz);
            cuMemcpyHtoD(partPos[i], Pointer.to(hostPartPos[i]), numPartsFloatSz);
            cuMemcpyHtoD(partVel[i], Pointer.to(hostPartVel[i]), numPartsFloatSz);
            for(int j=0;j<partElasticF[0].length;++j) {		//build identity mats for this
                cuMemAlloc(partElasticF[i][j], numPartsFloatSz);
                cuMemAlloc(partPlasticF[i][j], numPartsFloatSz); 
                if(i==j) {
                    cuMemcpyHtoD(partElasticF[i][j],	Pointer.to(h_part_eye), numPartsFloatSz);                    
                    cuMemcpyHtoD(partPlasticF[i][j],	Pointer.to(h_part_eye), numPartsFloatSz);             	       	
                } else {
                    cuMemsetD32(partElasticF[i][j], 0, numParts);                    
                    cuMemsetD32(partPlasticF[i][j], 0, numParts);             	
                }        	
            }    	        	
        }
        
    }//initCUDAMemPtrs_Parts
	
	/**
	 * allocate dev mem for all objects based on number of grid cells
	 */
	@Override
	protected final void initValues_Grids() {
		hostGridPos = new float[3][];
		hostGridVel = new float[3][];
		hostGridAccel = new float[3][];
		for(int i=0;i<hostPartPos.length;++i) {
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
        cuMemAlloc(gridMass, numGridFloatSz); 
		cuMemsetD32(gridMass, 0, ttlGridCount);
		for(int i=0;i<gridVel.length;++i) {
            cuMemAlloc(gridVel[i], numGridFloatSz);  
            cuMemAlloc(gridNewVel[i], numGridFloatSz);
            cuMemAlloc(gridForce[i], numGridFloatSz);
            
            cuMemsetD32(gridVel[i], 0, ttlGridCount);  
            cuMemsetD32(gridNewVel[i], 0, ttlGridCount);
            cuMemsetD32(gridForce[i], 0, ttlGridCount);           	
        }
	}//initValues_Grids
	
	/**
	 * Only performed from ctor
	 */
	private final void buildCudaDeviceConstructs() {
		partMass = new CUdeviceptr();   		partVolume = new CUdeviceptr();   
		gridMass = new CUdeviceptr();    
		partPos = new CUdeviceptr[3];			partVel = new CUdeviceptr[3];
		gridVel = new CUdeviceptr[3];			gridNewVel = new CUdeviceptr[3];			gridForce = new CUdeviceptr[3];
		partElasticF = new CUdeviceptr[3][];
		partPlasticF = new CUdeviceptr[3][];
		for(int i=0;i<partPos.length;++i) {			
			partPos[i] = new CUdeviceptr(); 
			partVel[i] = new CUdeviceptr();		
		
			gridVel[i] = new CUdeviceptr(); 		   
			gridNewVel[i] = new CUdeviceptr();  	   
			gridForce[i] = new CUdeviceptr();
			partElasticF[i] = new CUdeviceptr[3];
			partPlasticF[i] = new CUdeviceptr[3];
			
			for(int j=0;j<partPos.length;++j) {
				partElasticF[i][j] = new CUdeviceptr();
				partPlasticF[i][j] = new CUdeviceptr();
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
	 * @return
	 */
	//protected final int[] createSphere(TreeMap<String, ArrayList<float[]>> partVals, float ballRad, int numParts, float[] ctr) {   		 
	protected final int[] createSphere(TreeMap<String, ArrayList<float[]>> partVals, float ballRad, int numParts, myPointf ctr) {   		 
		float[] minVals = partVals.get("minMaxVals").get(0);
		float[] maxVals = partVals.get("minMaxVals").get(1); 
		ArrayList<float[]> posMap = partVals.get("pos");
		int[] returnIdxs = new int[2];
		//start at beginning of current posMap
		returnIdxs[0] = posMap.size();
		for (int i=0;i<numParts;++i) {
			//float[] posVals = getRandPosInSphereAra(ballRad, ctr); 
			float[] posVals = Base_DispWindow.AppMgr.getRandPosInSphere(ballRad, ctr).asArray(); 
			//find min/max values for all sphere particles
			for (int v = 0; v < 3; ++v) {
				if (posVals[v] < minVals[v]) {					minVals[v] = posVals[v];				} 
				else if (posVals[v] > maxVals[v]) {				maxVals[v] = posVals[v];				}	
			}
			posMap.add(posVals);
        }
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "createSphere",
				"Created a sphere of radius " + ballRad + " with "+numParts+" particles, centered at [" +ctr.toStrBrf() + "].");
		//Ending at final size of posMap
		returnIdxs[1] = posMap.size();
		return returnIdxs;
	}//createSphere
	
	protected final void setPartInitVelocities(TreeMap<String, ArrayList<float[]>> partVals, int stIdx, int endIdx, float[] initVel) {
		ArrayList<float[]> velMap = partVals.get("vel");
		for (int i=stIdx; i<endIdx;++i) {			velMap.add(initVel);		}
	}//setPartInitVelocities

	
	/**
	 * launch the kernel specified by the string key
	 * @param key string key of kernel to launch and kernel function grid, block and mem dims to access
	 */
	protected void launchKernel(String key) {
		//Set context properly before launching kernels
		JCudaDriver.cuCtxSetCurrent(context);
		int[][] kernelDims = funcGridDimAndMemSize.get(key);
        cuLaunchKernel(cuFuncs.get(key), 										// Kernel function
        		kernelDims[0][0],kernelDims[0][1],kernelDims[0][2],           	// Grid XYZ dimensions
        		kernelDims[1][0],kernelDims[1][1],kernelDims[1][2],				// Block XYZ dimensions
        		kernelDims[2][0],												// Shared memory size
        		null, 															// Stream 
                kernelParams.get(key), null);									// Kernel- and extra parameters
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
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Start CUDA Init.");
	 	//Re initialize maps of parameters and functions
        kernelParams = new TreeMap<String, Pointer>();
        funcGridDimAndMemSize = new HashMap<String, int[][]>();
        Pointer numPartsPtr = Pointer.to(new int[] {numParts});
        Pointer partMassPtr = Pointer.to(partMass);
        Pointer partVolPtr = Pointer.to(partVolume);
        Pointer gridMassPtr = Pointer.to(gridMass);
        Pointer[] partPosPtrAra = new Pointer[]{Pointer.to(partPos[0]),Pointer.to(partPos[1]),Pointer.to(partPos[2])};
        Pointer[] partVelPtrAra = new Pointer[]{Pointer.to(partVel[0]),Pointer.to(partVel[1]),Pointer.to(partVel[2])};
        
        Pointer cellSizePtr = Pointer.to(new float[] {cellSize});
        Pointer numGridSidePtr = Pointer.to(new int[]{gridSideCount});
        Pointer ttlGridCountPtr = Pointer.to(new int[] {ttlGridCount});
        Pointer[] gridVelPtrAra = new Pointer[]{Pointer.to(gridVel[0]),Pointer.to(gridVel[1]),Pointer.to(gridVel[2])};
        Pointer[] gridNewVelPtrAra = new Pointer[]{Pointer.to(gridNewVel[0]),Pointer.to(gridNewVel[1]),Pointer.to(gridNewVel[2])};
        Pointer[] gridFrcPtrAra = new Pointer[]{Pointer.to(gridForce[0]),Pointer.to(gridForce[1]),Pointer.to(gridForce[2])};
        Pointer[] gravPtrAra = new Pointer[] {Pointer.to(new float[] {gravity[0]}),Pointer.to(new float[] {gravity[1]}),Pointer.to(new float[] {gravity[2]})};
        
        Pointer minSimBndsPtr = Pointer.to(new float[] {minSimBnds});
        Pointer maxSimBndsPtr = Pointer.to(new float[] {maxSimBnds});
        Pointer deltaTPtr = Pointer.to(new float[] {deltaT});
        Pointer wallFricPtr = Pointer.to(new float[] {wallFric});
        

        kernelParams.put("projectToGridandComputeForces",Pointer.to(
        		numPartsPtr, numGridSidePtr, cellSizePtr, minSimBndsPtr,
        		Pointer.to(mat.getLambda0Ptr()), Pointer.to(mat.getMu0Ptr()), Pointer.to(mat.getHardeningCoeffPtr()),
        		partMassPtr, partVolPtr,
				partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
				partVelPtrAra[0], partVelPtrAra[1], partVelPtrAra[2],
    			//elastic matrix
				Pointer.to(partElasticF[0][0]), Pointer.to(partElasticF[0][1]), Pointer.to(partElasticF[0][2]),
				Pointer.to(partElasticF[1][0]), Pointer.to(partElasticF[1][1]), Pointer.to(partElasticF[1][2]),
				Pointer.to(partElasticF[2][0]), Pointer.to(partElasticF[2][1]), Pointer.to(partElasticF[2][2]),
				//plastic matrix
				Pointer.to(partPlasticF[0][0]), Pointer.to(partPlasticF[0][1]), Pointer.to(partPlasticF[0][2]),
				Pointer.to(partPlasticF[1][0]), Pointer.to(partPlasticF[1][1]), Pointer.to(partPlasticF[1][2]),
				Pointer.to(partPlasticF[2][0]), Pointer.to(partPlasticF[2][1]), Pointer.to(partPlasticF[2][2]),

				gridMassPtr,
				gridVelPtrAra[0], gridVelPtrAra[1], gridVelPtrAra[2],
				gridFrcPtrAra[0], gridFrcPtrAra[1], gridFrcPtrAra[2]));
        putFuncGridMemSize("projectToGridandComputeForces", part4GridDims, blkThdDims, shrdMemSize);
		
        kernelParams.put("projectToGridInit", Pointer.to(
				numPartsPtr, numGridSidePtr, cellSizePtr, minSimBndsPtr,
				partMassPtr,
				partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
				gridMassPtr));
        putFuncGridMemSize("projectToGridInit", partGridDims, blkThdDims, shrdMemSize);
        
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
		
        kernelParams.put("updDeformationGradient", Pointer.to(
        		numPartsPtr, numGridSidePtr,deltaTPtr, cellSizePtr, minSimBndsPtr,
        		Pointer.to(mat.getCriticalCompressionPtr()), Pointer.to(mat.getCriticalStretchPtr()),
				partPosPtrAra[0], partPosPtrAra[1], partPosPtrAra[2],
    			//elastic matrix
				Pointer.to(partElasticF[0][0]), Pointer.to(partElasticF[0][1]), Pointer.to(partElasticF[0][2]),
				Pointer.to(partElasticF[1][0]), Pointer.to(partElasticF[1][1]), Pointer.to(partElasticF[1][2]),
				Pointer.to(partElasticF[2][0]), Pointer.to(partElasticF[2][1]), Pointer.to(partElasticF[2][2]),
				//plastic matrix
				Pointer.to(partPlasticF[0][0]), Pointer.to(partPlasticF[0][1]), Pointer.to(partPlasticF[0][2]),
				Pointer.to(partPlasticF[1][0]), Pointer.to(partPlasticF[1][1]), Pointer.to(partPlasticF[1][2]),
				Pointer.to(partPlasticF[2][0]), Pointer.to(partPlasticF[2][1]), Pointer.to(partPlasticF[2][2]),
				//results
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
   
        win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Finished CUDA Init | Launch first MPM Pass.");
	 	//launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {launchKernel(initStepFuncKeys[j]);	}

    	simFlags.setSimIsBuilt(true);
    	win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Finished first MPM Pass.");
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
	 * @param krnFileName Cuda source code kernel file name
	 * @param ptxFileName File name for output Cuda ptx file.
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
		win.getMsgObj().dispInfoMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","Executing\n" + command);
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
			win.getMsgObj().dispErrorMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","nvcc process error : exitValue : " + exitValue);
			win.getMsgObj().dispErrorMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","errorMessage :\n" + errorMessage);
			win.getMsgObj().dispErrorMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","outputMessage :\n" + outputMessage);
			throw new IOException("Could not create .ptx file: " + errorMessage);
		}
		win.getMsgObj().dispInfoMessage("Base_MPMCudaSim("+simName+")", "compilePtxFile","Finished compiling PTX file : "+ ptxFileName);
	}//compilePtxFile
	
	/**
	 * Instance-specific per-sim cycle simulation execution code
	 * @param modAmtMillis
	 * @return
	 */
	@Override
	protected final boolean simMe_Indiv(float modAmtMillis) {
		//for every sim step, launch each kernel by key specified in simStepFuncKeys
		for (int j=0;j<simStepFuncKeys.length;++j) {	        	launchKernel(simStepFuncKeys[j]);		}
		return false;
	}
	
	/**
	 * Instance-specific post-sim cyle code
	 * @param modAmtMillis
	 * @return
	 */
	protected final boolean simMePost_Indiv(float modAmtMillis) {
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "simMePost_Indiv","Start copy relevant device buffers to host.");
		//copy from device data to host particle position or velocity arrays
 		if(simFlags.getShowParticles() || simFlags.getShowPartVels()) {
 			for(int i=0;i<hostPartPos.length;++i) {			cuMemcpyDtoH(Pointer.to(hostPartPos[i]), partPos[i], numPartsFloatSz);}
 		} 			
 		if(simFlags.getShowPartVels()) {
 			for(int i=0;i<hostPartVel.length;++i) { 		cuMemcpyDtoH(Pointer.to(hostPartVel[i]), partVel[i], numPartsFloatSz);}
 		}
 		//copy from device data to host grid velocity, accel or mass arrays
		if(simFlags.getShowGridVel()) {
 			for(int i=0;i<hostGridVel.length;++i) { 		cuMemcpyDtoH(Pointer.to(hostGridVel[i]), gridNewVel[i], numGridFloatSz);}
		}
		if(simFlags.getShowGridAccel()) {
 			for(int i=0;i<hostGridAccel.length;++i) { 		cuMemcpyDtoH(Pointer.to(hostGridAccel[i]), gridForce[i], numGridFloatSz);}		
		}		
		if(simFlags.getShowGridMass()) {					cuMemcpyDtoH(Pointer.to(hostGridMass), gridMass, numGridFloatSz);}
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "simMePost_Indiv","End copy relevant device buffers to host.");
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
	protected final void _drawParts(float animTimeMod, boolean showLocColors) {
		pa.pushMatState();
		if (showLocColors) {
			pa.drawPointCloudWithColors(hostPartPos[0].length, drawPointIncr, hostPartClrAra, hostPartPos[0], hostPartPos[1], hostPartPos[2]);
		} else {
			pa.drawPointCloudWithColors(hostPartPos[0].length, drawPointIncr, hostPartGreyAra, hostPartPos[0], hostPartPos[1], hostPartPos[2]);
		}
		pa.popMatState();
	}//_drawParts
	
	/**
	 * Draw instance class particle velocities
	 * @param animTimeMod
	 * @param pincr
	 */
	@Override
	protected final void _drawPartVel(float animTimeMod, int pincr) {
		pa.pushMatState();
		float minMag = MyMathUtils.EPS_F/vecLengthScale;
		for(int i=0;i<=hostPartVel[0].length-pincr;i+=pincr) {					
			if(		(Math.abs(hostPartVel[0][i]) > minMag) || 
					(Math.abs(hostPartVel[1][i]) > minMag) || 
					(Math.abs(hostPartVel[2][i]) > minMag)) {
				pa.pushMatState();
				pa.setStroke(hostPartClrAra[i], 255);
				pa.translate(hostPartPos[0][i], hostPartPos[1][i], hostPartPos[2][i]);
				pa.drawLine(0,0,0, vecLengthScale*hostPartVel[0][i],vecLengthScale*hostPartVel[1][i],vecLengthScale*hostPartVel[2][i]);
				pa.popMatState();
			}
		}			
		pa.popMatState();		
	}//_drawPartVel
	
	/**
	 * Draw instance class grid velocities - use _drawGridVec method
	 * @param animTimeMod
	 * @param pincr
	 */
	@Override
	protected final void _drawGridVel(float animTimeMod) {		_drawGridVec(gridVelClr, hostGridVel, hostGridPos);}

	/**
	 * Draw instance class grid accelerations - use _drawGridVec method
	 * @param animTimeMod
	 * @param pincr
	 */
	protected final void _drawGridAccel(float animTimeMod) {	_drawGridVec(gridAccelClr, hostGridAccel, hostGridPos);}
	
	/**
	 * Draw instance class grid masses - use _drawGridScalar method
	 * @param animTimeMod
	 * @param pincr
	 */
	protected final void _drawGridMass(float animTimeMod) {		_drawGridScalar(gridMassClr, hostGridMass, hostGridPos);}	
	
	/**
	 * Draw any colliders if they exist
	 * @param animTimeMod
	 */
	@Override
	protected final void _drawColliders(float animTimeMod) {
		pa.pushMatState();	
		drawColliders_Indiv(animTimeMod);
		pa.popMatState();
	}
	/**
	 * draw internal-to-sim colliders, if they exist
	 * @param animTimeMod
	 */
	protected abstract void drawColliders_Indiv(float animTimeMod);
	
}//class Base_MPMCudaSim 




