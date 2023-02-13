package MPM_CudaSim.sim.base;

import static jcuda.driver.JCudaDriver.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

import MPM_SimMain.sim.SimResetProcess;
import MPM_SimMain.sim.Base_MPMSim;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Render_Interface.IRenderInterface;
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
   
	protected HashMap<String, int[]> funcGridDimAndMemSize;
	
	////////////////////////////////////////////////////
    // CUDA references
    protected CUdevice dev;
    protected CUcontext context;
    protected CUmodule module;
    protected CUgraphicsResource pCudaResource;
	// CUDA Device ptr constructions
	protected CUdeviceptr part_mass, part_vol;
	protected CUdeviceptr[] part_pos, part_vel; 
    protected CUdeviceptr[] grid_vel, grid_newvel, grid_force;
	protected CUdeviceptr[][] part_fe, part_fp;    
    protected CUdeviceptr grid_mass;
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
    protected final int numCUDAThreads=512;
    
    /**
     * Raw initial particle values
     */
    TreeMap<String, ArrayList<float[]>> partVals;
    
    ////////////////////////////////////////////////////
    //representations for rendering
    /**
     * local representation of particle vector quantities for rendering
     */
    protected float[][] h_part_pos, h_part_vel;
    /**
     * local representation of grid vector quantities for rendering
     */
    protected float [][] h_grid_pos, h_grid_vel, h_grid_accel;
    /**
     * local rep of grid scalars for rendering
     */
	protected float[] h_grid_mass;     
	/**
	 * particle colors based on initial location
	 */
    protected int[][] h_part_clr_int, h_part_grey_int;
        
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
		JCuda.cudaFree(part_mass);		  
        JCuda.cudaFree(part_vol); 
        JCuda.cudaFree(grid_mass); 
        for(int i=0;i<part_pos.length;++i) {
        	JCuda.cudaFree(part_pos[i]);
        	JCuda.cudaFree(part_vel[i]);
	        JCuda.cudaFree(grid_vel[i]); 
	        JCuda.cudaFree(grid_newvel[i]);
	        JCuda.cudaFree(grid_force[i]);
	        
	        for(int j=0;j<part_fe[0].length;++j) {
	        	JCuda.cudaFree(part_fe[i][j]);                    
	        	JCuda.cudaFree(part_fp[i][j]);		        	
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
		//float size of particle arrays
		numPartsFloatSz = numParts * Sizeof.FLOAT; 
		//# cuda blocks for particle functions
		numBlocksParticles = numParts/numCUDAThreads+1;
        //init grid ptrs        
        numGridFloatSz = ttlGridCount * Sizeof.FLOAT;
        //# cuda blocks for grid based functions			
		numBlocksGrid = ttlGridCount/numCUDAThreads+1;			
		
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
        h_part_pos = new float[3][];
        h_part_vel = new float[3][];
        for(int i=0;i<h_part_pos.length;++i) {
        	h_part_pos[i] = new float[numParts];
        	h_part_vel[i] = new float[numParts];
        }
       
        
        float[] minVals = partVals.get("minMaxVals").get(0);
        float[] maxVals = partVals.get("minMaxVals").get(1);       
        float[] posAra, velAra;
        
        h_part_clr_int = new int[numParts][3]; 
        h_part_grey_int = new int[numParts][3]; 
        ArrayList<float[]> allPosAra = partVals.get("pos");
        ArrayList<float[]> allVelAra = partVals.get("vel");
        
        for(int i = 0; i < numParts; ++i){
        	h_part_mass[i] = particleMass;
        	posAra = allPosAra.get(i);
        	velAra = allVelAra.get(i);        	
        	for(int j=0;j<h_part_pos.length;++j) {
        		h_part_pos[j][i] = posAra[j];
        		h_part_vel[j][i] = velAra[j];
        	}
        	
        	h_part_clr_int[i] = getClrValInt(posAra,minVals,maxVals);
        	h_part_grey_int[i] = getGreyValInt(posAra,minVals,maxVals);
        	h_part_eye[i]=1.0f;
        }
        
        //Allocate memory for particle-related cuda pointers
        cuMemAlloc(part_vol, numPartsFloatSz); 
        cuMemsetD32(part_vol, 0, numParts);	//part_vol is a calculated quantity
        cuMemAlloc(part_mass, numPartsFloatSz);  
        cuMemcpyHtoD(part_mass, Pointer.to(h_part_mass), numPartsFloatSz);
        
        for(int i=0;i<part_pos.length;++i) {
           	cuMemAlloc(part_pos[i], numPartsFloatSz); 
           	cuMemAlloc(part_vel[i], numPartsFloatSz);
            cuMemcpyHtoD(part_pos[i], Pointer.to(h_part_pos[i]), numPartsFloatSz);
            cuMemcpyHtoD(part_vel[i], Pointer.to(h_part_vel[i]), numPartsFloatSz);
            for(int j=0;j<part_fe[0].length;++j) {		//build identity mats for this
                cuMemAlloc(part_fe[i][j], numPartsFloatSz);
                cuMemAlloc(part_fp[i][j], numPartsFloatSz); 
                if(i==j) {
                    cuMemcpyHtoD(part_fe[i][j],	Pointer.to(h_part_eye), numPartsFloatSz);                    
                    cuMemcpyHtoD(part_fp[i][j],	Pointer.to(h_part_eye), numPartsFloatSz);             	       	
                } else {
                    cuMemsetD32(part_fe[i][j], 0, numParts);                    
                    cuMemsetD32(part_fp[i][j], 0, numParts);             	
                }        	
            }    	        	
        }
        
    }//initCUDAMemPtrs_Parts
	
	/**
	 * allocate dev mem for all objects based on number of grid cells
	 */
	@Override
	protected final void initValues_Grids() {
		h_grid_pos = new float[3][];
		h_grid_vel = new float[3][];
		h_grid_accel = new float[3][];
		for(int i=0;i<h_part_pos.length;++i) {
			h_grid_pos[i] = new float[ttlGridCount];
			h_grid_vel[i] = new float[ttlGridCount];
	       	h_grid_accel[i] = new float[ttlGridCount];
        }
		
		h_grid_mass = new float[ttlGridCount];		
		//build grid locations
		int gridDim=0;
		for(int i=0;i<gridCount;++i) {
			float xPos = (i+.5f)*cellSize;
			for(int j=0;j<gridCount;++j) {		
				float yPos = (j+.5f)*cellSize;
				for(int k=0;k<gridCount;++k) {
					float zPos = (k+.5f)*cellSize;
					//weird orientation stuffs
					//TODO fix grid orientation
					h_grid_pos[0][gridDim] = zPos;
					h_grid_pos[1][gridDim] = yPos;
					h_grid_pos[2][gridDim] = xPos;					
					++gridDim;
				}
			}
		}			
		
		//Allocate memory and initialize values for cuda grid pointers
        cuMemAlloc(grid_mass, numGridFloatSz); 
		cuMemsetD32(grid_mass, 0, ttlGridCount);
		for(int i=0;i<grid_vel.length;++i) {
            cuMemAlloc(grid_vel[i], numGridFloatSz);  
            cuMemAlloc(grid_newvel[i], numGridFloatSz);
            cuMemAlloc(grid_force[i], numGridFloatSz);
            
            cuMemsetD32(grid_vel[i], 0, ttlGridCount);  
            cuMemsetD32(grid_newvel[i], 0, ttlGridCount);
            cuMemsetD32(grid_force[i], 0, ttlGridCount);           	
        }
	}//initValues_Grids
	
	/**
	 * Only performed from ctor
	 */
	private final void buildCudaDeviceConstructs() {
		part_mass = new CUdeviceptr();   		part_vol = new CUdeviceptr();   
		grid_mass = new CUdeviceptr();    
		part_pos = new CUdeviceptr[3];			part_vel = new CUdeviceptr[3];
		grid_vel = new CUdeviceptr[3];			grid_newvel = new CUdeviceptr[3];			grid_force = new CUdeviceptr[3];
		part_fe = new CUdeviceptr[3][];
		part_fp = new CUdeviceptr[3][];
		for(int i=0;i<part_pos.length;++i) {			
			part_pos[i] = new CUdeviceptr(); 
			part_vel[i] = new CUdeviceptr();		
		
			grid_vel[i] = new CUdeviceptr(); 		   
			grid_newvel[i] = new CUdeviceptr();  	   
			grid_force[i] = new CUdeviceptr();
			part_fe[i] = new CUdeviceptr[3];
			part_fp[i] = new CUdeviceptr[3];
			
			for(int j=0;j<part_pos.length;++j) {
				part_fe[i][j] = new CUdeviceptr();
				part_fp[i][j] = new CUdeviceptr();
			}
		}		
	}//buildCudaDeviceConstructs	
	
	private static final double lcl_third = 1.0/3.0;
	/**
	 * Find a random position in a sphere centered at ctr of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @param ctr
	 * @return
	 */
	public final myPointf getRandPosInSphere(double rad){ return getRandPosInSphere(rad, new float[3]);}
	public final myPointf getRandPosInSphere(double rad, float[] ctr){
		myPointf pos = new myPointf();
		double u = ThreadLocalRandom.current().nextDouble(0,1),	
			cosTheta = ThreadLocalRandom.current().nextDouble(-1,1),
			phi = ThreadLocalRandom.current().nextDouble(0,MyMathUtils.TWO_PI_F),
			r = rad * Math.pow(u, lcl_third),
			rSinTheta = r * Math.sin(Math.acos(cosTheta));			
		pos.set(rSinTheta * Math.cos(phi) + ctr[0], rSinTheta * Math.sin(phi) + ctr[1],cosTheta*r + ctr[2]);
		//pos._add(ctr[0],ctr[1],ctr[2]);
		return pos;
	}
	
	//return a float array of random positions within a sphere of radius rad at ctr 
	protected float[] getRandPosInSphereAra(float rad, float[] ctr){
//		myVectorf posV = new myVectorf();
//		ThreadLocalRandom rnd = ThreadLocalRandom.current();
//		double u = rnd.nextDouble(0,1);
//		Float r = (float) (rad * Math.pow(u, lcl_third));
//		do{
//			posV.set(rnd.nextDouble(-1,1), rnd.nextDouble(-1,1),rnd.nextDouble(-1,1));
//		} while (posV.sqMagn > 1.0f);
//		float[] posRes = posV.asArray();
//		for(int i=0;i<3;++i) {
//			posRes[i] *= r;
//			posRes[i] += ctr[i];
//		}
//		return posRes;
		myPointf pos = getRandPosInSphere(rad,ctr);
		return pos.asArray();
	}//getRandPosInSphereAra
		
	/**
	 * create a sphere with given center, with passed # of particles -0 returns array of [start IDX,end IDX] within posVals array for this sphere
	 * @param partVals map holding all appropriate particle values
	 * @param ballRad desired radius of sphere
	 * @param numParts
	 * @param initVel
	 * @param ctr
	 * @return
	 */
	protected final int[] createSphere(TreeMap<String, ArrayList<float[]>> partVals, float ballRad, int numParts, float[] ctr) {   		 
		float[] minVals = partVals.get("minMaxVals").get(0);
		float[] maxVals = partVals.get("minMaxVals").get(1); 
		ArrayList<float[]> posMap = partVals.get("pos");
		int[] returnIdxs = new int[2];
		//start at beginning of current posMap
		returnIdxs[0] = posMap.size();
		for (int i=0;i<numParts;++i) {
			float[] posVals = getRandPosInSphereAra(ballRad, ctr); 
			//find min/max values for all sphere particles
			for (int v = 0; v < 3; ++v) {
				if (posVals[v] < minVals[v]) {					minVals[v] = posVals[v];				} 
				else if (posVals[v] > maxVals[v]) {				maxVals[v] = posVals[v];				}	
			}
			posMap.add(posVals);
        }
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "createSphere",
				"Created a sphere of radius " + ballRad + " with "+numParts+" particles, centered at [" +ctr[0] +"," +ctr[1] +"," +ctr[2] + "].");
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
	 * @param kernelVal values for specific kernal 
	 * @param key string key of kernel to launch
	 */
	private void launchKernel(String key) {
		//Set context properly before launching kernels
		JCudaDriver.cuCtxSetCurrent(context);
		int[] kernalDims = funcGridDimAndMemSize.get(key);
        cuLaunchKernel(cuFuncs.get(key), 
        		kernalDims[0], 1, 1,           // Grid dimension 
                numCUDAThreads, 1, 1,  // Block dimension
                kernalDims[1], null,           // Shared memory size and stream 
                kernelParams.get(key), null);// Kernel- and extra parameters
        //Allow kernel to complete
        cuCtxSynchronize();
	}
	/**
	 * Set up all essential cuda kernels and launch them for initial pass
	 */
	private void cudaSetup() {    
		win.getMsgObj().dispDebugMessage("Base_MPMCudaSim("+simName+")", "cudaSetup","Start CUDA Init.");
	 	//Re initialize maps of parameters and functions
        kernelParams = new TreeMap<String, Pointer>();
        funcGridDimAndMemSize = new HashMap<String, int[]>();
        
        kernelParams.put("projectToGridandComputeForces",Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
        		Pointer.to(mat.getLambda0Ptr()), Pointer.to(mat.getMu0Ptr()), Pointer.to(mat.getHardeningCoeffPtr()),
        		Pointer.to(part_mass), Pointer.to(part_vol),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(part_vel[0]),Pointer.to(part_vel[1]),Pointer.to(part_vel[2]),
    			//elastic matrix
				Pointer.to(part_fe[0][0]), Pointer.to(part_fe[0][1]), Pointer.to(part_fe[0][2]),
				Pointer.to(part_fe[1][0]), Pointer.to(part_fe[1][1]), Pointer.to(part_fe[1][2]),
				Pointer.to(part_fe[2][0]), Pointer.to(part_fe[2][1]), Pointer.to(part_fe[2][2]),
				//plastic matrix
				Pointer.to(part_fp[0][0]), Pointer.to(part_fp[0][1]), Pointer.to(part_fp[0][2]),
				Pointer.to(part_fp[1][0]), Pointer.to(part_fp[1][1]), Pointer.to(part_fp[1][2]),
				Pointer.to(part_fp[2][0]), Pointer.to(part_fp[2][1]), Pointer.to(part_fp[2][2]),

				Pointer.to(grid_mass),
				Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),
				Pointer.to(grid_force[0]),Pointer.to(grid_force[1]),Pointer.to(grid_force[2])));
        funcGridDimAndMemSize.put("projectToGridandComputeForces", new int[] {numBlocksParticles*4, 0});
		
        kernelParams.put("projectToGridInit", Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
				Pointer.to(part_mass),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(grid_mass)));
        funcGridDimAndMemSize.put("projectToGridInit", new int[] {numBlocksParticles, 0});
		
		kernelParams.put("computeVol", Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
				Pointer.to(part_mass),Pointer.to(part_vol),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(grid_mass)));
		funcGridDimAndMemSize.put("computeVol", new int[] {numBlocksParticles, 0});
		
		kernelParams.put("updPartVelocities",Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
				Pointer.to(mat.getAlphaPicFlipPtr()),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(part_vel[0]),Pointer.to(part_vel[1]),Pointer.to(part_vel[2]),
				Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),		
				Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2])));  
		funcGridDimAndMemSize.put("updPartVelocities", new int[] {numBlocksParticles*4, numCUDAThreads*6*Sizeof.FLOAT});
		
	    kernelParams.put("compGridVelocities", Pointer.to(
				Pointer.to(new int[] {ttlGridCount}), Pointer.to(new float[] {gravity[0]}),Pointer.to(new float[] {gravity[1]}),Pointer.to(new float[] {gravity[2]}),Pointer.to(new float[] {deltaT}),
				Pointer.to(grid_mass),		
				Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),			
				Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2]),
				Pointer.to(grid_force[0]),Pointer.to(grid_force[1]),Pointer.to(grid_force[2])));
		funcGridDimAndMemSize.put("compGridVelocities", new int[] {numBlocksGrid, 0});
	    
	    kernelParams.put("clearGrid", Pointer.to(
	    		Pointer.to(new int[] {ttlGridCount}), Pointer.to(grid_mass),		
	    		Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),			
	    		Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2]),
				Pointer.to(grid_force[0]),Pointer.to(grid_force[1]),Pointer.to(grid_force[2])));
		funcGridDimAndMemSize.put("clearGrid", new int[] {numBlocksGrid, 0});
		//Only supports wall collisions
	    kernelParams.put("gridCollisions", Pointer.to(
        		Pointer.to(new int[] {ttlGridCount}),Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {cellSize}), 
        		Pointer.to(new float[] {minSimBnds}),Pointer.to(new float[] {maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {deltaT}),Pointer.to(grid_mass),
        		Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2])));
		funcGridDimAndMemSize.put("gridCollisions", new int[] {numBlocksGrid, 0});

        kernelParams.put("updDeformationGradient", Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {deltaT}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
        		Pointer.to(mat.getCriticalCompressionPtr()), Pointer.to(mat.getCriticalStretchPtr()),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
    			//elastic matrix
				Pointer.to(part_fe[0][0]), Pointer.to(part_fe[0][1]), Pointer.to(part_fe[0][2]),
				Pointer.to(part_fe[1][0]), Pointer.to(part_fe[1][1]), Pointer.to(part_fe[1][2]),
				Pointer.to(part_fe[2][0]), Pointer.to(part_fe[2][1]), Pointer.to(part_fe[2][2]),
				//plastic matrix
				Pointer.to(part_fp[0][0]), Pointer.to(part_fp[0][1]), Pointer.to(part_fp[0][2]),
				Pointer.to(part_fp[1][0]), Pointer.to(part_fp[1][1]), Pointer.to(part_fp[1][2]),
				Pointer.to(part_fp[2][0]), Pointer.to(part_fp[2][1]), Pointer.to(part_fp[2][2]),
				//results
				Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2])));     
        funcGridDimAndMemSize.put("updDeformationGradient", new int[] {numBlocksParticles, 0});  
        
        kernelParams.put("partCollAndUpdPos", Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new float[] {minSimBnds}),Pointer.to(new float[] {maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {deltaT}),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(part_vel[0]),Pointer.to(part_vel[1]),Pointer.to(part_vel[2])				
        		));   
        funcGridDimAndMemSize.put("partCollAndUpdPos", new int[] {numBlocksParticles, 0});        
   
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
 			for(int i=0;i<h_part_pos.length;++i) {		cuMemcpyDtoH(Pointer.to(h_part_pos[i]),part_pos[i], numPartsFloatSz);}
 		} 			
 		if(simFlags.getShowPartVels()) {
 			for(int i=0;i<h_part_vel.length;++i) { 		cuMemcpyDtoH(Pointer.to(h_part_vel[i]),part_vel[i], numPartsFloatSz);}
 		}
 		//copy from device data to host grid velocity, accel or mass arrays
		if(simFlags.getShowGridVel()) {
 			for(int i=0;i<h_grid_vel.length;++i) { 		cuMemcpyDtoH(Pointer.to(h_grid_vel[i]),grid_newvel[i], numGridFloatSz);}
		}
		if(simFlags.getShowGridAccel()) {
 			for(int i=0;i<h_grid_accel.length;++i) { 	cuMemcpyDtoH(Pointer.to(h_grid_accel[i]),grid_force[i], numGridFloatSz);}		
		}		
		if(simFlags.getShowGridMass()) {
			cuMemcpyDtoH(Pointer.to(h_grid_mass),grid_mass, numGridFloatSz);	
		}
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
			pa.drawPointCloudWithColors(h_part_pos[0].length, drawPointIncr, h_part_clr_int, h_part_pos[0], h_part_pos[1], h_part_pos[2]);
		} else {
			pa.drawPointCloudWithColors(h_part_pos[0].length, drawPointIncr, h_part_grey_int, h_part_pos[0], h_part_pos[1], h_part_pos[2]);
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
		for(int i=0;i<=h_part_vel[0].length-pincr;i+=pincr) {					
			if(		(Math.abs(h_part_vel[0][i]) > minMag) || 
					(Math.abs(h_part_vel[1][i]) > minMag) || 
					(Math.abs(h_part_vel[2][i]) > minMag)) {
				pa.pushMatState();
				pa.setStroke(h_part_clr_int[i], 255);
				pa.translate(h_part_pos[0][i], h_part_pos[1][i], h_part_pos[2][i]);
				pa.drawLine(0,0,0, vecLengthScale*h_part_vel[0][i],vecLengthScale*h_part_vel[1][i],vecLengthScale*h_part_vel[2][i]);
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
	protected final void _drawGridVel(float animTimeMod) {		_drawGridVec(gridVelClr, h_grid_vel, h_grid_pos);}

	/**
	 * Draw instance class grid accelerations - use _drawGridVec method
	 * @param animTimeMod
	 * @param pincr
	 */
	protected final void _drawGridAccel(float animTimeMod) {	_drawGridVec(gridAccelClr, h_grid_accel, h_grid_pos);}
	
	/**
	 * Draw instance class grid masses - use _drawGridScalar method
	 * @param animTimeMod
	 * @param pincr
	 */
	protected final void _drawGridMass(float animTimeMod) {		_drawGridScalar(gridMassClr, h_grid_mass, h_grid_pos);}	
	
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




