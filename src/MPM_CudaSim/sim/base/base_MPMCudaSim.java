package MPM_CudaSim.sim.base;

import static jcuda.driver.JCudaDriver.*;

import java.io.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;

import MPM_CudaSim.material.myMaterial;
import MPM_CudaSim.sim.SimResetProcess;
import MPM_CudaSim.ui.MPM_SimWindow;
import MPM_CudaSim.utils.MPM_SimUpdateFromUIData;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_UI_Objects.windowUI.base.myDispWindow;
import jcuda.*;
import jcuda.driver.*;
import jcuda.runtime.JCuda;

/**
 * abstract class describing a simulation world. Called by sim window to 
 * execute simulation and render results. Instancing classes can hold 
 * different configurations/initialization setups
 */
public abstract class base_MPMCudaSim{	

	public static IRenderInterface pa;
	//name of instancing sim
	public final String simName;
	//owning window
	protected MPM_SimWindow win;	
	//cuda kernel file name
	private String ptxFileName = "MPM_CUDA_Sim_New.ptx";
	//material quantities of particle matter
	public myMaterial mat;	
		
	/**
	 * current ui values describing variables used in the simulation
	 */
	public MPM_SimUpdateFromUIData currUIVals;	
	
	//const matrix for calculations - z is up/down
	protected final float[] gravity = new float[] {0, 0, -9.8f};
	//scale amount for visualization to fill cube frame in 3d world; particle radius, scaled by different visual scales
	protected float sclAmt;
	
	//flags relevant to simulator execution
	protected int[] simFlags;	
	public static final int
		debugSimIDX 			= 0,
		simIsBuiltIDX			= 1,			//simulation is finished being built
		showLocColors			= 2,			//display particles by color of their initial location
		showCollider			= 3,			//display visual rep of collider object, if one exists
		showParticles			= 4,			//display visual rep of material points
		showParticleVelArrows 	= 5,			//plot velocity arrows for each particle	
		showGrid				= 6,			//plot the computational grid
		showGridVelArrows 		= 7,			//plot velocity arrows for each gridNode
		showGridAccelArrows 	= 8,			//plot acceleration arrows for each gridNode
		showGridMass  			= 9,			//plot variable sized spheres proportional to gridnode mass
		showActiveNodes		  	= 10,			//show the grid nodes influenced by each particle
		CUDADevInit				= 11;			//if 1 time cuda device and kernel file init is complete
	protected static final int numSimFlags = 12;

	//time of current process start, from initial construction of mapmgr - TODO use this to monitor specific process time elapsed.  set to 0 at beginning of a particular process, then measure time elapsed in process
	protected long simStartTime;
	//time mapMgr built, in millis - used as offset for instant to provide smaller values for timestamp
	protected final long expMgrBuiltTime;	
	//constants for collider calcs
	protected static float cyl_da = (float) (Math.PI/18.0f);	
	
	////////////////////////////////////////////////////
    // Maps holding CUDA function points and parameter pointers
	// This facilitates CUDA calcs and access to appropriately configured CUDA args
	protected TreeMap<String, Pointer> kernelParams;
	protected TreeMap<String, CUfunction> cuFuncs;
	/**
	 * Lists of names of kernel functions to perform for MPM algorithm.   
	 */
	protected String[] CUFileFuncNames = new String[] {"projectToGridandComputeForces","projectToGridInit", "computeVol", "updPartVelocities",  "compGridVelocities", "partCollAndUpdPos", "gridCollisions", "clearGrid", "updDeformationGradient" };
	/**
	 * These are the individual kernel functions to execute in order for initialization of simulation
	 */
	protected String[] initStepFuncKeys = new String[] {"clearGrid","projectToGridInit", "computeVol","compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "partCollAndUpdPos"}; 
	/**
	 * These are the individual kernel functions to execute in order for each sim step
	 */
	protected String[] simStepFuncKeys = new String[] {"clearGrid", "projectToGridandComputeForces", "compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "partCollAndUpdPos"}; 
   
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
    //total # of grid cells in grid
    protected int ttlGridCount;
    //# of cuda blocks to use for particles
    protected int numBlocksParticles;
    //# of cuda blocks to use for grid
    protected int numBlocksGrid;
	//# parts * size of float & num grid cells * size float 
	protected long numPartsFloatSz, numGridFloatSz;
	//# of cuda threads
    protected final int numCUDAThreads=128;
    
    //Raw initial particle values
    TreeMap<String, ArrayList<float[]>> partVals;
    
    ////////////////////////////////////////////////////
    //representations for rendering
    //local representation of particle vector quantities for rendering
    protected float[][] h_part_pos, h_part_vel;
    //local representation of grid vector quantities for rendering
    protected float [][] h_grid_pos, h_grid_vel, h_grid_accel;
    //local rep of grid scalars for rendering
	protected float[] h_grid_mass;     
    //particle colors based on initial location
    protected int[][] h_part_clr_int;
    
    ////////////////////////////////////////////////////
    //Sim instance variables, populated from currUIVals structure on creation/ui update
    
	//# of snowballs - all particles evenly distrbuted amongst this many snowballs
	protected int numSnowballs;
	//# of particles total in the sim
	protected int numParts;
	//sim iterations per frame
    protected int simStepsPerFrame;	
	//# of cells per side - cube so same in all 3 dims; # of particles in sim
	protected int gridCount;
	//Fraction of points to draw (every x'th point will be drawn)
	protected int drawPointIncr;
	
	//timestep of simulation - 
	protected float deltaT;
	//Initial particle velocities
	protected float initVel;
	//mass of particles
	protected float particleMass;
	//size of single dimension of grid cell
	protected float cellSize;
	//simulation boundaries - symmetric cube, only need min and max, grid length per dim
	protected float minSimBnds, maxSimBnds, gridDim;
	//friction coefficients of wall and TODO:colliders
    protected float wallFric;
    //TODO support this
    protected float collFric;
	//Scaling value to scale drawn vectors
	protected float vecLengthScale;
    	
	//grid count per side - center grid always in display; grid cell dim per side
	//@SuppressWarnings("unchecked")
	public base_MPMCudaSim(IRenderInterface _pa, MPM_SimWindow _win, String _simName, MPM_SimUpdateFromUIData _currUIVals) {		
		pa=_pa;win=_win;simName = _simName;		
		currUIVals = new MPM_SimUpdateFromUIData(win);
		//initialize cuda device pointers
		buildCudaDeviceConstructs();
	
		//for display of time since experiment was built  
		Instant now = Instant.now();
		expMgrBuiltTime = now.toEpochMilli();//milliseconds since 1/1/1970 when this exec was built.
		//mat's quantities are managed by UI - only need to instance once
		mat = new myMaterial(_currUIVals);
		//initialize active nodes set - array of sets, array membership is node ID % numThreadsAvail
		//setup flag array
		initSimFlags();		
		//redundant, but placed to specify that cuda kernel needs to be loaded
		setSimFlags(CUDADevInit,false);
		//Hold sim setup particle values
		partVals = new TreeMap<String, ArrayList<float[]>>();
		
		//set up grid and initialize sim with UI values and reset sim
		updateSimVals_FromUI(_currUIVals);
	}//MPM_ABS_Sim
	
	/**
	 * Update the simulator values whenever UI values change
	 * @param upd
	 */
	public final void updateSimVals_FromUI(MPM_SimUpdateFromUIData upd) {
		
		//Determine level of reset required
		//Some UI changes should not require any changes to sim, others might require only remaking the kernels
		SimResetProcess procToDo = checkValuesForChanges(upd);
		///////////////////////
		// Require complete rebuild
		// SimResetProcess.RebuildSim
		
		//# of grid cells per side of cube
		gridCount = upd.getGridCellsPerSide();
		//# of snowballs to make
		numSnowballs = upd.getNumSnowballs();
		//cell size
		cellSize = upd.getGridCellSize();	
		
		////////////////////////////
		// Require remaking objects, but can retain # of objects, centers and base velocities
		// SimResetProcess.ResetSim
		//# of particles requested to have in simulator
		numParts = upd.getNumParticles();
		//Initial velocity of spheres
		initVel = upd.getInitVel();
				
		////////////////////////////
		// Requires kernels to be remade but no changes to sim environment
		// SimResetProcess.RemakeKernel
		//time step
		deltaT = upd.getTimeStep();		
		//particle mass to use
		particleMass = upd.getPartMass();	
		// wall friction
		wallFric = upd.getWallFricCoeff();
		// non-wall collider friction : TODO not yet supported
		collFric = upd.getCollFricCoeff();
		// update materials to match ui values
		mat.updateMatVals_FromUI(upd);				
				
		////////////////////////////
		// Requires no special handling or sim modification
		// SimResetProcess.DoNothing
		//# of simulation steps to perform between renders
		simStepsPerFrame = upd.getSimStepsPerFrame();
		//every x'th point to draw - use higher values to draw fewer points
		drawPointIncr = upd.getDrawPtIncr();
		//drawn vector/scalar quantities scaling
		vecLengthScale = upd.getDrawnVecScale();		
		
		////////////////////////////
		// calculated derived values dependent on UI values
		maxSimBnds = (gridCount*cellSize)/2.0f;
		minSimBnds = -maxSimBnds;
		gridDim = maxSimBnds - minSimBnds;	
		//scale amount to fill 1500 x 1500 x 1500 visualization cube
		sclAmt = myDispWindow.AppMgr.gridDimX/(gridCount * cellSize);
		//float size of particle arrays
		numPartsFloatSz = numParts * Sizeof.FLOAT; 
		//# cuda blocks for particle functions
		numBlocksParticles = numParts/numCUDAThreads+1;
		//total grid size       
        ttlGridCount=gridCount*gridCount*gridCount;
        //init grid ptrs        
        numGridFloatSz = ttlGridCount * Sizeof.FLOAT;
        //# cuda blocks for grid based functions			
		numBlocksGrid = ttlGridCount/numCUDAThreads+1;			
		
		//update instancing sim values
		updateSimVals_FromUI_Indiv(upd);
				
		//copy UI data to local var - copy to local last so that values that have changed can be observed
		currUIVals.setAllVals(upd);
		
		//UI changes forced reset of the simulator
		resetSim(procToDo);		
		
	}//updateSimVals_FromUI
	
	protected SimResetProcess checkValuesForChanges(MPM_SimUpdateFromUIData upd) {
		HashMap<Integer,Integer> IntIdxsToCheck = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> FloatIdxsToCheck = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> BoolIdxsToCheck = new HashMap<Integer,Integer>();

		IntIdxsToCheck.put(MPM_SimWindow.gIDX_GridCount, MPM_SimWindow.gIDX_GridCount);
		IntIdxsToCheck.put(MPM_SimWindow.gIDX_NumSnowballs, MPM_SimWindow.gIDX_NumSnowballs);
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_GridCellSize, MPM_SimWindow.gIDX_GridCellSize);
		
		if(upd.checkPassedValuesChanged(currUIVals, IntIdxsToCheck, FloatIdxsToCheck, BoolIdxsToCheck)) {
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.RebuildSim");
			return SimResetProcess.RebuildSim;
		}
		
		IntIdxsToCheck.clear();
		FloatIdxsToCheck.clear();

		IntIdxsToCheck.put(MPM_SimWindow.gIDX_NumParticles, MPM_SimWindow.gIDX_NumParticles);
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_InitVel, MPM_SimWindow.gIDX_InitVel);		
		if(upd.checkPassedValuesChanged(currUIVals, IntIdxsToCheck, FloatIdxsToCheck, BoolIdxsToCheck)) {
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.ResetSim");
			return SimResetProcess.ResetSim;
		}
		
		
		boolean matsHaveChanged = upd.haveMaterialValsChanged(currUIVals);
		if(matsHaveChanged) {
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "checkValuesForChanges","Materials have changed; Specifying SimResetProcess.RemakeKernel");
			return SimResetProcess.RemakeKernel;			
		}
		IntIdxsToCheck.clear();
		FloatIdxsToCheck.clear();
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_TimeStep, MPM_SimWindow.gIDX_TimeStep);		
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_PartMass, MPM_SimWindow.gIDX_PartMass);		
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_wallFricCoeff, MPM_SimWindow.gIDX_wallFricCoeff);		
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_CollFricCoeff, MPM_SimWindow.gIDX_CollFricCoeff);		
		
		if(upd.checkPassedValuesChanged(currUIVals, IntIdxsToCheck, FloatIdxsToCheck, BoolIdxsToCheck)) {
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.RemakeKernel");
			return SimResetProcess.RemakeKernel;
		}
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.DoNothing");
		return SimResetProcess.DoNothing;
	}
	
	/**
	 * Update instancing class variables based on potential UI changes
	 * @param upd
	 */
	protected abstract void updateSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd);
	
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
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "initOnceCUDASetup","\tRegistering Kernel Function Key : " + key);
			CUfunction c = new CUfunction();			
			cuModuleGetFunction(c, module, key);
			cuFuncs.put(key,  c);
		}
    
        setSimFlags(CUDADevInit, true);
	}//loadModuleAndSetFuncPtrs
		
	private int getClrValInt(Float val, Float min, Float max) {
		Float denom = (max-min);
		if(denom <=0) return 255;
		return (int) (55.0f + 200.0f * (val - min)/denom);	
	}
	
	/**
	 * allocate dev mem for all objects based on number of particles
	 */
	private void initCUDAMemPtrs_Parts() {
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
        for(int i = 0; i < numParts; ++i){
        	h_part_mass[i] = particleMass;
        	posAra = partVals.get("pos").get(i);
        	velAra = partVals.get("vel").get(i);        	
        	for(int j=0;j<h_part_pos.length;++j) {
        		h_part_pos[j][i] = posAra[j];
        		h_part_vel[j][i] = velAra[j];
        	}
        	for(int j=0;j<h_part_clr_int[i].length;++j) {
        		h_part_clr_int[i][j] = getClrValInt(h_part_pos[j][i],minVals[j],maxVals[j]);
        	}
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
	private void initCUDAMemPtrs_Grids() {   		
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
	}//initCUDAMemPtrs_Grids
	
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
	//return a float array of random positions within a sphere of radius rad at ctr 
	protected float[] getRandPosInSphereAra(float rad, float[] ctr){
		myVectorf pos = new myVectorf();
		ThreadLocalRandom rnd = ThreadLocalRandom.current();
		double u = rnd.nextDouble(0,1);
		Float r = (float) (rad * Math.pow(u, lcl_third));
		do{
			pos.set(rnd.nextDouble(-1,1), rnd.nextDouble(-1,1),rnd.nextDouble(-1,1));
		} while (pos.sqMagn > 1.0f);
		float[] posRes = pos.asArray();
		for(int i=0;i<3;++i) {
			posRes[i] *= r;
			posRes[i] += ctr[i];
		}
		return posRes;
	}//getRandPosInSphereAra
	
	/**
	 * Get a random location within +/- bound to place sphere. bound should 
	 * take into account desired sphere's radius, so as to not breach collider box.
	 * @param bound furthest +/- value per axis to locate center of sphere.
	 * @return a 3D center coord
	 */
	protected myVectorf getRandSphereCenter(float bound) {
		myVectorf ctr = new myVectorf(
				ThreadLocalRandom.current().nextDouble(-1,1), 
				ThreadLocalRandom.current().nextDouble(-1,1),
				ThreadLocalRandom.current().nextDouble(-1,1));
		ctr._mult(bound);
		return ctr;
	}//getRandSphereCenter
		
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
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "createSphere",
				"Created a sphere of radius " + ballRad + " with "+numParts+" particles, centered at [" +ctr[0] +"," +ctr[1] +"," +ctr[2] + "].");
		returnIdxs[1] = posMap.size();
		return returnIdxs;
	}//createSphere
	
	protected final void setPartInitVelocities(TreeMap<String, ArrayList<float[]>> partVals, int stIdx, int endIdx, float[] initVel) {
		ArrayList<float[]> velMap = partVals.get("vel");
		for (int i=stIdx; i<endIdx;++i) {
			velMap.add(initVel);
		}
	}//setPartInitVelocities

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
	 * Reinitialize partVals map
	 */
	protected final void initPartValsMap() {
		partVals.clear();
        partVals.put("pos",new ArrayList<float[]>());
        partVals.put("vel",new ArrayList<float[]>());
        //Used to map colors of particles
        partVals.put("minMaxVals",new ArrayList<float[]>());
        //initialize min and max values
        partVals.get("minMaxVals").add(new float[] {100000.0f,100000.0f,100000.0f});
        partVals.get("minMaxVals").add(new float[] {-100000.0f,-100000.0f,-100000.0f}); 
	}
	
	private void cudaSetup() {    
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
   
        win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "cudaSetup","Finished CUDA Init | Launch first MPM Pass.");
	 	//launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {
        	//win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "cudaSetup","\tinitStepFuncKeys Handle : "+initStepFuncKeys[j]);
        	launchKernel(initStepFuncKeys[j]);	}

    	setSimFlags(simIsBuiltIDX, true);
    	win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "cudaSetup","Finished first MPM Pass.");
	}//cudaSetup
	
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
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "compilePtxFile","Executing\n" + command);
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
			win.getMsgObj().dispErrorMessage("MPM_Abs_CUDASim:"+simName, "compilePtxFile","nvcc process error : exitValue : " + exitValue);
			win.getMsgObj().dispErrorMessage("MPM_Abs_CUDASim:"+simName, "compilePtxFile","errorMessage :\n" + errorMessage);
			win.getMsgObj().dispErrorMessage("MPM_Abs_CUDASim:"+simName, "compilePtxFile","outputMessage :\n" + outputMessage);
			throw new IOException("Could not create .ptx file: " + errorMessage);
		}
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "compilePtxFile","Finished compiling PTX file : "+ ptxFileName);
	}//compilePtxFile

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
    
    private final void buildSimEnvironment(SimResetProcess rebuildSim) {
        //initialize all particle values
		//Build the particle layout for this simulation
		initPartValsMap();     
		
		//determine sim-specific particle layouts
		if (rebuildSim == SimResetProcess.RebuildSim) {
			//call this if we wish to rebuild simulation layout
			buildPartLayoutMap(partVals);
		} else {
			//call this if we want to reinitialize existing simulation configuration
			reinitSimObjects(partVals);
		}      	
    }//handlePartValsMap
    
    
    /**
     * Reset simulation environment
     * @param rebuildSim Should entire sim environment be rebuilt?
     */
	public final void resetSim(SimResetProcess rebuildSim) {
		if (rebuildSim == SimResetProcess.DoNothing) {
			return;
		}
		//stop simulation and reset
		myDispWindow.AppMgr.setSimIsRunning(false);	
		//setSimFlags(CUDADevInit,false);
		setSimFlags(simIsBuiltIDX, false);	
		//if only remaking kernel, don't rebuild simulation environment
		if (rebuildSim != SimResetProcess.RemakeKernel) {
			buildSimEnvironment(rebuildSim);
		}       
		
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
		simStartTime = getCurTime();	
		
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "resetSim","Start resetting sim");
		
		if (!getSimFlags(CUDADevInit)) {
            //init cuda device and kernel file if not done already - only do 1 time
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "resetSim","CUDA Module load/init");
            initCUDAModuleSetup();
        }

		//Set context
		JCudaDriver.cuCtxSetCurrent(context);    
        //init ptrs to particle-based arrays - numparts and numPartsFloatSz need to be initialized by here
        initCUDAMemPtrs_Parts();		
        //initialize pointers to grid-based arrays
        initCUDAMemPtrs_Grids();		
		//rebuild cuda kernel configurations
		cudaSetup();

		setSimFlags(simIsBuiltIDX, true);
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "resetSim","Finished resetting sim");
	}//resetSim

	
	//boolean flag handling
	protected void initSimFlags(){simFlags = new int[1 + numSimFlags/32]; for(int i = 0; i<numSimFlags; ++i){setSimFlags(i,false);}}
	public boolean getSimFlags(int idx){int bitLoc = 1<<(idx%32);return (simFlags[idx/32] & bitLoc) == bitLoc;}	
	public void setSimFlags(int idx, boolean val) {
		boolean curVal = getSimFlags(idx);
		if(val == curVal) {return;}
		int flIDX = idx/32, mask = 1<<(idx%32);
		simFlags[flIDX] = (val ?  simFlags[flIDX] | mask : simFlags[flIDX] & ~mask);
		switch(idx){
			case debugSimIDX 			: {break;}	
			case simIsBuiltIDX 			: {break;}			
			case showLocColors 			: {break;}
			case showCollider 			: {break;}
			case showParticles			: {break;}
			case showParticleVelArrows 	: {break;}
			case showGrid				: {break;}				
			case showGridVelArrows 		: {break;}		
			case showGridAccelArrows 	: {break;}
			case showGridMass  			: {break;}		
			case showActiveNodes  		: {break;}
			default :{}			
		}			
	}//setExecFlags
	
	/**
	 * Execute simStepsPerFrame step of simulation
	 * @param modAmtMillis is in milliseconds, counting # of ms since last sim call 
	 * @return true when simulation run is complete - true turns run sim flag off
	 */	
	public final boolean simMe(float modAmtMillis) {
 		for (int i=0;i<simStepsPerFrame;++i) { 	        
 			//for every sim step, launch each kernel by key specified in simStepFuncKeys
			for (int j=0;j<simStepFuncKeys.length;++j) {	        	launchKernel(simStepFuncKeys[j]);		}
 		}
 		
		//copy from device data to host particle position or velocity arrays
 		if(getSimFlags(showParticles) || getSimFlags(showParticleVelArrows)) {
 			for(int i=0;i<h_part_pos.length;++i) {		cuMemcpyDtoH(Pointer.to(h_part_pos[i]),part_pos[i], numPartsFloatSz);}
 		} 		
		
 		if(getSimFlags(showParticleVelArrows)) {
 			for(int i=0;i<h_part_vel.length;++i) { 		cuMemcpyDtoH(Pointer.to(h_part_vel[i]),part_vel[i], numPartsFloatSz);}
 		}
 		//copy from device data to host grid velocity, accel or mass arrays
		if(getSimFlags(showGridVelArrows)) {
 			for(int i=0;i<h_grid_vel.length;++i) { 		cuMemcpyDtoH(Pointer.to(h_grid_vel[i]),grid_newvel[i], numGridFloatSz);}
		}
		if(getSimFlags(showGridAccelArrows)) {
 			for(int i=0;i<h_grid_accel.length;++i) { 	cuMemcpyDtoH(Pointer.to(h_grid_accel[i]),grid_force[i], numGridFloatSz);}		
		}		
		if(getSimFlags(showGridMass)) {
			cuMemcpyDtoH(Pointer.to(h_grid_mass),grid_mass, numGridFloatSz);	
		}

		return false;
	}//simMe
	
	/**
	 * sim method to show execution time and debug information for each sim step
	 * @param modAmtMillis
	 * @return
	 */
	public abstract boolean simMeDebug(float modAmtMillis);	//simMeDebug	

	public void showTimeMsgSimStart(String _str) {win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "showTimeMsgSimStart",_str+" Time Now : "+(getCurTime() - simStartTime));}
	//display message and time now
	public void showTimeMsgNow(String _str, long stTime) {	win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim:"+simName, "showTimeMsgNow",_str+" Time Now : "+(getCurTime() - stTime)+" ms");}
	
	/////////////////////////////
	// utility : time stamp; display messages; state flags
	
	//get time from "start time" (ctor run for map manager)
	protected long getCurTime() {			
		Instant instant = Instant.now();
		return instant.toEpochMilli() - expMgrBuiltTime;//milliseconds since 1/1/1970, subtracting when sim was built to keep millis low		
	}//getCurTime() 	
	//returns a positive int value in millis of current world time since sim start
	protected long getCurSimTime() {	return getCurTime() - simStartTime;}
	protected String getTimeStrFromProcStart() {return  getTimeStrFromPassedMillis(getCurSimTime());}
	//get a decent display of passed milliseconds elapsed
	//	long msElapsed = getCurRunTimeForProc();
	protected String getTimeStrFromPassedMillis(long msElapsed) {
		long ms = msElapsed % 1000, sec = (msElapsed / 1000) % 60, min = (msElapsed / 60000) % 60, hr = (msElapsed / 3600000) % 24;	
		String res = String.format("%02d:%02d:%02d.%03d", hr, min, sec, ms);
		return res;
	}//getTimeStrFromPassedMillis	
	

	///////////////////////////
	// end message display functionality

	
	private final int[] gridVecClr = new int[] {250,0,0}, gridAccelClr = new int[] {0,0,250}, gridMassClr = new int[] {0,150,50};
	private final int[] whitePoints = new int[] {255,255,255}; 

	//draw 1 frame of results	//animTimeMod is in seconds, counting # of seconds since last draw
	public final void drawMe(float animTimeMod) {
		if(!getSimFlags(simIsBuiltIDX)) {return;}//if not built yet, don't try to draw anything
		//render all particles - TODO determine better rendering method
		pa.pushMatState();
		//set stroke values and visual scale
			pa.setStrokeWt(2.0f/sclAmt);
			pa.scale(sclAmt);	
			
			//draw material points
			if(getSimFlags(showParticles)){	
				pa.pushMatState();
				if (getSimFlags(showLocColors)) {
					pa.drawPointCloudWithColors(h_part_pos[0].length, drawPointIncr, h_part_clr_int, h_part_pos[0], h_part_pos[1], h_part_pos[2]);
				} else {
					pa.drawPointCloudWithColor(h_part_pos[0].length, drawPointIncr, whitePoints, h_part_pos[0], h_part_pos[1], h_part_pos[2]);
				}
				pa.popMatState();
			} 
			
			if(getSimFlags(showParticleVelArrows)){
				_drawPointVel(drawPointIncr);
			}
			
			//draw colliders, if exist
			if(getSimFlags(showCollider)){
				_drawColliders(animTimeMod);
			}
			//if desired, draw grid
			if(getSimFlags(showGrid)) {	_drawGrid();}
			if (getSimFlags(showGridVelArrows)) {	_drawGridVec(gridVecClr, h_grid_vel[0], h_grid_vel[1], h_grid_vel[2]);}
			if (getSimFlags(showGridAccelArrows)){	_drawGridVec(gridAccelClr, h_grid_accel[0], h_grid_accel[1], h_grid_accel[2]);}
			if(getSimFlags(showGridMass)) {			_drawGridScalar(gridMassClr, h_grid_mass);}
		pa.popMatState();
	}//drawMe
	
	private void _drawPointVel(int pincr) {
		pa.pushMatState();
		float minMag = MyMathUtils.EPS_F/vecLengthScale;
		for(int i=0;i<=h_part_vel[0].length-pincr;i+=pincr) {					
			if(		(Math.abs(h_part_vel[0][i]) > minMag) || 
					(Math.abs(h_part_vel[1][i]) > minMag) || 
					(Math.abs(h_part_vel[2][i]) > minMag)) {
				pa.pushMatState();
				pa.setStroke(h_part_clr_int[i], 100);
				pa.translate(h_part_pos[0][i], h_part_pos[1][i], h_part_pos[2][i]);
				pa.drawLine(0,0,0, vecLengthScale*h_part_vel[0][i],vecLengthScale*h_part_vel[1][i],vecLengthScale*h_part_vel[2][i]);
				pa.popMatState();
			}
		}			
		pa.popMatState();
	}//drawPointVel
	
	/**
	 * Every 10th grid line should be drawn
	 */
	private final int gridIncr = 10;
	private void _drawGrid() {
		pa.pushMatState();		
			pa.setStroke(255,255,255,20);
			pa.translate(minSimBnds,minSimBnds,minSimBnds);
			//shows every "incr" gridcells
			for (int i=0; i<=gridCount;i+=gridIncr) {
				float iLoc = i*cellSize;
				for(int j=0;j<=gridCount;j+=gridIncr) {
					myVectorf startPos=new myVectorf(iLoc,j*cellSize,0.0f);
					myVectorf endPos=new myVectorf(iLoc, startPos.y,gridDim);
					pa.drawLine(startPos,endPos);
				}
				for(int k=0;k<=gridCount;k+=gridIncr) {
					myVectorf startPos=new myVectorf(iLoc,0.0f, k*cellSize);
					myVectorf endPos=new myVectorf(iLoc,gridDim,startPos.z);
					pa.drawLine(startPos,endPos);
				}
			}
			for(int j=0;j<=gridCount;j+=gridIncr) {
				float jLoc = j*cellSize;
				for(int k=0;k<=gridCount;k+=gridIncr) {
					myVectorf startPos=new myVectorf(0.0f,jLoc,k*cellSize);
					myVectorf endPos=new myVectorf(gridDim,jLoc,startPos.z);
					pa.drawLine(startPos,endPos);
				}
			}
		pa.popMatState();		
	}//_drawGrid()

	private void _drawGridVec(int[] clr, float[] xVal, float[] yVal, float[] zVal) {
		float minMag = MyMathUtils.EPS_F/vecLengthScale;
		pa.pushMatState();	
		pa.setStroke(clr,20);
		pa.translate(minSimBnds,minSimBnds,minSimBnds);
		for (int i=0; i<ttlGridCount;++i) {			
			if(		(Math.abs(xVal[i]) > minMag) || 
					(Math.abs(yVal[i]) > minMag) || 
					(Math.abs(zVal[i]) > minMag)) {
				pa.pushMatState();	
				pa.translate(h_grid_pos[0][i], h_grid_pos[1][i], h_grid_pos[2][i]);
				pa.drawLine(0,0,0, vecLengthScale*xVal[i],vecLengthScale*yVal[i],vecLengthScale*zVal[i]);
				pa.popMatState();
			}
		}
		pa.popMatState();
	}
	
	private void _drawGridScalar(int[] clr, float[] xVal) {
		float minMag = MyMathUtils.EPS_F/vecLengthScale;
		pa.pushMatState();	
		pa.setSphereDetail(4);
		pa.setStroke(clr,20);
		pa.translate(minSimBnds,minSimBnds,minSimBnds);
		for (int i=0; i<ttlGridCount;++i) {			
			if(		(Math.abs(xVal[i]) > minMag)) {
				pa.pushMatState();	
				pa.translate(h_grid_pos[0][i], h_grid_pos[1][i], h_grid_pos[2][i]);
				pa.drawSphere(xVal[i]*vecLengthScale);
				pa.popMatState();
			}
		}
		pa.popMatState();		
	}
	
	/**
	 * Draw any colliders if they exist
	 * @param animTimeMod
	 */
	private void _drawColliders(float animTimeMod) {
		pa.pushMatState();	
		drawColliders_Indiv(animTimeMod);
		pa.popMatState();
	}
	//draw internal-to-sim colliders, if they exist
	protected abstract void drawColliders_Indiv(float animTimeMod);
	
	//some utility functions and constants
	protected myPointf Pf(myPointf O, float x, myVectorf I, double y, myVectorf J, double z, myVectorf K) {return new myPointf(O.x+x*I.x+y*J.x+z*K.x,O.y+x*I.y+y*J.y+z*K.y,O.z+x*I.z+y*J.z+z*K.z);}  // O+xI+yJ+zK
	//returns arrays of start and end points for cylinder of passed radius r going from point a to point b
	@SuppressWarnings("unchecked")
	protected ArrayList<myPointf>[] buildCylinder(myPointf A, myPointf B, float r){
		myPointf P = A;
		myVectorf V = new myVectorf(A, B);
		myVectorf I = myVectorf.UP;
		myVectorf Nvec = new myVectorf(I.y*V.z-I.z*V.y, I.z*V.x-I.x*V.z, I.x*V.y-I.y*V.x);
		if(Math.abs(Nvec.magn) < 0.000001) {//singular - cylinder wanting to go up
			System.out.println("vec singlr : "+ Nvec.magn + " V:"+V.toStrBrf()+" | I : " + I.toStrBrf());
			I = myVectorf.RIGHT;
			Nvec = new myVectorf(I.y*V.z-I.z*V.y, I.z*V.x-I.x*V.z, I.x*V.y-I.y*V.x);// Nf(I,V);
		}
		myVectorf J = Nvec._normalize();
		ArrayList<myPointf>[] res = new ArrayList[2];
		for(int idx =0;idx<res.length;++idx) {res[idx]=new ArrayList<myPointf>();}
		float rcA, rsA;
		for(float a=0; a<=MyMathUtils.TWO_PI_F+cyl_da; a+=cyl_da) {
			rcA = (float) (r*Math.cos(a)); 
			rsA = (float) (r*Math.sin(a));
			res[0].add(Pf(P,rcA,I,rsA,J,0.0,V)); 
			res[1].add(Pf(P,rcA,I,rsA,J,1.0,V));
		}
		return res;
	}
	
}//class MPM_ABS_Sim 




