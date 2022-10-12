package MPM_CudaSim.sim.base;

import static jcuda.driver.JCudaDriver.*;

import java.io.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;

import MPM_CudaSim.material.myMaterial;
import MPM_CudaSim.ui.MPM_SimWindow;
import MPM_CudaSim.utils.MPM_SimUpdateFromUIData;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
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
	
	/**
	 * current ui values describing variables used in the simulation
	 */
	public MPM_SimUpdateFromUIData currUIVals;	

	//these values are just to set initial values for simulation
	public static int numGridCellsDefault = 200;
	public static float cellSizeDefault = .10f;
	public static int numPartsUI_Init = 100000;
	//snow density varies from 50 to ~800
	//powdery snow is about 100
	//wet firm compacted snow is about 600
	//public static float materialDensityDefault = 100.0f;
	
	//# of particles to have in sim
	protected int numParts = numPartsUI_Init;
	//grid dim - cube so same in all 3 dims; # of particles in sim
	protected int gridCount = numGridCellsDefault;
	protected float cellSize = cellSizeDefault;
	//# parts * size of float & num grid cells * size float 
	protected long numPartsFloatSz, numGridFloatSz;
	
	//simulation boundaries - symmetric cube, only need min and max, grid length per dim
	protected float minSimBnds, maxSimBnds, gridDim;
	//dimension of grid cells
	protected float particleMass;
	//whether current OS supports ansi terminal color settings
	public static boolean supportsANSITerm = false;
	//parameters
	//timestep of simulation - 
	protected static float deltaT = 4e-4f;

	//const matrix for calculations - z is up/down
	protected final float[] gravity = new float[] {0, 0, -9.8f};
	//scale amount for visualization to fill cube frame in 3d world; particle radius, scaled by different visual scales
	protected float sclAmt;
	
	//material quantities of particle matter
	public myMaterial mat;	
	
	//friction coefficients of colliders
	public static float wallFric = 1.0f;
	public static float collFric = 1.0f;
	//iterations per frame
	public static int simStepsPerFrame = 2;

	//flags relevant to simulator execution
	protected int[] simFlags;	
	public static final int
		debugSimIDX 			= 0,
		simIsBuiltIDX			= 1,			//simulation is finished being built - set in thread manager
		gridIsBuiltIDX			= 2,			//sim grid is built
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
	
	protected TreeMap<String, Pointer> kernelParams;
	protected TreeMap<String, CUfunction> cuFuncs;
	protected String[] CUFileFuncNames = new String[] {"projectToGridandComputeForces","projectToGridInit", "computeVol", "updPartVelocities",  "compGridVelocities", "partCollAndUpdPos", "gridCollisions", "clearGrid", "updDeformationGradient" };
	protected String[] initStepFuncKeys = new String[] {"clearGrid","projectToGridInit", "computeVol","compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "partCollAndUpdPos"}; 
	protected String[] simStepFuncKeys = new String[] {"clearGrid", "projectToGridandComputeForces", "compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "partCollAndUpdPos"}; 
   
	protected HashMap<String, int[]> funcGridDimAndMemSize;
    
	protected CUdeviceptr part_mass, part_vol;
	protected CUdeviceptr[] part_pos, part_vel; 
    protected CUdeviceptr[] grid_vel, grid_newvel, grid_force;

	protected CUdeviceptr[][] part_fe, part_fp;
    
    protected CUdeviceptr grid_mass;
    
    //local representation of particle quantities for rendering
    protected float[][] h_part_pos, h_part_vel;
    //local representation of grid quantities for rendering
    protected float [][] h_grid_pos, h_grid_vel, h_grid_accel;
    
	protected float[] h_grid_mass; 
    
    //colors based on initial location
    protected int[][] h_part_clr_int;
    
    protected CUcontext pctx;
    protected CUdevice dev;
    protected CUmodule module;
    protected CUgraphicsResource pCudaResource;
    
    protected int gridSize;
    protected int numBlocksParticles;
    protected int numBlocksGrid;
	//# of cuda threads
    protected  int numCUDAThreads=128;
	
	//grid count per side - center grid always in display; grid cell dim per side
	//@SuppressWarnings("unchecked")
	public base_MPMCudaSim(IRenderInterface _pa,MPM_SimWindow _win, String _simName, int _gridCount, float _h, int _numParts, float _density) {		
		pa=_pa;win=_win;simName = _simName;		
		
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

		//whether this system supports an ansi terminal or not
		supportsANSITerm = (System.console() != null && System.getenv().get("TERM") != null);	
		
		//for display of time since experiment was built  
		Instant now = Instant.now();
		expMgrBuiltTime = now.toEpochMilli();//milliseconds since 1/1/1970 when this exec was built.
		//mat's quantities are managed by UI - only need to instance once
		mat = new myMaterial();
		//initialize active nodes set - array of sets, array membership is node ID % numThreadsAvail
		//setup flag array
		
		initSimFlags();
		
		//set up grid and initialize sim
		setGridValsAndInit(_gridCount, _h, _numParts, _density);
		//resetSim(false);
	}//MPM_ABS_Sim
	
	public final void updateMapMorphVals_FromUI(MPM_SimUpdateFromUIData upd) {
		//TODO Check if values need to be reset
		currUIVals.setAllVals(upd);//can't use the same mapUpdFromUIData everywhere because we compare differences
	}
	
	//run 1 time to load kernel and assign function pointers to functions
	private void initCUDAModuleSetup() {
		// Enable exceptions and omit all subsequent error checks
        JCudaDriver.setExceptionsEnabled(true); 
        //Initialize the driver and create a context for the first device. (device 0
        //build maps that hold values for cuda grid dim and shared mem size
        
        cuInit(0);
        pctx = new CUcontext();
        dev = new CUdevice();
        cuDeviceGet(dev, 0);
        cuCtxCreate(pctx, 0, dev);
		// Load the ptx file.
		module = new CUmodule();
//		System.out.println("Working Directory = " +  System.getProperty("user.dir"));
//		try {
//			compilePtxFile("src\\Cuda\\MPM_ABS_Sim.cu","MPM_ABS_Sim.ptx");
//			ptxFileName = "MPM_ABS_Sim.ptx";
//			//ptxFileName = MPM_ABS_Sim.preparePtxFile("src\\Cuda\\MPM_ABS_Sim.cu");
//		} catch (Exception e) {
//			System.out.println(e.getMessage());
//		}
		
		cuModuleLoad(module, ptxFileName);    	
		
        // Obtain a function pointer to each function in cuda kernel file
		cuFuncs = new TreeMap<String, CUfunction>();
		for (int i =0;i<CUFileFuncNames.length; ++i) {
			String key = CUFileFuncNames[i];			
			win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "initOnceCUDASetup","\tRegistering Kernel Function Key : " + key);
			CUfunction c = new CUfunction();			
			cuModuleGetFunction(c, module, key);
			cuFuncs.put(key,  c);
		}
    
        setSimFlags(CUDADevInit, true);
	}//loadModuleAndSetFuncPtrs
	
	//returns a color value (0.0f -> 255.0f) in appropriate location in span of min to max
//	private Float getClrVal(Float val, Float min, Float max) {
//		Float denom = (max-min);
//		if(denom ==0) return 255.0f;
//		return 55.0f + 200.0f * (val - min)/denom;	
//	}
	
	private int getClrValInt(Float val, Float min, Float max) {
		Float denom = (max-min);
		if(denom <=0) return 255;
		return (int) (55.0f + 200.0f * (val - min)/denom);	
	}
	

	
	/**
	 * call whenever grid dimensions change
	 * @param _gridCount # of grid cells per dim
	 * @param _h cell size
	 * @param _numParts # of material points to draw
	 */
	public final void setGridValsAndInit(int _gridCount, float _cellSize, int _numParts, float _partMass) {
		//# of grid cells per side of cube
		gridCount = _gridCount;
		//# of particles
		numParts = _numParts;
		//cell size
		cellSize = _cellSize;	
		
		particleMass = _partMass;		
		maxSimBnds = (gridCount*cellSize)/2.0f;
		minSimBnds = -maxSimBnds;
		gridDim = maxSimBnds - minSimBnds;		
		//scale amount to fill 1500 x 1500 x 1500 visualization cube
		sclAmt = 1500.0f/(gridCount * cellSize);

		setSimFlags(gridIsBuiltIDX, false);
		setSimFlags(simIsBuiltIDX, false);
		resetSim(true);
	}//setGridValsAndInit

	
	//allocate dev mem for all objects based on number of particles
	private void initCUDAMemPtrs_Parts(TreeMap<String, ArrayList<Float[]>> partVals) {
        float h_part_mass[] =new float[numParts];
        float h_part_eye[] =new float[numParts];
        //making class variables so can be rendered
        h_part_pos = new float[3][];
        h_part_vel = new float[3][];
        for(int i=0;i<h_part_pos.length;++i) {
        	h_part_pos[i] = new float[numParts];
        	h_part_vel[i] = new float[numParts];
        }
       
        //h_part_clr = new float[numParts][3];
        h_part_clr_int = new int[numParts][3];
        
		Float[] minVals = partVals.get("minMaxVals").get(0);
		Float[] maxVals = partVals.get("minMaxVals").get(1);       
        Float[] posAra, velAra;
        
        
        for(int i = 0; i < numParts; ++i){
        	h_part_mass[i] = particleMass;
        	posAra = partVals.get("pos").get(i);
        	velAra = partVals.get("vel").get(i);
        	
        	for(int j=0;j<h_part_pos.length;++j) {
        		h_part_pos[j][i]=posAra[j];
        		h_part_vel[j][i]=velAra[j];
        	}
        	int[] clrAra = new int[3];
        	for(int j=0;j<clrAra.length;++j) {
        		clrAra[j] = getClrValInt(h_part_pos[j][i],minVals[j],maxVals[j]);
        	}
        	h_part_clr_int[i] = clrAra;
        	h_part_eye[i]=1.0f;

        }
        cuMemAlloc(part_mass, numPartsFloatSz);  
        cuMemcpyHtoD(part_mass, 	Pointer.to(h_part_mass), numPartsFloatSz);
        cuMemAlloc(part_vol, numPartsFloatSz); 
        cuMemsetD32(part_vol, 0, numParts);	//part_vol is a calculated quantity

        for(int i=0;i<part_pos.length;++i) {
           	cuMemAlloc(part_pos[i], numPartsFloatSz); 
           	cuMemAlloc(part_vel[i], numPartsFloatSz); 

            cuMemcpyHtoD(part_pos[i],	Pointer.to(h_part_pos[i]), numPartsFloatSz);
            cuMemcpyHtoD(part_vel[i],	Pointer.to(h_part_vel[i]), numPartsFloatSz);
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
    }
	
	//allocate dev mem for all objects based on number of grid cells
	private void initCUDAMemPtrs_Grids() {   		
		h_grid_pos = new float[3][];
		h_grid_vel = new float[3][];
		h_grid_accel = new float[3][];
		for(int i=0;i<h_part_pos.length;++i) {
			h_grid_pos[i] = new float[gridSize];
			h_grid_vel[i] = new float[gridSize];
	       	h_grid_accel[i] = new float[gridSize];
        }
		
		h_grid_mass = new float[gridSize];		
		//build grid locations
		int gridDim=0;
		for(int i=0;i<gridCount;++i) {
			float xPos = (i+.5f)*cellSize;
			for(int j=0;j<gridCount;++j) {		
				float yPos = (j+.5f)*cellSize;
				for(int k=0;k<gridCount;++k) {
					float zPos = (k+.5f)*cellSize;
					//weird orientation stuffs
					h_grid_pos[0][gridDim] = zPos;
					h_grid_pos[1][gridDim] = yPos;
					h_grid_pos[2][gridDim] = xPos;					
					++gridDim;
				}
			}
		}				
		
        cuMemAlloc(grid_mass, numGridFloatSz); 
		cuMemsetD32(grid_mass, 0, gridSize);

		for(int i=0;i<grid_vel.length;++i) {
            cuMemAlloc(grid_vel[i], numGridFloatSz);  
            cuMemAlloc(grid_newvel[i], numGridFloatSz);
            cuMemAlloc(grid_force[i], numGridFloatSz);
            
            cuMemsetD32(grid_vel[i], 0, gridSize);  
            cuMemsetD32(grid_newvel[i], 0, gridSize);
            cuMemsetD32(grid_force[i], 0, gridSize);           	
        }
	}//initCUDAMemPtrs_Grids	
	
	private static final double lcl_third = 1.0/3.0;
	//return a float array of random positions within a sphere of radius rad at ctr 
	private Float[] getRandPosInSphereAra(float rad, float[] ctr){
		myVectorf pos = new myVectorf();
		double u = ThreadLocalRandom.current().nextDouble(0,1);
		Float r = (float) (rad * Math.pow(u, lcl_third));
		do{
			pos.set(ThreadLocalRandom.current().nextDouble(-1,1), ThreadLocalRandom.current().nextDouble(-1,1),ThreadLocalRandom.current().nextDouble(-1,1));
		} while (pos.sqMagn > 1.0f);
		Float[] posRes = new Float[] {pos.x, pos.y, pos.z};
		for(int i=0;i<3;++i) {
			posRes[i] *= r;
			posRes[i] += ctr[i];
		}
		return posRes;
	}//getRandPosInSphereAra
	
	/**
	 * create a sphere with given center, with passed # of particles -0 returns ball radius
	 * @param partVals map holding all appropriate particle values
	 * @param ballRad desired radius of sphere
	 * @param numParts
	 * @param initVel
	 * @param ctr
	 * @return
	 */
	protected final float createSphere(TreeMap<String, ArrayList<Float[]>> partVals, float ballRad, int numParts, float[] initVel, float[] ctr) {
    
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "createSphere","Start creating a sphere with "+numParts+" particles in a sphere of radius " + ballRad +".");
		 
		Float[] minVals = partVals.get("minMaxVals").get(0);
		Float[] maxVals = partVals.get("minMaxVals").get(1); 
		
		for (int i=0;i<numParts;++i) {
			Float[] posVals = getRandPosInSphereAra(ballRad, ctr);        				
			for (int v = 0; v < 3; ++v) {
				minVals[v] = (posVals[v] < minVals[v] ?posVals[v] : minVals[v] );
				maxVals[v] = (posVals[v] > maxVals[v] ?posVals[v] : maxVals[v] );
			}
			partVals.get("pos").add(posVals);
			//init vel
			partVals.get("vel").add(new Float[] {initVel[0],initVel[1],initVel[2]});
        }
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "createSphere","Finished creating a sphere with "+numParts+" particles.");
		
		return ballRad;
	}//createSphere

	/**
	 * build initial layout for particles for this simulation
	 * @param partVals initialized string-keyed map holding arrays of values
	 * @return partVals array with sim-appropriate values
	 */
	protected abstract void buildPartLayoutMap(TreeMap<String, ArrayList<Float[]>> partVals);		
	
	/**
	 * initialize all particle values
	 */
	private TreeMap<String, ArrayList<Float[]>> initAllParticles() {
		//Build the particle layout for this simulation
        //create base particle layout
        TreeMap<String, ArrayList<Float[]>> partVals = new TreeMap<String, ArrayList<Float[]>>();
        partVals.put("pos",new ArrayList<Float[]>());
        partVals.put("vel",new ArrayList<Float[]>());
        partVals.put("minMaxVals",new ArrayList<Float[]>());
        //initialize min and max values
        partVals.get("minMaxVals").add(new Float[] {100000.0f,100000.0f,100000.0f});
        partVals.get("minMaxVals").add(new Float[] {-100000.0f,-100000.0f,-100000.0f});      
		
		//determine sim-specific particle layouts
        buildPartLayoutMap(partVals);
        return partVals;
	}//initAllParticles
	
	private void cudaSetup() {   
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Start CUDA Init");
		if (!getSimFlags(CUDADevInit)) {
            //init cuda device and kernel file if not done already - only do 1 time
            initCUDAModuleSetup();
        }
		//total grid size       
        gridSize=gridCount*gridCount*gridCount;
        //initialize all particle values
        TreeMap<String, ArrayList<Float[]>> partVals = initAllParticles();

        numParts = partVals.get("pos").size();
        numPartsFloatSz = numParts * Sizeof.FLOAT;

        win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Total # of particles in simulation : " + numParts);

        //init ptrs to particle-based arrays - numparts and numPartsFloatSz need to be initialized here
        initCUDAMemPtrs_Parts(partVals);		

        //init grid ptrs        
        numGridFloatSz = gridSize * Sizeof.FLOAT;
        initCUDAMemPtrs_Grids();
        
        float delT = getDeltaT();        
        numBlocksParticles=numParts/numCUDAThreads+1;
        numBlocksGrid=gridSize/numCUDAThreads+1;
        
        kernelParams = new TreeMap<String, Pointer>();
        funcGridDimAndMemSize = new HashMap<String, int[]>();
        
        kernelParams.put("projectToGridandComputeForces",Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
        		Pointer.to(new float[] {mat.getLambda0()}), Pointer.to(new float[] {mat.getMu0()}), Pointer.to(new float[] {mat.hardeningCoeff}),
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
				Pointer.to(new float[] {mat.alphaPicFlip}),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(part_vel[0]),Pointer.to(part_vel[1]),Pointer.to(part_vel[2]),
				Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),		
				Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2])));  
		funcGridDimAndMemSize.put("updPartVelocities", new int[] {numBlocksParticles*4, numCUDAThreads*6*Sizeof.FLOAT});
		
	    kernelParams.put("compGridVelocities", Pointer.to(
				Pointer.to(new int[] {gridSize}), Pointer.to(new float[] {gravity[0]}),Pointer.to(new float[] {gravity[1]}),Pointer.to(new float[] {gravity[2]}),Pointer.to(new float[] {delT}),
				Pointer.to(grid_mass),
		
				Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),			
				Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2]),
				Pointer.to(grid_force[0]),Pointer.to(grid_force[1]),Pointer.to(grid_force[2])));
		funcGridDimAndMemSize.put("compGridVelocities", new int[] {numBlocksGrid, 0});
	    
	    kernelParams.put("clearGrid", Pointer.to(
	    		Pointer.to(new int[] {gridSize}), Pointer.to(grid_mass),		
	    		Pointer.to(grid_vel[0]),Pointer.to(grid_vel[1]),Pointer.to(grid_vel[2]),			
	    		Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2]),
				Pointer.to(grid_force[0]),Pointer.to(grid_force[1]),Pointer.to(grid_force[2])));
		funcGridDimAndMemSize.put("clearGrid", new int[] {numBlocksGrid, 0});

	    kernelParams.put("gridCollisions", Pointer.to(
        		Pointer.to(new int[] {gridSize}),Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {cellSize}), 
        		Pointer.to(new float[] {minSimBnds}),Pointer.to(new float[] {maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {delT}),Pointer.to(grid_mass),
        		Pointer.to(grid_newvel[0]),Pointer.to(grid_newvel[1]),Pointer.to(grid_newvel[2])));
		funcGridDimAndMemSize.put("gridCollisions", new int[] {numBlocksGrid, 0});

        kernelParams.put("updDeformationGradient", Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {delT}), Pointer.to(new float[] {cellSize}), Pointer.to(new float[] {minSimBnds}),
        		Pointer.to(new float[] {mat.criticalCompression}), Pointer.to(new float[] {mat.criticalStretch}),
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
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {delT}),
				Pointer.to(part_pos[0]),Pointer.to(part_pos[1]),Pointer.to(part_pos[2]),
				Pointer.to(part_vel[0]),Pointer.to(part_vel[1]),Pointer.to(part_vel[2])				
        		));   
        funcGridDimAndMemSize.put("partCollAndUpdPos", new int[] {numBlocksParticles, 0});        
   
        win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Finished CUDA Init | Launch first MPM Pass.");
	 	//launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {        	launchKernel(initStepFuncKeys[j]);	}
    	setSimFlags(gridIsBuiltIDX, true);
    	setSimFlags(simIsBuiltIDX, true);
    	win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Finished first MPM Pass.");
	}//cudaSetup
	
	/**
	 * launch the kernel specified by the string key
	 * @param kernelVal values for specific kernal 
	 * @param key string key of kernel to launch
	 */
	private void launchKernel(String key) {
		int[] kernelVal = funcGridDimAndMemSize.get(key);
		//dispMessage("MPM_Abs_CUDASim : "+simName, "launchKernel","Launching kernel : " + cuFuncs.get(key) + " with kernelVal[0] : " +kernelVal[0] + " kernelVal[1] : " + kernelVal[1]+" kernelParams.get(key) : "+  kernelParams.get(key) + " key : "+key , MsgCodes.info1);
		
        cuLaunchKernel(cuFuncs.get(key), 
        		kernelVal[0], 1, 1,           // Grid dimension 
                numCUDAThreads, 1, 1,  // Block dimension
                kernelVal[1], null,           // Shared memory size and stream 
                kernelParams.get(key), null);// Kernel- and extra parameters		
	}
    
	//compiles Ptx file from file in passed file name -> cuFileName needs to have format "xxxxx.cu"
	public final void compilePtxFile(String krnFileName, String ptxFileName) throws IOException {
		File cuFile = new File(krnFileName);
		if (!cuFile.exists()) {
			throw new IOException("Kernel file not found: " + krnFileName);
		}
		String modelString = "-m" + System.getProperty("sun.arch.data.model");
		//build compilation command
		String command = "nvcc " + modelString + " -ptx " + cuFile.getPath() + " -o " + ptxFileName;
		//execute compilation
		System.out.println("Executing\n" + command);
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
			System.out.println("nvcc process error : exitValue : " + exitValue);
			System.out.println("errorMessage :\n" + errorMessage);
			System.out.println("outputMessage :\n" + outputMessage);
			throw new IOException("Could not create .ptx file: " + errorMessage);
		}
		System.out.println("Finished compiling PTX file : "+ ptxFileName);
	}

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

	public final void resetSim(boolean freeCuda) {
		if(freeCuda) {
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
		}
		//sim start time - time from when sim object was first instanced
		simStartTime = getCurTime();	
		
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "resetSim","Start resetting sim");
		// Only send original UI values to this sim if it has finished being constructed
//		if ((resetAllVals) && (win.currSim != null)) {
//			win.resetUIVals(false);
//		}
		setSimFlags(CUDADevInit,false);
		setSimFlags(simIsBuiltIDX, false);
		//rebuild cuda kernel configurations
		cudaSetup();

		setSimFlags(simIsBuiltIDX, true);
		win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "resetSim","Finished resetting sim");
	}//initOnce


	//only called after particles are built
	public void setSimPartsAreBuilt() {setSimFlags(simIsBuiltIDX, true);}
	public void setSimGridIsBuilt() {setSimFlags(gridIsBuiltIDX, true);}	
	
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
			case gridIsBuiltIDX	 		: {break;}		
			case showCollider 			: {
				break;}
			case showParticles			: {				
				break;}
			case showParticleVelArrows 	: {
				break;}
			case showGrid				: {
				break;}				
			case showGridVelArrows 		: {
				break;}		
			case showGridAccelArrows 	: {
				break;}
			case showGridMass  			: {
				break;}		
			case showActiveNodes  : {
				break;}
			default :{}			
		}			
	}//setExecFlags
	
	/**
	 * Execute simStepsPerFrame step of simulation
	 * @param modAmtMillis is in milliseconds, counting # of ms since last sim call 
	 * @return true when simulation run is complete - true turns run sim flag off
	 */	
	public final boolean simMe(float modAmtMillis) {
	 	JCudaDriver.cuCtxSetCurrent(pctx);
 		for (int i=0;i<simStepsPerFrame;++i) { 	        
 			//for every sim step, launch each kernel by key specified in simStepFuncKeys
			for (int j=0;j<simStepFuncKeys.length;++j) {	        	launchKernel(simStepFuncKeys[j]);		}
 		}
 		
 		if((getSimFlags(showParticles)) || getSimFlags(showParticleVelArrows)) {
			//copy from device data to host particle position arrays
 			for(int i=0;i<h_part_pos.length;++i) {			cuMemcpyDtoH(Pointer.to(h_part_pos[i]),part_pos[i], numPartsFloatSz);			}
 		} 		
		
 		if(getSimFlags(showParticleVelArrows)) {
 			for(int i=0;i<h_part_vel.length;++i) { 				cuMemcpyDtoH(Pointer.to(h_part_vel[i]),part_vel[i], numPartsFloatSz);			}
 		}
 		
		if(getSimFlags(showGridVelArrows)) {
 			for(int i=0;i<h_grid_vel.length;++i) { 				cuMemcpyDtoH(Pointer.to(h_grid_vel[i]),grid_newvel[i], numGridFloatSz);			}
		}
		if(getSimFlags(showGridAccelArrows)) {
 			for(int i=0;i<h_grid_accel.length;++i) { 				cuMemcpyDtoH(Pointer.to(h_grid_accel[i]),grid_force[i], numGridFloatSz);			}		
		}
		
		if(getSimFlags(showGridMass)) {
			cuMemcpyDtoH(Pointer.to(h_grid_mass),grid_mass, numGridFloatSz);	
		}

		return false;
	}//simMe
	
	//sim method to show execution time for each step
	public abstract boolean simMeDebug(float modAmtMillis);	//simMeDebug	

	public static float getDeltaT() {return base_MPMCudaSim.deltaT;}
	public void setDeltaT(float _delT) {base_MPMCudaSim.deltaT = _delT;}

	public void showTimeMsgSimStart(String _str) {win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "showTimeMsgNow",_str+" Time Now : "+(getCurTime() - simStartTime));}
	//display message and time now
	public void showTimeMsgNow(String _str, long stTime) {	win.getMsgObj().dispInfoMessage("MPM_Abs_CUDASim : "+simName, "showTimeMsgNow",_str+" Time Now : "+(getCurTime() - stTime)+" ms");}
	
	/////////////////////////////
	// utility : time stamp; display messages; state flags
	
	//get time from "start time" (ctor run for map manager)
	protected long getCurTime() {			
		Instant instant = Instant.now();
		return instant.toEpochMilli() - expMgrBuiltTime;//milliseconds since 1/1/1970, subtracting when mapmgr was built to keep millis low		
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

	//draw 1 frame of results	//animTimeMod is in seconds, counting # of seconds since last draw
	public final void drawMe(float animTimeMod) {
		if(!getSimFlags(simIsBuiltIDX)) {return;}//if not built yet, don't try to draw anything
		//render all particles - TODO determine better rendering method
		pa.pushMatState();
		//set stroke values and visual scale
			pa.setStrokeWt(2.0f/sclAmt);
			pa.scale(sclAmt);	
			
			//TODO control via UI - draw ever pincr points in point array
			int pincr = 1;
			//draw material points
			if(getSimFlags(showParticles)){	
				pa.pushMatState();
				pa.drawPointCloudWithColors(numParts, pincr, h_part_clr_int, h_part_pos[0], h_part_pos[1], h_part_pos[2]);
				pa.popMatState();
			} 
			
			if(getSimFlags(showParticleVelArrows)){
				drawPointVel(pincr);
			}
			
			//draw colliders, if exist
			drawCollider(animTimeMod);
			
			//if desired, draw grid
			if(getSimFlags(showGrid)) {	drawGrid();}
			//TODO control via UI - scale size of vectors
			float mult = .01f;
			if (getSimFlags(showGridVelArrows)) {	_drawGridVec(mult, gridVecClr, h_grid_vel[0], h_grid_vel[1], h_grid_vel[2]);}
			if (getSimFlags(showGridAccelArrows)){	_drawGridVec(mult, gridAccelClr, h_grid_accel[0], h_grid_accel[1], h_grid_accel[2]);}
			if(getSimFlags(showGridMass)) {			_drawGridScalar(mult, gridMassClr, h_grid_mass);}
		pa.popMatState();
	}//drawMe
	
	private void drawPointVel(int pincr) {
		pa.pushMatState();
		float mult = .01f;
		float minMag = MyMathUtils.EPS_F/mult;
		for(int i=0;i<=numParts-pincr;i+=pincr) {					
			if(		(Math.abs(h_part_vel[0][i]) > minMag) || 
					(Math.abs(h_part_vel[1][i]) > minMag) || 
					(Math.abs(h_part_vel[2][i]) > minMag)) {
				pa.pushMatState();
				pa.setStroke(h_part_clr_int[i], 100);
				pa.translate(h_part_pos[0][i], h_part_pos[1][i], h_part_pos[2][i]);
				pa.drawLine(0,0,0, mult*h_part_vel[0][i],mult*h_part_vel[1][i],mult*h_part_vel[2][i]);
				pa.popMatState();
			}
		}			
		pa.popMatState();
	}//drawPointVel

	private void drawGrid() {
		int incr = 10;
		pa.pushMatState();		
			pa.setStroke(0,0,0,20);
			pa.translate(minSimBnds,minSimBnds,minSimBnds);
			//shows every "incr" gridcells
			for (int i=0; i<=gridCount;i+=incr) {
				float iLoc = i*cellSize;
				for(int j=0;j<=gridCount;j+=incr) {
					myVectorf startPos=new myVectorf(iLoc,j*cellSize,0.0f);
					myVectorf endPos=new myVectorf(iLoc, startPos.y,gridDim);
					pa.drawLine(startPos,endPos);
				}
				for(int k=0;k<=gridCount;k+=incr) {
					myVectorf startPos=new myVectorf(iLoc,0.0f, k*cellSize);
					myVectorf endPos=new myVectorf(iLoc,gridDim,startPos.z);
					pa.drawLine(startPos,endPos);
				}
			}
			for(int j=0;j<=gridCount;j+=incr) {
				float jLoc = j*cellSize;
				for(int k=0;k<=gridCount;k+=incr) {
					myVectorf startPos=new myVectorf(0.0f,jLoc,k*cellSize);
					myVectorf endPos=new myVectorf(gridDim,jLoc,startPos.z);
					pa.drawLine(startPos,endPos);
				}
			}
		pa.popMatState();		
	}//drawGrid()
	
	private void _drawGridVec(float mult, int[] clr, float[] xVal, float[] yVal, float[] zVal) {
		float minMag = MyMathUtils.EPS_F/mult;
		pa.pushMatState();	
		pa.setStroke(clr,20);
		pa.translate(minSimBnds,minSimBnds,minSimBnds);
		for (int i=0; i<gridSize;++i) {			
			if(		(Math.abs(xVal[i]) > minMag) || 
					(Math.abs(yVal[i]) > minMag) || 
					(Math.abs(zVal[i]) > minMag)) {
				pa.pushMatState();	
				pa.translate(h_grid_pos[0][i], h_grid_pos[1][i], h_grid_pos[2][i]);
				pa.drawLine(0,0,0, mult*xVal[i],mult*yVal[i],mult*zVal[i]);
				pa.popMatState();
			}
		}
		pa.popMatState();
	}
	
	private void _drawGridScalar(float mult, int[] clr, float[] xVal) {
		float minMag = MyMathUtils.EPS_F/(.01f * mult);
		pa.pushMatState();	
		pa.setSphereDetail(4);
		pa.setStroke(clr,20);
		pa.translate(minSimBnds,minSimBnds,minSimBnds);
		for (int i=0; i<gridSize;++i) {			
			if(		(Math.abs(xVal[i]) > minMag)) {
				pa.pushMatState();	
				pa.translate(h_grid_pos[0][i], h_grid_pos[1][i], h_grid_pos[2][i]);
				pa.drawSphere(xVal[i]*mult);
				pa.popMatState();
			}
		}
		pa.popMatState();		
	}
	
	//draw internal-to-sim colliders, if they exist
	protected abstract void drawCollider(float animTimeMod);
	
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




