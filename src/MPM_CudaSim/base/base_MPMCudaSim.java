package MPM_CudaSim.base;

import static jcuda.driver.JCudaDriver.*;

import java.io.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;

//import org.jblas.FloatMatrix;

import MPM_CudaSim.myMaterial;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;
import base_Utils_Objects.io.ConsoleCLR;
import base_Utils_Objects.io.MsgCodes;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import jcuda.*;
import jcuda.driver.*;
import jcuda.runtime.JCuda;

//abstract class describing a simulation world.  called by sim window to executed simulation and to render results
//instancing classes can hold different configurations/initialization setups
public abstract class base_MPMCudaSim{	
	public static IRenderInterface pa;
	//name of instancing sim
	public final String simName;
	
	//cuda kernel file name
	private String ptxFileName = "MPM_CUDA_Sim_New.ptx";	
	//vs c++ compiler location
	//private final static String VSCppCompLoc = "C:\\Program Files (x86)\\Microsoft Visual Studio 14.0\\VC\\bin\\x86_amd64\\";
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
	protected  static float deltaT = 4e-4f;


	//const matrix for calculations - z is down
	//protected final FloatMatrix gravity = new FloatMatrix(new float[] {0, 0, -9.8f});
	protected final float[] gravity = new float[] {0, 0, -9.8f};
	//scale amount for visualization to fill cube frame in 3d world; particle radius, scaled by different visual scales
	protected float sclAmt;
	
	//material quantities of particle matter
	public myMaterial mat;	
	
	//friction coefficients of colliders
	public static float wallFric = 1.0f, collFric = 1.0f;

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

	//iterations per frame
	public static int simStepsPerFrame = 2;

	//time of current process start, from initial construction of mapmgr - TODO use this to monitor specific process time elapsed.  set to 0 at beginning of a particular process, then measure time elapsed in process
	protected long simStartTime;
	//time mapMgr built, in millis - used as offset for instant to provide smaller values for timestamp
	protected final long expMgrBuiltTime;	
	//constants for collider calcs
	protected static float cyl_da = (float) (Math.PI/18.0f), TWO_PI = (float) (2.0 * Math.PI);	
	
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
    //float[][] h_part_clr;
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
	public base_MPMCudaSim(IRenderInterface _pa,String _simName, int _gridCount, float _h, int _numParts, float _density) {		
		pa=_pa;simName = _simName;	
		
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
	//run 1 time to load kernel and assign function pointers to functions
	private void initOnceCUDASetup() {
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
			System.out.println("Registering Kernel Function Key : " + key);
			CUfunction c = new CUfunction();			
			cuModuleGetFunction(c, module, key);
			cuFuncs.put(key,  c);
		}
    
        setSimFlags(CUDADevInit, true);
	}//loadModuleAndSetFuncPtrs
	
	//returns a color value (0.0f -> 255.0f) in appropriate location in span of min to max
	private Float getClrVal(Float val, Float min, Float max) {
		Float denom = (max-min);
		if(denom ==0) return 255.0f;
		return 55.0f + 200.0f * (val - min)/denom;	
	}
	
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
	public final void setGridValsAndInit(int _gridCount, float _h, int _numParts, float _partMass) {
		//# of grid cells per side of cube
		gridCount = _gridCount;
		//# of particles
		numParts = _numParts;
		//cell size
		cellSize = _h;	
		
		particleMass = _partMass;		
		maxSimBnds = (gridCount*cellSize)/2.0f;
		minSimBnds = -maxSimBnds;
		gridDim = maxSimBnds - minSimBnds;		
		//scale amount to fill 1500 x 1500 x 1500 visualization cube
		sclAmt = 1500.0f/(gridCount * cellSize);

		//notify gridbuilder that we have a new grid
		//gridThdMgr.setNewSimVals();
		setSimFlags(gridIsBuiltIDX, false);
		//partThdMgr.setNewSimVals();
		setSimFlags(simIsBuiltIDX, false);
		//build grid - upon completion the gridbuilder will call resetSim to build particles
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
	
	//create a sphere with given center, with passed # of particles -0 returns ball radius
	protected final float createSphere(TreeMap<String, ArrayList<Float[]>> partVals, int numParts, float[] initVel, float[] ctr) {
    
		//build sphere of particles - scale volume of sphere based on cuberoot of # of particles, with 1000 particles being baseline sphere - radius will be function of how many particles are built
        float ballRad = (float) (3.0*Math.cbrt(numParts)/sclAmt);		
		dispMessage("MPM_Abs_CUDASim : "+simName, "createSphere","Start creating a sphere with "+numParts+" particles in a sphere of radius " + ballRad +".", MsgCodes.info1);
		 
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
		dispMessage("MPM_Abs_CUDASim : "+simName, "createSphere","Finished creating a sphere with "+numParts+" particles.", MsgCodes.info1);
		
		return ballRad;
	}//createSphere

	protected abstract TreeMap<String, ArrayList<Float[]>> buildPartLayout();		
	
	private void cudaSetup() {//throws IOException {     
		dispMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Start CUDA Init", MsgCodes.info1);
		if (!getSimFlags(CUDADevInit)) {
            //init cuda device and kernel file if not done already - only do 1 time
            initOnceCUDASetup();
        }
		//Build the particle layout for this simulation
		TreeMap<String, ArrayList<Float[]>> partVals = buildPartLayout();
		//total grid size       
        gridSize=gridCount*gridCount*gridCount;
        
        numParts = partVals.get("pos").size();
        System.out.println("# of particles : " + numParts);
        numPartsFloatSz = numParts * Sizeof.FLOAT;

        //init ptrs to particle-based arrays - numparts and numPartsFloatSz need to be initialized here
        initCUDAMemPtrs_Parts(partVals);

        //init grid ptrs
        float delT = getDeltaT();

        numGridFloatSz = gridSize * Sizeof.FLOAT;
        initCUDAMemPtrs_Grids();
        
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
   
        dispMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Finished CUDA Init | Launch first MPM Pass.", MsgCodes.info1);
	 	//launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {        	launchKernel(initStepFuncKeys[j]);	}
    	setSimFlags(gridIsBuiltIDX, true);
    	setSimFlags(simIsBuiltIDX, true);
    	dispMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Finished first MPM Pass.", MsgCodes.info1);
	}//cudaSetup
	
	/**
	 * launch the kernel specified by the string key
	 * @param kernelVal values for specific kernal 
	 * @param key string key of kernel to launch
	 */
	private void launchKernel(String key) {
		int[] kernelVal = funcGridDimAndMemSize.get(key);
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


	//call whenever setting/resetting simulation world - no changes to particle count/grid dimensions, just resetting to initial state
	//does not rebuild thread worker list, just re-executes existing thread workers
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
		dispMessage("MPM_Abs_CUDASim : "+simName, "resetSim","Start resetting sim", MsgCodes.info1);
		setSimFlags(CUDADevInit,false);
		setSimFlags(simIsBuiltIDX, false);
		
		cudaSetup();

		setSimFlags(simIsBuiltIDX, true);
		dispMessage("MPM_Abs_CUDASim : "+simName, "resetSim","Finished resetting sim", MsgCodes.info1);
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
//		
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
		//Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z)
		return false;
	}//simMe}
	
	//sim method to show execution time for each step
	public abstract boolean simMeDebug(float modAmtMillis);	//simMeDebug	

	public static float getDeltaT() {return base_MPMCudaSim.deltaT;}
	public void setDeltaT(float _delT) {base_MPMCudaSim.deltaT = _delT;}

	public void showTimeMsgSimStart(String _str) {dispMessage("MPM_Abs_CUDASim : "+simName, "showTimeMsgNow",_str+" Time Now : "+(getCurTime() - simStartTime), MsgCodes.info1);}
	//display message and time now
	public void showTimeMsgNow(String _str, long stTime) {	dispMessage("MPM_Abs_CUDASim : "+simName, "showTimeMsgNow",_str+" Time Now : "+(getCurTime() - stTime)+" ms", MsgCodes.info1);}
	
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
	// start message display functionality
	
	private String buildClrStr(ConsoleCLR bk, ConsoleCLR clr, String str) {return bk.toString() + clr.toString() + str + ConsoleCLR.RESET.toString();	}
	private String _processMsgCode(String src, MsgCodes useCode) {
		if (!supportsANSITerm) {return src;}
		switch(useCode) {//add background + letter color for messages
			//info messages
			case info1 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.WHITE, src);}				//basic informational printout
			case info2 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.CYAN, src);}
			case info3 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.YELLOW, src);}				//informational output from external source
			case info4 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.GREEN, src);}
			case info5 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.CYAN_BOLD, src);}			//beginning or ending of processing
			//warning messages                                                 , 
			case warning1 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLACK_BOLD, src);}
			case warning2 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLUE_BOLD, src);}			//warning info re: ui does not exist
			case warning3 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLACK_UNDERLINED, src);}
			case warning4 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLUE_UNDERLINED, src);}		//info message about unexpected behavior
			case warning5 : {	return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND, ConsoleCLR.BLUE_BRIGHT, src);}
			//error messages                                                   , 
			case error1 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.RED_UNDERLINED, src);}		//try/catch error
			case error2 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.RED_BOLD, src);}			//code-based error
			case error3 : {		return  buildClrStr(ConsoleCLR.RED_BACKGROUND_BRIGHT, ConsoleCLR.BLACK_BOLD, src);}		//file load error
			case error4 : {		return  buildClrStr(ConsoleCLR.WHITE_BACKGROUND_BRIGHT, ConsoleCLR.RED_BRIGHT, src);}	//error message thrown by external process
			case error5 : {		return  buildClrStr(ConsoleCLR.BLACK_BACKGROUND, ConsoleCLR.RED_BOLD_BRIGHT, src);}
		}
		return src;
	}//_processMsgCode	

	public void dispMessageAra(String[] _sAra, String _callingClass, String _callingMethod, int _perLine, MsgCodes useCode) {dispMessageAra( _sAra,  _callingClass, _callingMethod, _perLine,  useCode, true);}
	//show array of strings, either just to console or to applet window
	public void dispMessageAra(String[] _sAra, String _callingClass, String _callingMethod, int _perLine, MsgCodes useCode, boolean onlyConsole) {
		String callingClassPrfx = getTimeStrFromProcStart() +"|" + _callingClass;		 
		for(int i=0;i<_sAra.length; i+=_perLine){
			String s = "";
			for(int j=0; j<_perLine; ++j){	
				if((i+j >= _sAra.length)) {continue;}
				s+= _sAra[i+j]+ "\t";}
			_dispMessage_base(callingClassPrfx,_callingMethod,s, useCode,onlyConsole);
		}
	}//dispMessageAra

	public void dispMessage(String srcClass, String srcMethod, String msgText, MsgCodes useCode){_dispMessage_base(getTimeStrFromProcStart() +"|" + srcClass,srcMethod,msgText, useCode,true);	}	
	public void dispMessage(String srcClass, String srcMethod, String msgText, MsgCodes useCode, boolean onlyConsole) {_dispMessage_base(getTimeStrFromProcStart() +"|" + srcClass,srcMethod,msgText, useCode,onlyConsole);	}	
	private void _dispMessage_base(String srcClass, String srcMethod, String msgText, MsgCodes useCode, boolean onlyConsole) {		
		String msg = _processMsgCode(srcClass + "::" + srcMethod + " : " + msgText, useCode);
		if((onlyConsole) || (pa == null)) {		System.out.println(msg);	} else {		pa.outStr2Scr(msg);	}
	}//dispMessage
	
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
			
			//draw material points
			if(getSimFlags(showParticles)){drawPoints(true,getSimFlags(showParticleVelArrows));}
			else if(getSimFlags(showParticleVelArrows)){drawPoints(false,true);}
			
			//draw colliders, if exist
			drawCollider(animTimeMod);
			
			//if desired, draw grid
			if(getSimFlags(showGrid)) {							drawGrid(true, getSimFlags(showGridVelArrows),getSimFlags(showGridAccelArrows));	}
			else if (getSimFlags(showGridVelArrows)) {			drawGrid(false, true,getSimFlags(showGridAccelArrows));	}
			else if (getSimFlags(showGridAccelArrows)){			drawGrid(false, false,true);	}
			
			if(getSimFlags(showGridMass)) {
				float mult = 0.1f;
				_drawGridScalar(mult, gridMassClr, h_grid_mass);
				
			}
		pa.popMatState();
	}//drawMe
	
	private void drawPoints(boolean showPoints, boolean showPointVel) {
		pa.pushMatState();
		//draw the points
		int pincr = 1;
		//draw all points as a shape object
		if(showPoints) {pa.drawPointCloudWithColors(numParts, pincr, h_part_clr_int, h_part_pos[0], h_part_pos[1], h_part_pos[2]);}
		
		if(showPointVel) {
			float mult = .01f;
			float minMag = MyMathUtils.eps_f/mult;
			for(int i=0;i<=numParts-pincr;i+=pincr) {					
				if(		(Math.abs(h_part_vel[0][i]) > minMag) || 
						(Math.abs(h_part_vel[1][i]) > minMag) || 
						(Math.abs(h_part_vel[2][i]) > minMag)) {
					pa.pushMatState();
					pa.setStroke(h_part_clr_int[i], 100);
					pa.translate(h_part_pos[0][i], h_part_pos[1][i], h_part_pos[2][i]);
					//((my_procApplet)pa).stroke(h_part_clr_int[i][0], h_part_clr_int[i][1], h_part_clr_int[i][2]);
					pa.drawLine(0,0,0, mult*h_part_vel[0][i],mult*h_part_vel[1][i],mult*h_part_vel[2][i]);
					//((my_procApplet)pa).vertex(h_part_pos_x[i], h_part_pos_y[i], h_part_pos_z[i]);
					pa.popMatState();
				}
			}			
		}
		pa.popMatState();
	}

	private void drawGrid(boolean showGrid, boolean showGridVel, boolean showGridForce) {
		int incr = 10;
		pa.pushMatState();		
		if(showGrid) {	
			pa.pushMatState();	
			pa.setStroke(0,0,0,20);
			pa.translate(minSimBnds,minSimBnds,minSimBnds);
			//shows every "incr" gridcells
			for (int i=0; i<=gridCount;i+=incr) {
				float iLoc = i*cellSize;
				for(int j=0;j<=gridCount;j+=incr) {
					myVectorf startPos=new myVectorf(iLoc,j*cellSize,0.0f);
					myVectorf endPos=new myVectorf(iLoc,j*cellSize,gridDim);
					pa.drawLine(startPos,endPos);
				}
				for(int k=0;k<=gridCount;k+=incr) {
					myVectorf startPos=new myVectorf(iLoc,0.0f, k*cellSize);
					myVectorf endPos=new myVectorf(iLoc,gridDim,k*cellSize);
					pa.drawLine(startPos,endPos);
				}
			}
			for(int j=0;j<=gridCount;j+=incr) {
				float jLoc = j*cellSize;
				for(int k=0;k<=gridCount;k+=incr) {
					myVectorf startPos=new myVectorf(0.0f,jLoc,k*cellSize);
					myVectorf endPos=new myVectorf(gridDim,jLoc,k*cellSize);
					pa.drawLine(startPos,endPos);
				}
			}
			pa.popMatState();
		}
		if(showGridVel) {
			float mult = .01f;
			_drawGridVec(mult, gridVecClr, h_grid_vel[0], h_grid_vel[1], h_grid_vel[2]);
		}
		
		if(showGridForce) {
			float mult = .01f;
			_drawGridVec(mult, gridAccelClr, h_grid_accel[0], h_grid_accel[1], h_grid_accel[2]);
		}
		
		pa.popMatState();
		
	}
	
	private void _drawGridVec(float mult, int[] clr, float[] xVal, float[] yVal, float[] zVal) {
		float minMag = MyMathUtils.eps_f/mult;
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
		float minMag = MyMathUtils.eps_f/(.01f * mult);
		pa.pushMatState();	
		pa.setSphereDetail(4);
		pa.setStroke(clr,20);
		pa.translate(minSimBnds,minSimBnds,minSimBnds);
		for (int i=0; i<gridSize;++i) {			
			if(		(Math.abs(xVal[i]) > minMag)) {
				pa.pushMatState();	
				//pa.translate(h_grid_pos_x[i], h_grid_pos_y[i], h_grid_pos_z[i]);
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
		for(float a=0; a<=TWO_PI+cyl_da; a+=cyl_da) {
			rcA = (float) (r*Math.cos(a)); 
			rsA = (float) (r*Math.sin(a));
			res[0].add(Pf(P,rcA,I,rsA,J,0.0,V)); 
			res[1].add(Pf(P,rcA,I,rsA,J,1.0,V));
		}
		return res;
	}
	
}//class MPM_ABS_Sim 




