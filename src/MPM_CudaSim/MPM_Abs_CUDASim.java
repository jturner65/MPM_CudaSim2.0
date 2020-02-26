package MPM_CudaSim;

import static jcuda.driver.JCudaDriver.*;

import java.io.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;

import org.jblas.FloatMatrix;

import base_UI_Objects.my_procApplet;
import base_Utils_Objects.io.ConsoleCLR;
import base_Utils_Objects.io.MsgCodes;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import jcuda.*;
import jcuda.driver.*;
import jcuda.runtime.JCuda;
import processing.core.PConstants;

//abstract class describing a simulation world.  called by sim window to executed simulation and to render results
//instancing classes can hold different configurations/initialization setups
public abstract class MPM_Abs_CUDASim{	
	public static my_procApplet pa;
	//name of instancing sim
	public final String simName;
	
	//cuda kernel file name
	private String ptxFileName = "MPM_CUDA_Sim_New.ptx";	
	//vs c++ compiler location
	//private final static String VSCppCompLoc = "C:\\Program Files (x86)\\Microsoft Visual Studio 14.0\\VC\\bin\\x86_amd64\\";
	//these values are just to set initial values for simulation
	public static int numGridCells = 200;
	public static float cellSize = .10f;
	public static int numPartsUI_Init = 100000;
	
	//# of particles to have in sim
	protected int numParts = numPartsUI_Init;
	//grid dim - cube so same in all 3 dims; # of particles in sim
	protected int gridCount = numGridCells;
	protected float h = cellSize;
	//# parts * size of float - currently always builds same # of particles
	protected long numPartsFloatSz;
	
	//simulation boundaries - symmetric cube, only need min and max, grid length per dim
	protected float minSimBnds, maxSimBnds, gridDim;
	//dimension of grid cells
	protected float half_h, h3;
	//whether current OS supports ansi terminal color settings
	public static boolean supportsANSITerm = false;
	//parameters
	//timestep of simulation - 
	private static float deltaT = 4e-4f;
	//snow density varies from 50 to ~800
	//powdery snow is about 100
	//wet firm compacted snow is about 600
	private static float initDensity = 100.0f;

	//const matrix for calculations - z is down
	final FloatMatrix gravity = new FloatMatrix(new float[] {0, 0, -9.8f});

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
		showParticleVelArrows 	= 4,			//plot velocity arrows for each particle	
		showGrid				= 5,			//plot the computational grid
		showGridVelArrows 		= 6,			//plot velocity arrows for each gridNode
		showGridAccelArrows 	= 7,			//plot acceleration arrows for each gridNode
		showGridMass  			= 8,			//plot variable sized spheres proportional to gridnode mass
		showActiveNodes		  	= 9,			//show the grid nodes influenced by each particle
		CUDADevInit				= 10;			//if 1 time cuda device and kernel file init is complete
	protected static final int numSimFlags = 11;

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
	protected String[] initStepFuncKeys = new String[] {"projectToGridInit", "computeVol","compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "clearGrid","partCollAndUpdPos"}; 
	protected String[] simStepFuncKeys = new String[] {"projectToGridandComputeForces", "compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "clearGrid", "partCollAndUpdPos"}; 
   
	protected HashMap<String, int[]> funcGridDimAndMemSize;
    
    CUdeviceptr part_mass = new CUdeviceptr();
    CUdeviceptr part_vol = new CUdeviceptr();
    CUdeviceptr part_pos_x = new CUdeviceptr();
    CUdeviceptr part_pos_y = new CUdeviceptr();
    CUdeviceptr part_pos_z = new CUdeviceptr();    
    CUdeviceptr part_vel_x = new CUdeviceptr();
    CUdeviceptr part_vel_y = new CUdeviceptr();
    CUdeviceptr part_vel_z = new CUdeviceptr();    
    CUdeviceptr part_fe_11 = new CUdeviceptr();
    CUdeviceptr part_fe_12 = new CUdeviceptr();
    CUdeviceptr part_fe_13 = new CUdeviceptr();
    CUdeviceptr part_fe_21 = new CUdeviceptr();
    CUdeviceptr part_fe_22 = new CUdeviceptr();
    CUdeviceptr part_fe_23 = new CUdeviceptr();
    CUdeviceptr part_fe_31 = new CUdeviceptr();
    CUdeviceptr part_fe_32 = new CUdeviceptr();
    CUdeviceptr part_fe_33 = new CUdeviceptr();    
    CUdeviceptr part_fp_11 = new CUdeviceptr();
    CUdeviceptr part_fp_12 = new CUdeviceptr();
    CUdeviceptr part_fp_13 = new CUdeviceptr();
    CUdeviceptr part_fp_21 = new CUdeviceptr();
    CUdeviceptr part_fp_22 = new CUdeviceptr();
    CUdeviceptr part_fp_23 = new CUdeviceptr();
    CUdeviceptr part_fp_31 = new CUdeviceptr();
    CUdeviceptr part_fp_32 = new CUdeviceptr();
    CUdeviceptr part_fp_33 = new CUdeviceptr();    
    
    CUdeviceptr grid_mass = new CUdeviceptr();
    CUdeviceptr grid_pos_x = new CUdeviceptr();
    CUdeviceptr grid_pos_y = new CUdeviceptr();
    CUdeviceptr grid_pos_z = new CUdeviceptr();    
    CUdeviceptr grid_vel_x = new CUdeviceptr();
    CUdeviceptr grid_vel_y = new CUdeviceptr();
    CUdeviceptr grid_vel_z = new CUdeviceptr();    
    CUdeviceptr grid_newvel_x = new CUdeviceptr();
    CUdeviceptr grid_newvel_y = new CUdeviceptr();
    CUdeviceptr grid_newvel_z = new CUdeviceptr();    
    CUdeviceptr grid_force_x = new CUdeviceptr();
    CUdeviceptr grid_force_y = new CUdeviceptr();
    CUdeviceptr grid_force_z = new CUdeviceptr();

    float[] h_part_pos_x,h_part_pos_y, h_part_pos_z;
    //colors based on initial location
    //float[][] h_part_clr;
    int[][] h_part_clr_int;
    
    CUcontext pctx;
    CUdevice dev;
    CUmodule module;
    CUgraphicsResource pCudaResource;
    
    public int gridSize;
	public int numBlocksParticles;
	public int numBlocksGrid;
	//# of cuda threads
	public int numCUDAThreads=128;
	
	//grid count per side - center grid always in display; grid cell dim per side
	//@SuppressWarnings("unchecked")
	public MPM_Abs_CUDASim(my_procApplet _pa,String _simName, int _gridCount, float _h, int _numParts) {		
		pa=_pa;simName = _simName;
//		//for multithreading - do not use instanced version in PApplet - we may not use processing-based build to run simulation
//		th_exec = Executors.newCachedThreadPool();		
//		int numThreadsTtl = Runtime.getRuntime().availableProcessors();
//		//# of threads to use to build structures - doesn't change
//		numThreadsAvail = (numThreadsTtl > 2 ? numThreadsTtl-2 : 1);
		
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
		setGridValsAndInit(_gridCount, _h, _numParts);
		resetSim(false);
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
    
        this.setSimFlags(CUDADevInit, true);
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
	
	//allocate dev mem for all objects based on number of particles
	private void initCUDAMemPtrs_Parts(TreeMap<String, ArrayList<Float[]>> partVals) {
        float h_part_mass[] =new float[numParts];
        float h_part_vel_x[] =new float[numParts];
        float h_part_vel_y[] =new float[numParts];
        float h_part_vel_z[] =new float[numParts];
        float h_part_eye[] =new float[numParts];
        //making class variables so can be rendered
        h_part_pos_x = new float[numParts];
        h_part_pos_y =new float[numParts];
        h_part_pos_z =new float[numParts];
        //h_part_clr = new float[numParts][3];
        h_part_clr_int = new int[numParts][3];
        
		Float[] minVals = partVals.get("minMaxVals").get(0);
		Float[] maxVals = partVals.get("minMaxVals").get(1);       
        Float[] posAra, velAra;
        
        float initMass = h*h*h*initDensity;//myParticle.density;
        for(int i = 0; i < numParts; ++i){
        	h_part_mass[i]=initMass;
        	posAra = partVals.get("pos").get(i);
        	h_part_pos_x[i]=posAra[0];
        	h_part_pos_y[i]=posAra[1];
        	h_part_pos_z[i]=posAra[2];
//        	h_part_clr[i] = new float[] {
//    			getClrVal(h_part_pos_x[i],minVals[0],maxVals[0]), 
//				getClrVal(h_part_pos_y[i],minVals[1],maxVals[1]), 
//				getClrVal(h_part_pos_z[i],minVals[2],maxVals[2]) 
//        	};  
        	h_part_clr_int[i] = new int[] {
        			getClrValInt(h_part_pos_x[i],minVals[0],maxVals[0]), 
        			getClrValInt(h_part_pos_y[i],minVals[1],maxVals[1]), 
        			getClrValInt(h_part_pos_z[i],minVals[2],maxVals[2]) 
        	};  
        	velAra = partVals.get("vel").get(i);
        	//h_part_clr[i] = new float[] {(.5f + h_part_pos_x[i])*128, (.5f + h_part_pos_y[i])*128, (.5f + h_part_pos_z[i])*128};
        	h_part_vel_x[i]=velAra[0];
        	h_part_vel_y[i]=velAra[1];
        	h_part_vel_z[i]=velAra[2];
        	h_part_eye[i]=1.0f;

        }
        cuMemAlloc(part_mass, numPartsFloatSz);  
        cuMemAlloc(part_vol, numPartsFloatSz); 
        cuMemAlloc(part_pos_x, numPartsFloatSz);
        cuMemAlloc(part_pos_y, numPartsFloatSz);
        cuMemAlloc(part_pos_z, numPartsFloatSz); 
        cuMemAlloc(part_vel_x, numPartsFloatSz);
        cuMemAlloc(part_vel_y, numPartsFloatSz); 
        cuMemAlloc(part_vel_z, numPartsFloatSz);        

        cuMemAlloc(part_fe_11, numPartsFloatSz);
        cuMemAlloc(part_fe_12, numPartsFloatSz);
        cuMemAlloc(part_fe_13, numPartsFloatSz);  
        cuMemAlloc(part_fe_21, numPartsFloatSz);  
        cuMemAlloc(part_fe_22, numPartsFloatSz);
        cuMemAlloc(part_fe_23, numPartsFloatSz); 
        cuMemAlloc(part_fe_31, numPartsFloatSz);   
        cuMemAlloc(part_fe_32, numPartsFloatSz);
        cuMemAlloc(part_fe_33, numPartsFloatSz);        
    
        cuMemAlloc(part_fp_11, numPartsFloatSz); 
        cuMemAlloc(part_fp_12, numPartsFloatSz);
        cuMemAlloc(part_fp_13, numPartsFloatSz);   
        cuMemAlloc(part_fp_21, numPartsFloatSz);  
        cuMemAlloc(part_fp_22, numPartsFloatSz);   
        cuMemAlloc(part_fp_23, numPartsFloatSz);   
        cuMemAlloc(part_fp_31, numPartsFloatSz);  
        cuMemAlloc(part_fp_32, numPartsFloatSz);    
        cuMemAlloc(part_fp_33, numPartsFloatSz);
        
        cuMemcpyHtoD(part_mass, 	Pointer.to(h_part_mass), numPartsFloatSz);
        cuMemcpyHtoD(part_pos_x,	Pointer.to(h_part_pos_x), numPartsFloatSz);
        cuMemcpyHtoD(part_pos_y,	Pointer.to(h_part_pos_y), numPartsFloatSz);
        cuMemcpyHtoD(part_pos_z,	Pointer.to(h_part_pos_z), numPartsFloatSz);
        cuMemcpyHtoD(part_vel_x,	Pointer.to(h_part_vel_x), numPartsFloatSz);
        cuMemcpyHtoD(part_vel_y,	Pointer.to(h_part_vel_y), numPartsFloatSz);
        cuMemcpyHtoD(part_vel_z,	Pointer.to(h_part_vel_z), numPartsFloatSz);
      
        cuMemcpyHtoD(part_fe_11,	Pointer.to(h_part_eye), numPartsFloatSz);
        cuMemcpyHtoD(part_fe_22,	Pointer.to(h_part_eye), numPartsFloatSz);
        cuMemcpyHtoD(part_fe_33,	Pointer.to(h_part_eye), numPartsFloatSz);
        
        cuMemcpyHtoD(part_fp_11,	Pointer.to(h_part_eye), numPartsFloatSz);
        cuMemcpyHtoD(part_fp_22,	Pointer.to(h_part_eye), numPartsFloatSz);
        cuMemcpyHtoD(part_fp_33,	Pointer.to(h_part_eye), numPartsFloatSz);
        
        cuMemsetD32(part_fe_12, 0, numParts);
        cuMemsetD32(part_fe_13, 0, numParts);
        cuMemsetD32(part_fe_21, 0, numParts);
        cuMemsetD32(part_fe_23, 0, numParts);
        cuMemsetD32(part_fe_31, 0, numParts);
        cuMemsetD32(part_fe_32, 0, numParts);
        
        cuMemsetD32(part_fp_12, 0, numParts);
        cuMemsetD32(part_fp_13, 0, numParts);
        cuMemsetD32(part_fp_21, 0, numParts);
        cuMemsetD32(part_fp_23, 0, numParts);
        cuMemsetD32(part_fp_31, 0, numParts);
        cuMemsetD32(part_fp_32, 0, numParts);
    }
	
	//allocate dev mem for all objects based on number of grid cells
	private void initCUDAMemPtrs_Grids(int gridSzFltSz) {        
        cuMemAlloc(grid_mass, gridSzFltSz); 
        cuMemAlloc(grid_pos_x, gridSzFltSz);     
        cuMemAlloc(grid_pos_y, gridSzFltSz);   
        cuMemAlloc(grid_pos_z, gridSzFltSz); 
        cuMemAlloc(grid_vel_x, gridSzFltSz);  
        cuMemAlloc(grid_vel_y, gridSzFltSz);     
        cuMemAlloc(grid_vel_z, gridSzFltSz);        
 
        cuMemAlloc(grid_newvel_x, gridSzFltSz);
        cuMemAlloc(grid_newvel_y, gridSzFltSz);    
        cuMemAlloc(grid_newvel_z, gridSzFltSz);        
 
        cuMemAlloc(grid_force_x, gridSzFltSz);   
        cuMemAlloc(grid_force_y, gridSzFltSz);
        cuMemAlloc(grid_force_z, gridSzFltSz);		
        
		cuMemsetD32(grid_mass, 0, gridSize);
		cuMemsetD32(grid_vel_x, 0, gridSize);
		cuMemsetD32(grid_vel_y, 0, gridSize);
		cuMemsetD32(grid_vel_z, 0, gridSize);
		cuMemsetD32(grid_newvel_x, 0, gridSize);
		cuMemsetD32(grid_newvel_y, 0, gridSize);
		cuMemsetD32(grid_newvel_z, 0, gridSize);
		cuMemsetD32(grid_force_x, 0, gridSize);
		cuMemsetD32(grid_force_y, 0, gridSize);
		cuMemsetD32(grid_force_z, 0, gridSize);	
	}
	
	
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
	protected float createSphere(TreeMap<String, ArrayList<Float[]>> partVals, int numParts, float[] initVel, float[] ctr) {
    
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

        int gridSzFltSz = gridSize * Sizeof.FLOAT;
        initCUDAMemPtrs_Grids(gridSzFltSz);
        
        numBlocksParticles=numParts/numCUDAThreads+1;
        numBlocksGrid=gridSize/numCUDAThreads+1;
        kernelParams = new TreeMap<String, Pointer>();
        funcGridDimAndMemSize = new HashMap<String, int[]>();
        
        kernelParams.put("projectToGridandComputeForces",Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
        		Pointer.to(new float[] {this.mat.getLambda0()}), Pointer.to(new float[] {this.mat.getMu0()}), Pointer.to(new float[] {this.mat.hardeningCoeff}),
        		Pointer.to(part_mass), Pointer.to(part_vol),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(part_vel_x),Pointer.to(part_vel_y),Pointer.to(part_vel_z),

				Pointer.to(part_fe_11), Pointer.to(part_fe_12), Pointer.to(part_fe_13),
				Pointer.to(part_fe_21), Pointer.to(part_fe_22), Pointer.to(part_fe_23),
				Pointer.to(part_fe_31), Pointer.to(part_fe_32), Pointer.to(part_fe_33),

				Pointer.to(part_fp_11), Pointer.to(part_fp_12), Pointer.to(part_fp_13),
				Pointer.to(part_fp_21), Pointer.to(part_fp_22), Pointer.to(part_fp_23),
				Pointer.to(part_fp_31), Pointer.to(part_fp_32), Pointer.to(part_fp_33),

				Pointer.to(grid_mass),
				Pointer.to(grid_vel_x),Pointer.to(grid_vel_y),Pointer.to(grid_vel_z),
				Pointer.to(grid_force_x),Pointer.to(grid_force_y),Pointer.to(grid_force_z)));
        funcGridDimAndMemSize.put("projectToGridandComputeForces", new int[] {numBlocksParticles*4, 0});
		
        kernelParams.put("projectToGridInit", Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
				Pointer.to(part_mass),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(grid_mass)));
        funcGridDimAndMemSize.put("projectToGridInit", new int[] {numBlocksParticles, 0});
		
		kernelParams.put("computeVol", Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
				Pointer.to(part_mass),Pointer.to(part_vol),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(grid_mass)));
		funcGridDimAndMemSize.put("computeVol", new int[] {numBlocksParticles, 0});
		
		kernelParams.put("updPartVelocities",Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
				Pointer.to(new float[] {this.mat.alphaPicFlip}),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(part_vel_x),Pointer.to(part_vel_y),Pointer.to(part_vel_z),
				Pointer.to(grid_vel_x),Pointer.to(grid_vel_y),Pointer.to(grid_vel_z),				
				Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z)));  
		funcGridDimAndMemSize.put("updPartVelocities", new int[] {numBlocksParticles*4, numCUDAThreads*6*Sizeof.FLOAT});
		
	    kernelParams.put("compGridVelocities", Pointer.to(
				Pointer.to(new int[] {gridSize}), Pointer.to(new float[] {gravity.data[0]}),Pointer.to(new float[] {gravity.data[1]}),Pointer.to(new float[] {gravity.data[2]}),Pointer.to(new float[] {delT}),
				Pointer.to(grid_mass),
		
				Pointer.to(grid_vel_x),Pointer.to(grid_vel_y),Pointer.to(grid_vel_z),				
				Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z),
				Pointer.to(grid_force_x),Pointer.to(grid_force_y),Pointer.to(grid_force_z)));
		funcGridDimAndMemSize.put("compGridVelocities", new int[] {numBlocksGrid, 0});
	    
	    kernelParams.put("clearGrid", Pointer.to(
	    		Pointer.to(new int[] {gridSize}), Pointer.to(grid_mass),		
				Pointer.to(grid_vel_x),Pointer.to(grid_vel_y),Pointer.to(grid_vel_z),				
				Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z),
				Pointer.to(grid_force_x),Pointer.to(grid_force_y),Pointer.to(grid_force_z)));
		funcGridDimAndMemSize.put("clearGrid", new int[] {numBlocksGrid, 0});

	    kernelParams.put("gridCollisions", Pointer.to(
        		Pointer.to(new int[] {gridSize}),Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {h}), 
        		Pointer.to(new float[] {this.minSimBnds}),Pointer.to(new float[] {this.maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {delT}),Pointer.to(grid_mass),
        		Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z)));
		funcGridDimAndMemSize.put("gridCollisions", new int[] {numBlocksGrid, 0});

 
        kernelParams.put("partCollAndUpdPos", Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new float[] {this.minSimBnds}),Pointer.to(new float[] {this.maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {delT}),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(part_vel_x),Pointer.to(part_vel_y),Pointer.to(part_vel_z)));   
        funcGridDimAndMemSize.put("partCollAndUpdPos", new int[] {numBlocksParticles, 0});
        
        
        kernelParams.put("updDeformationGradient", Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {delT}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
        		Pointer.to(new float[] {this.mat.criticalCompression}), Pointer.to(new float[] {this.mat.criticalStretch}),
    			Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
    			
				Pointer.to(part_fe_11), Pointer.to(part_fe_12), Pointer.to(part_fe_13),
				Pointer.to(part_fe_21), Pointer.to(part_fe_22), Pointer.to(part_fe_23),
				Pointer.to(part_fe_31), Pointer.to(part_fe_32), Pointer.to(part_fe_33),

				Pointer.to(part_fp_11), Pointer.to(part_fp_12), Pointer.to(part_fp_13),
				Pointer.to(part_fp_21), Pointer.to(part_fp_22), Pointer.to(part_fp_23),
				Pointer.to(part_fp_31), Pointer.to(part_fp_32), Pointer.to(part_fp_33),
				
        		Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z)));     
        funcGridDimAndMemSize.put("updDeformationGradient", new int[] {numBlocksParticles, 0});     
        dispMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Finished CUDA Init | Launch first MPM Pass.", MsgCodes.info1);
	 	//launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {
        	String key = initStepFuncKeys[j];    
        	//System.out.println("Building Init Kernels Key : " + key);
        	int[] kernelVal = funcGridDimAndMemSize.get(key);
            cuLaunchKernel(cuFuncs.get(key), 
            		kernelVal[0], 1, 1,           // Grid dimension 
                    numCUDAThreads, 1, 1,  // Block dimension
                    kernelVal[1], null,           // Shared memory size and stream 
                    kernelParams.get(key), null);// Kernel- and extra parameters
        }
    	setSimFlags(gridIsBuiltIDX, true);
    	setSimFlags(simIsBuiltIDX, true);
    	dispMessage("MPM_Abs_CUDASim : "+simName, "cudaSetup","Finished first MPM Pass.", MsgCodes.info1);
	}//cudaSetup
    
	//compiles Ptx file from file in passed file name -> cuFileName needs to have format "xxxxx.cu"
	public void compilePtxFile(String krnFileName, String ptxFileName) throws IOException {
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

		String errorMessage = new String(toByteArray(process.getErrorStream())), outputMessage = new String(toByteArray(process.getInputStream()));
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
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        byte buffer[] = new byte[8192];
        while (true){
            int read = inputStream.read(buffer);
            if (read == -1){       break;     }
            baos.write(buffer, 0, read);
        }
        return baos.toByteArray();
    }

	
	//call whenever grid dimensions change
	public void setGridValsAndInit(int _gridCount, float _h, int _numParts) {
		//# of grid cells per side of cube
		gridCount = _gridCount;
		//# of particles
		numParts = _numParts;
		//cell size
		h = _h;	
		half_h = h*.5f;
		h3 = h*h*h;
		maxSimBnds = (gridCount*h)/2.0f;
		minSimBnds = -maxSimBnds;
		gridDim = maxSimBnds - minSimBnds;		
		//scale amount to fill 1500 x 1500 x 1500 visualization cube
		sclAmt = 1500.0f/(gridCount * h);

		//notify gridbuilder that we have a new grid
		//gridThdMgr.setNewSimVals();
		setSimFlags(gridIsBuiltIDX, false);
		//partThdMgr.setNewSimVals();
		setSimFlags(simIsBuiltIDX, false);
		//build grid - upon completion the gridbuilder will call resetSim to build particles
		resetSim(true);
	}//setGridValsAndInit	


	//call whenever setting/resetting simulation world - no changes to particle count/grid dimensions, just resetting to initial state
	//does not rebuild thread worker list, just re-executes existing thread workers
	public void resetSim(boolean freeCuda) {
		if(freeCuda) {
			//reset active ptrs
			JCuda.cudaFree(part_mass);		  
	        JCuda.cudaFree(part_vol); 
	        JCuda.cudaFree(part_pos_x); JCuda.cudaFree(part_pos_y);JCuda.cudaFree(part_pos_z);        
	
	        JCuda.cudaFree(part_vel_x); JCuda.cudaFree(part_vel_y);JCuda.cudaFree(part_vel_z);        
	
	        JCuda.cudaFree(part_fe_11);JCuda.cudaFree(part_fe_12);JCuda.cudaFree(part_fe_13);  
	        JCuda.cudaFree(part_fe_21);JCuda.cudaFree(part_fe_22);JCuda.cudaFree(part_fe_23); 
	        JCuda.cudaFree(part_fe_31);JCuda.cudaFree(part_fe_32);JCuda.cudaFree(part_fe_33);        
	    
	        JCuda.cudaFree(part_fp_11);JCuda.cudaFree(part_fp_12);JCuda.cudaFree(part_fp_13);
	        JCuda.cudaFree(part_fp_21);JCuda.cudaFree(part_fp_22);JCuda.cudaFree(part_fp_23);   
	        JCuda.cudaFree(part_fp_31);JCuda.cudaFree(part_fp_32);JCuda.cudaFree(part_fp_33);       
	 
	        JCuda.cudaFree(grid_mass); 
	        JCuda.cudaFree(grid_pos_x);JCuda.cudaFree(grid_pos_y);JCuda.cudaFree(grid_pos_z); 	  
	        JCuda.cudaFree(grid_vel_x);JCuda.cudaFree(grid_vel_y);JCuda.cudaFree(grid_vel_z); 	 
	        JCuda.cudaFree(grid_newvel_x);JCuda.cudaFree(grid_newvel_y);JCuda.cudaFree(grid_newvel_z); 
	        JCuda.cudaFree(grid_force_x);JCuda.cudaFree(grid_force_y);JCuda.cudaFree(grid_force_z);
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
	
	//called by partThdMgr runnable after particles are all built
	public void resetSimEnd() {
		System.out.println("Finished rebuilding particles of sim : " + getCurTime());
		//initial "sim" run - sets initial particle volume and density
		resetSimEnd_Priv();		
	}//resetSimEnd

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

	//execute simStepsPerFrame step of simulation
	//modAmtMillis is in milliseconds, counting # of ms since last sim call
	//returns true when simulation run is complete - true turns run sim flag off
	@SuppressWarnings("unused")
	public boolean simMe(float modAmtMillis) {
		/*cuDeviceGet(dev, 0);
	    cuCtxCreate(pctx, 0, dev);*/
	 	JCudaDriver.cuCtxSetCurrent(pctx);
	 	int updPartVelMemSz =  numCUDAThreads*6*Sizeof.FLOAT;
	 	int gridSzFltSz = gridSize * Sizeof.FLOAT;
 		for (int i=0;i<simStepsPerFrame;++i) { 	        
			for (int j=0;j<simStepFuncKeys.length;++j) {
				String key = simStepFuncKeys[j];
				//System.out.println("Running kernel : " + key);				
	        	int[] kernelVal = funcGridDimAndMemSize.get(key);
	            cuLaunchKernel(cuFuncs.get(key), 
	            		kernelVal[0], 1, 1,           // Grid dimension 
	                    numCUDAThreads, 1, 1,  // Block dimension
	                    kernelVal[1], null,           // Shared memory size and stream 
	                    kernelParams.get(key), null);// Kernel- and extra parameters
			}
		}
 		
		//copy from device data to host particle position arrays
		cuMemcpyDtoH(Pointer.to(h_part_pos_x),part_pos_x, numPartsFloatSz);
		cuMemcpyDtoH(Pointer.to(h_part_pos_y),part_pos_y, numPartsFloatSz);
		cuMemcpyDtoH(Pointer.to(h_part_pos_z),part_pos_z, numPartsFloatSz);
		return false;
	}//simMe}
	
	//sim method to show execution time for each step
	public abstract boolean simMeDebug(float modAmtMillis);	//simMeDebug	

	public static float getDeltaT() {return MPM_Abs_CUDASim.deltaT;}
	public void setDeltaT(float _delT) {MPM_Abs_CUDASim.deltaT = _delT;}

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
		
	//draw 1 frame of results	//animTimeMod is in seconds, counting # of seconds since last draw
	public void drawMe(float animTimeMod) {
		if(!getSimFlags(simIsBuiltIDX)) {return;}//if not built yet, don't try to draw anything
		//render all particles - TODO determine better rendering method
		pa.pushMatrix();pa.pushStyle();
		pa.strokeWeight(3.0f/sclAmt);
		pa.scale(sclAmt);	

		pa.pushMatrix();pa.pushStyle();
		//draw the points
		int pincr = 1;
		pa.beginShape(PConstants.POINTS);
		for(int i=0;i<=numParts-pincr;i+=pincr) {				
			//pa.stroke(h_part_clr[i][0], h_part_clr[i][1], h_part_clr[i][2]);
			pa.stroke(h_part_clr_int[i][0], h_part_clr_int[i][1], h_part_clr_int[i][2]);
			//pa.point(h_part_pos_x[i], h_part_pos_y[i], h_part_pos_z[i]);
			pa.vertex(h_part_pos_x[i], h_part_pos_y[i], h_part_pos_z[i]);
		}
		pa.endShape();
		pa.popStyle();pa.popMatrix();
		
		if(getSimFlags(showGrid)) {
			pa.pushMatrix();pa.pushStyle();			
			pa.stroke(0,0,0,20);
			pa.translate(minSimBnds,minSimBnds,minSimBnds);
			int incr = 10;
			//shows every "incr" gridcells
			for (int i=0; i<=gridCount;i+=incr) {
				float iLoc = i*h;
				for(int j=0;j<=gridCount;j+=incr) {
					myVectorf startPos=new myVectorf(iLoc,j*h,0.0f);
					myVectorf endPos=new myVectorf(iLoc,j*h,gridDim);
					pa.line(startPos,endPos);
				}
				for(int k=0;k<=gridCount;k+=incr) {
					myVectorf startPos=new myVectorf(iLoc,0.0f, k*h);
					myVectorf endPos=new myVectorf(iLoc,gridDim,k*h);
					pa.line(startPos,endPos);
				}
			}
			for(int j=0;j<=gridCount;j+=incr) {
				float jLoc = j*h;
				for(int k=0;k<=gridCount;k+=incr) {
					myVectorf startPos=new myVectorf(0.0f,jLoc,k*h);
					myVectorf endPos=new myVectorf(gridDim,jLoc,k*h);
					pa.line(startPos,endPos);
				}
			}
			pa.popStyle();pa.popMatrix();
		}
		pa.popStyle();pa.popMatrix();
	}//drawMe
	
	
	//reinitialize instanced Simulation- set up all resources - re call this every time new sim is being set up
	protected abstract void resetSimEnd_Priv();
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

//instance of sim world with 2 big snow boulders slamming into each other 
class MPM_Cuda2Balls extends MPM_Abs_CUDASim {
	//scale w/timestep
	private static float initVel = 30.0f;
	
	public MPM_Cuda2Balls(my_procApplet _pa,int _gridCount, float _h, int _numParts) {
		super(_pa,"2 Big Snowballs",_gridCount, _h,_numParts);
	}
	
	@Override
	protected void resetSimEnd_Priv() {

	}//resetSimPriv
	
	//build particle layout for cuda sim - use multiples of h as radius
	@Override
	protected TreeMap<String, ArrayList<Float[]>> buildPartLayout() {	
        //create particle layout
        TreeMap<String, ArrayList<Float[]>> partVals = new TreeMap<String, ArrayList<Float[]>>();
        partVals.put("pos",new ArrayList<Float[]>());
        partVals.put("vel",new ArrayList<Float[]>());
        partVals.put("minMaxVals",new ArrayList<Float[]>());
        //initialize min and max values
        partVals.get("minMaxVals").add(new Float[] {100000.0f,100000.0f,100000.0f});
        partVals.get("minMaxVals").add(new Float[] {-100000.0f,-100000.0f,-100000.0f});      

        int numPartsPerSphere = numParts/2;
        //float ctrOfGrid = (minSimBnds + maxSimBnds)/2.0f, diff = maxSimBnds - minSimBnds; 
        //float qtrDiff = .25f*diff, halfDiff = 2*qtrDiff;
       
		//Float[] minVals = partVals.get("minMaxVals").get(0);
		//Float[] maxVals = partVals.get("minMaxVals").get(1);       
		//float hSq = h*h;
		//scale amt is width divided by # of grid cells and size of each cell
		//float scaledRad = 150.0f/this.sclAmt;		//was .25
		//float sphereRad = scaledRad/h;
		//float sphereSqRad = sphereRad*sphereRad;
		float offScl = 600.0f/this.sclAmt;

        float xVel = -1.4f*initVel, 
        	yVel = 0, 
    		zVel = .5f*-initVel;
		
		float xOff = .4f  *offScl, 
			yOff = .5f*offScl, 
			zOff = .4f *offScl;
		//int incr = 1;
		//lower ball
        //buildSphere(partVals,new float[] {xOff, yOff, zOff}, new float [] {xVel, yVel, zVel}, incr, sphereSqRad, minVals, maxVals);
		float sphereRad = createSphere(partVals,numPartsPerSphere, new float [] {xVel, yVel, zVel}, new float[] {xOff, yOff, zOff});
				
		//find min and max values for 1st built sphere
//		for(int i=0;i<3;++i) {
//			float rad = (maxVals[i] - minVals[i]) * .5f;
//			System.out.println("First built sphere idx "+i+" max : " + maxVals[i]+" and min : " + minVals[i] + " -> unscaled radius : " + rad + " Scaled (apparent) radius : "+ (rad*this.sclAmt));
//		}
		xVel *= -1;
        yVel *= -1;
        xOff = -.1f*sphereRad *offScl;
        yOff = .5f*offScl;
        zOff = (.4f*sphereRad)*offScl;
        //upper ball        
        //buildSphere(partVals,new float[] {xOff, yOff, zOff}, new float [] {xVel, yVel, zVel}, incr, sphereSqRad, minVals, maxVals);
		createSphere(partVals,numPartsPerSphere, new float [] {xVel, yVel, zVel}, new float[] {xOff, yOff, zOff});

         //end create particle layout	
		return partVals;
	}//buildPartLayout

	
	@Override
	//draw scene-specific collider, if it existss
	protected void drawCollider(float animTimeMod) {}
	
	//sim method to show execution time for each step
	public boolean simMeDebug(float modAmtMillis) {
		//TODO any debugging that might be supportable here
		return false;
	}//simMeDebug
}//class MPM_Cuda2Balls




