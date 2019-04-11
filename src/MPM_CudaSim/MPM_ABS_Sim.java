package MPM_CudaSim;

import static jcuda.driver.JCudaDriver.*;

import java.awt.Color;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadLocalRandom;

import org.jblas.FloatMatrix;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2ES2;
import com.jogamp.opengl.GL3ES3;

import jcuda.*;
import jcuda.driver.*;
import jcuda.runtime.JCuda;
import processing.opengl.PGL;
import processing.opengl.PJOGL;

//abstract class describing a simulation world.  called by sim window to executed simulation and to render results
//instancing classes can hold different configurations/initialization setups
public abstract class MPM_ABS_Sim{	
	//runnable to launch threads to manage particle minipulations
//	public myPartBuilder partThdMgr;
//	//runnable to build and manage grid manipulations
//	public myGridBuilder gridThdMgr;
	
	//cuda kernel file name
	private String ptxFileName = "MPM_ABS_Sim_New.ptx";	
	//vs c++ compiler location
	//private final static String VSCppCompLoc = "C:\\Program Files (x86)\\Microsoft Visual Studio 14.0\\VC\\bin\\x86_amd64\\";

	//grid dim - cube so same in all 3 dims; # of particles in sim
	protected int gridCount;

	public int numParts = 130234;
	public static int numPartsUI_Init = 130234;
	//# parts * size of float - currently always builds same # of particles
	protected long numPartsFloatSz;
	
	//simulation boundaries - symmetric cube, only need min and max, grid length per dim
	protected float minSimBnds, maxSimBnds, gridDim;
	//dimension of grid cells
	protected float h, half_h, h3;
	
	//particle sphere to demonstrate simulation
	myVectorf snoBallCtr;
	float snoBallRad;
	
	//parameters
	//timestep of simulation - 
	private static float deltaT = 3.5e-4f;
	//scale w/timestep
	private static float initVel = 30.0f;
	//snow density varies from 50 to ~800
	//powdery snow is about 100
	//wet firm compacted snow is about 600
	private static float initDensity = 100.0f;
	//gravity 
	public final float gGravity = -9.8f;
	//const matrix for calculations - z is down
	final FloatMatrix gravity = new FloatMatrix(new float[] {0, 0, gGravity});
	
	//particle mass
	//public static float pMass = 100.0f*(.01f*.01f*.01f)/4;
	
	//scale amount for visualization to fill cube frame in 3d world; particle radius, scaled by different visual scales
	protected float sclAmt;
	protected float partRad;
	
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
	
	//for multithreading 
	public ExecutorService th_exec;
	public int numThreadsAvail;	
	
	//iterations per frame
	public static int simStepsPerFrame = 1;
	//time of simulation start, in millis
	private long simStartTime;

	//time exec built, in millis - used as offset for instant to provide smaller values for timestamp
	private final long execBuiltTime;
	//constants for collider calcs
	protected static float cyl_da = (float) (Math.PI/18.0f), TWO_PI = (float) (2.0 * Math.PI);	
	
	protected TreeMap<String, Pointer> kernelParams;
	protected TreeMap<String, CUfunction> cuFuncs;
	protected String[] CUFileFuncNames = new String[] {"projectToGridandComputeForces","projectToGridInit", "computeVol", "updPartVelocities",  "compGridVelocities", "partCollAndUpdPos", "gridCollisions", "updDeformationGradient" };
	protected String[] initStepFuncKeys = new String[] {"projectToGridInit", "computeVol","compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "partCollAndUpdPos"}; 
	protected String[] simStepFuncKeys = new String[] {"projectToGridandComputeForces", "compGridVelocities", "gridCollisions", "updDeformationGradient", "updPartVelocities", "partCollAndUpdPos"}; 
   
    
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
    float[][] h_part_clr;
    
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
	@SuppressWarnings("unchecked")
	public MPM_ABS_Sim(int _gridCount, float _h) {		
		//for multithreading - do not use instanced version in PApplet - we may not use processing-based build to run simulation
		th_exec = Executors.newCachedThreadPool();		
		int numThreadsTtl = Runtime.getRuntime().availableProcessors();
		//# of threads to use to build structures - doesn't change
		numThreadsAvail = (numThreadsTtl > 2 ? numThreadsTtl-2 : 1);

		Instant now = Instant.now();
		execBuiltTime = now.toEpochMilli();//milliseconds since 1/1/1970 when this sim was built.
		//mat's quantities are managed by UI - only need to instance once
		mat = new myMaterial();
		//initialize active nodes set - array of sets, array membership is node ID % numThreadsAvail
		//might not be balanced, but worst case all nodes in single thread.  not many nodes, and minimal calculation is node-specific ConcurrentHashMap<myGridNode,Boolean>[]
//		activeNodes = new ConcurrentHashMap[numThreadsAvail];
//		for(int i=0;i<this.numThreadsAvail;++i) {activeNodes[i] = new ConcurrentHashMap<myGridNode, activeNodeAgg>();}	
		//set up threads to be used to build the grid
//		gridThdMgr = new myGridBuilder(this, numThreadsAvail);
//		//set up thread builders to be used to build particles
//		partThdMgr = new myPartBuilder(this, numThreadsAvail);		
		//setup flag array
		
		initSimFlags();
		
		//set up grid and initialize sim
		setGridValsAndInit(_gridCount, _h);
		try {
			cudaSetup();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}//MPM_ABS_Sim
	//run 1 time to load kernel and assign function pointers to functions
	private void initOnceCUDASetup() {
		// Enable exceptions and omit all subsequent error checks
        JCudaDriver.setExceptionsEnabled(true); 
        //Initialize the driver and create a context for the first device. (device 0
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
			CUfunction c = new CUfunction();
			System.out.println("Key : " + key);
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
        h_part_clr = new float[numParts][3];
        
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
        	h_part_clr[i] = new float[] {
    			getClrVal(h_part_pos_x[i],minVals[0],maxVals[0]), 
				getClrVal(h_part_pos_y[i],minVals[1],maxVals[1]), 
				getClrVal(h_part_pos_z[i],minVals[2],maxVals[2]) 
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
	}	
	
	//create sphere with given center, with passed # of particles
	private void createSphere(TreeMap<String, ArrayList<Float[]>> partVals, int numParts, float yVel, float[] ctr) {
		//radius value for particles for display - call after numParts specified
        float rawPartRad = 2.0f, partRad = calcScaledVal(rawPartRad);        
		//build sphere of particles - scale volume of sphere based on cuberoot of # of particles, with 1000 particles being baseline sphere 
        float ballRad = (float) (3.0*Math.cbrt(numParts) * partRad/2.0f);		
		System.out.println("part rad : "+ partRad+" ball rad : "+ ballRad);

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
			partVals.get("vel").add(new Float[] {0f,yVel,0f});
        }
	}//createSphere

	//build particle layout for cuda sim - use multiples of h as radius
	private TreeMap<String, ArrayList<Float[]>> buildPartLayout() {		
        //create particle layout
        TreeMap<String, ArrayList<Float[]>> partVals = new TreeMap<String, ArrayList<Float[]>>();
        partVals.put("pos",new ArrayList<Float[]>());
        partVals.put("vel",new ArrayList<Float[]>());
        partVals.put("minMaxVals",new ArrayList<Float[]>());
        //initialize min and max values
        partVals.get("minMaxVals").add(new Float[] {100000.0f,100000.0f,100000.0f});
        partVals.get("minMaxVals").add(new Float[] {-100000.0f,-100000.0f,-100000.0f});      


        float xVel = -initVel, 
        	yVel = initVel, 
        		zVel = 0;//-initVel/3.0f;
        int numPartsPerSphere = 130234/2;
        float ctrOfGrid = (minSimBnds + maxSimBnds)/2.0f, diff = maxSimBnds - minSimBnds; 
        float qtrDiff = .25f*diff, halfDiff = 2*qtrDiff;
        float[] ctr1 = new float[] {(ctrOfGrid+qtrDiff), (ctrOfGrid - halfDiff), ctrOfGrid};
        float[] ctr2 = new float[] {(ctrOfGrid+qtrDiff), (ctrOfGrid + halfDiff), (ctrOfGrid + qtrDiff)};
        //create sphere at certain location based on grid dimensions (so spheres are same size regardless of dimensions)
        //createSphere(partVals,numPartsPerSphere, yVel, );
       
		Float[] minVals = partVals.get("minMaxVals").get(0);
		Float[] maxVals = partVals.get("minMaxVals").get(1);       
		//float hSq = h*h;
		//scale amt is width divided by # of grid cells and size of each cell
		float scaledRad = 150.0f/this.sclAmt;		//was .25
		float sphereRad = scaledRad/h;
		float sphereSqRad = sphereRad*sphereRad;
		float offScl = 600.0f/this.sclAmt;
		float offsetBase = .25f;
		float xOff = offsetBase *offScl, 
				yOff = (-.25f + offsetBase)*offScl, 
				zOff = offsetBase*offScl;
		int incr = 1;
		//lower ball
        buildSphere(partVals,new float[] {xOff, yOff, zOff}, new float [] {xVel, yVel, zVel}, incr, sphereSqRad, minVals, maxVals);
    
        float rawPartRad = 2.0f;        
        float partRad = calcScaledVal(rawPartRad), ballRad = (float) (3.0*Math.cbrt(numParts) * partRad/2.0f);		
		System.out.println("Sphere Sq Rad : " + sphereSqRad+ "\t|Scale Amount : " + this.sclAmt+"|part rad : "+ partRad+" ball rad : "+ ballRad+ " * scl amt :  " + (this.sclAmt*ballRad));
		//find min and max values for 1st built sphere
		for(int i=0;i<3;++i) {
			float rad = (maxVals[i] - minVals[i]) * .5f;
			System.out.println("First built sphere idx "+i+" max : " + maxVals[i]+" and min : " + minVals[i] + " -> unscaled radius : " + rad + " Scaled (apparent) radius : "+ (rad*this.sclAmt));
		}
		xVel *= -1;
        yVel *= -1;
        xOff = (-.15f + offsetBase)*offScl;
        yOff = (.15f + offsetBase)*offScl;
        zOff = (.15f + offsetBase)*offScl;
        //upper ball        
        buildSphere(partVals,new float[] {xOff, yOff, zOff}, new float [] {xVel, yVel, zVel}, incr, sphereSqRad, minVals, maxVals);

         //end create particle layout	
		return partVals;
	}//buildPartLayout
	//this will build a single sphere of particles
	private void buildSphere(TreeMap<String, ArrayList<Float[]>> partVals,float[] offset, float [] vel, int incr, float sphereSqRad, Float[] minVals, Float[] maxVals) {
		
        for(int i = -25; i < 25; i+=incr){
        	float xpos=offset[0]+h*i, iSqM25 = i*i;
        	for(int j = -25; j < 25; j+=incr){
        		float ypos=offset[1]+h*j, ijSqM25 = iSqM25 + j*j;
        		for(int k = -25; k < 25; k+=incr){
        			float zpos=offset[2]+h*k;
        			float sqmagn=((ijSqM25+(k*k)));
        			if(sqmagn<sphereSqRad) {
        				Float[] vals = new Float[] {xpos+(float)Math.random()*h, ypos+(float)Math.random()*h,zpos+(float)Math.random()*h};
        				
        				for (int v = 0; v < 3; ++v) {
        					minVals[v] = (vals[v] < minVals[v] ?vals[v] : minVals[v] );
        					maxVals[v] = (vals[v] > maxVals[v] ?vals[v] : maxVals[v] );
        				}
        				partVals.get("pos").add(vals);
        				//init vel
        				partVals.get("vel").add(new Float[] {vel[0],vel[1],vel[2]});
        			}        			
                }        	
            }
        }     
	
	}//buildSphere
		
	public void cudaSetup() throws IOException {        
		if (!this.getSimFlags(CUDADevInit)) {
            //init cuda device and kernel file if not done already
            this.initOnceCUDASetup();
        }

		TreeMap<String, ArrayList<Float[]>> partVals = buildPartLayout();
		//total grid size       
        gridSize=gridCount*gridCount*gridCount;
        
        numParts = partVals.get("pos").size();
        System.out.println("# of particles : " + numParts);
        numPartsFloatSz = numParts * Sizeof.FLOAT;

        //init ptrs to particle-based arrays - numparts and numPartsFloatSz need to be initialized here
        this.initCUDAMemPtrs_Parts(partVals);

        //init grid ptrs
        float delT = this.getDeltaT();

        int gridSzFltSz = gridSize * Sizeof.FLOAT;
        initCUDAMemPtrs_Grids(gridSzFltSz);
        
        numBlocksParticles=numParts/numCUDAThreads+1;
        numBlocksGrid=gridSize/numCUDAThreads+1;
//	protected String[] CUFileFuncNames = new String[] {"projectToGridandComputeForces","projectToGridInit", "computeVol","getPos", "updPartVelocities", "updPartPositions", "compGridVelocities", "particleCollisions", "gridCollisions", "updDeformationGradient" };
        kernelParams = new TreeMap<String, Pointer>();
                
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
		
        kernelParams.put("projectToGridInit", Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
				Pointer.to(part_mass),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(grid_mass)));
		
		kernelParams.put("computeVol", Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
				Pointer.to(part_mass),Pointer.to(part_vol),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(grid_mass)));
		
		kernelParams.put("updPartVelocities",Pointer.to(
				Pointer.to(new int[] {numParts}), Pointer.to(new int[] {gridCount}), Pointer.to(new float[] {h}), Pointer.to(new float[] {this.minSimBnds}),
				Pointer.to(new float[] {this.mat.alphaPicFlip}),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(part_vel_x),Pointer.to(part_vel_y),Pointer.to(part_vel_z),
				Pointer.to(grid_vel_x),Pointer.to(grid_vel_y),Pointer.to(grid_vel_z),				
				Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z)));  
		
	    kernelParams.put("compGridVelocities", Pointer.to(
				Pointer.to(new int[] {gridSize}), Pointer.to(new float[] {gravity.data[0]}),Pointer.to(new float[] {gravity.data[1]}),Pointer.to(new float[] {gravity.data[2]}),Pointer.to(new float[] {delT}),
				Pointer.to(grid_mass),
		
				Pointer.to(grid_vel_x),Pointer.to(grid_vel_y),Pointer.to(grid_vel_z),				
				Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z),
				Pointer.to(grid_force_x),Pointer.to(grid_force_y),Pointer.to(grid_force_z)));

	    kernelParams.put("gridCollisions", Pointer.to(
        		Pointer.to(new int[] {gridSize}),Pointer.to(new int[] {gridCount}),Pointer.to(new float[] {h}), 
        		Pointer.to(new float[] {this.minSimBnds}),Pointer.to(new float[] {this.maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {delT}),Pointer.to(grid_mass),
        		Pointer.to(grid_newvel_x),Pointer.to(grid_newvel_y),Pointer.to(grid_newvel_z)));

 
        kernelParams.put("partCollAndUpdPos", Pointer.to(
        		Pointer.to(new int[] {numParts}), Pointer.to(new float[] {this.minSimBnds}),Pointer.to(new float[] {this.maxSimBnds}),
        		Pointer.to(new float[] {wallFric}),Pointer.to(new float[] {delT}),
				Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z),
				Pointer.to(part_vel_x),Pointer.to(part_vel_y),Pointer.to(part_vel_z)));   
        
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
        CUdeviceptr graphicsPointer=null;
//        kernelParams.put("getPos", Pointer.to(
//        		Pointer.to(new int[] {numParts}), Pointer.to(graphicsPointer),        		
//        		Pointer.to(part_pos_x),Pointer.to(part_pos_y),Pointer.to(part_pos_z)));        
 
	 	int[] initGridDims = new int[] {numBlocksParticles, numBlocksParticles,numBlocksGrid,numBlocksGrid,numBlocksParticles, numBlocksParticles*4,	numBlocksParticles,numBlocksParticles};
	 	int[] initSharedMemSz = new int[] {0,0,0,0,0,numCUDAThreads*6*Sizeof.FLOAT,0,0};
        
	 	//launch init functions
        for (int j=0;j<initStepFuncKeys.length;++j) {
        	String key = initStepFuncKeys[j];        	
            cuLaunchKernel(cuFuncs.get(key), 
            		initGridDims[j], 1, 1,           // Grid dimension 
                    numCUDAThreads, 1, 1,  // Block dimension
                    initSharedMemSz[j], null,           // Shared memory size and stream 
                    kernelParams.get(key), null // Kernel- and extra parameters
                );
        }
	 	
	   
        
    	setSimFlags(gridIsBuiltIDX, true);
    	setSimFlags(simIsBuiltIDX, true);
	}
    
	//compiles Ptx file from file in passed file name -> cuFileName needs to have format "xxxxx.cu"
	public void compilePtxFile(String krnFileName, String ptxFileName) throws IOException {
	//NOTE : using new version of CUDA (as of 8/7/18) w/vs2015 compiles this incorrectly/makes it hang. TODO need to investigate this
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
    
	//returns scaled value of passed value, given sclAmt has been calced - used for display purposes, scales by amount grid is increased to fill display cube
	protected float calcScaledVal(float val) {	return val/sclAmt;}
	
	//call whenever grid dimensions change
	public void setGridValsAndInit(int _gridCount, float _h) {
		//# of grid cells per side of cube
		gridCount = _gridCount;
		//cell size
		h = _h;	
		half_h = h*.5f;
		h3 = h*h*h;
		maxSimBnds = (gridCount*h)/2.0f;
		minSimBnds = -maxSimBnds;
		gridDim = maxSimBnds - minSimBnds;		
		//scale amount to fill 1500 x 1500 x 1500 visualization cube
		sclAmt = 1500.0f/(gridCount * h);
			
		//initialize list of particles
		//parts = new myRndrdPart[numParts];
		//initialize snowball - dependent on # of particles and grid dimensions
		//initSnowBall();

		//notify gridbuilder that we have a new grid
		//gridThdMgr.setNewSimVals();
		setSimFlags(gridIsBuiltIDX, false);
		//partThdMgr.setNewSimVals();
		setSimFlags(simIsBuiltIDX, false);
		//build grid - upon completion the gridbuilder will call resetSim to build particles
		//th_exec.execute(gridThdMgr);
	}//setGridValsAndInit	
	
	
	private void initSnowBall() {
		//radius value for particles for display - call after numParts specified
		partRad = calcScaledVal(10.0f);
		//build sphere of particles - scale volume of sphere based on cuberoot of # of particles, with 1000 particles being baseline sphere 
		snoBallRad = (float) (5.0f * Math.cbrt(numParts) * partRad/10.0);		
		float ctrOfGrid = (minSimBnds + maxSimBnds)/2.0f;
		snoBallCtr = new myVectorf(ctrOfGrid,ctrOfGrid,.5f * (maxSimBnds - minSimBnds) +  minSimBnds );		
	}//initSnowBall
		
	
	//call with same grid but with new particle count
	public void setNumPartsAndReset(int _numParts) {
		//rebuild sim with new # of particles
		resetSim();
	}//setNumPartsAndReset

	//call whenever setting/resetting simulation world - no changes to particle count/grid dimensions, just resetting to initial state
	//does not rebuild thread worker list, just re-executes existing thread workers
	public void resetSim() {
		//reset active ptrs
		JCuda.cudaFree(part_mass);		  
        JCuda.cudaFree(part_vol); 
        JCuda.cudaFree(part_pos_x);
        JCuda.cudaFree(part_pos_y);
        JCuda.cudaFree(part_pos_z);        

        JCuda.cudaFree(part_vel_x);
        JCuda.cudaFree(part_vel_y); 
        JCuda.cudaFree(part_vel_z);        

        JCuda.cudaFree(part_fe_11);
        JCuda.cudaFree(part_fe_12);
        JCuda.cudaFree(part_fe_13);  
        JCuda.cudaFree(part_fe_21);  
        JCuda.cudaFree(part_fe_22);
        JCuda.cudaFree(part_fe_23); 
        JCuda.cudaFree(part_fe_31);   
        JCuda.cudaFree(part_fe_32);
        JCuda.cudaFree(part_fe_33);        
    
        JCuda.cudaFree(part_fp_11); 
        JCuda.cudaFree(part_fp_12);
        JCuda.cudaFree(part_fp_13);   
        JCuda.cudaFree(part_fp_21);  
        JCuda.cudaFree(part_fp_22);   
        JCuda.cudaFree(part_fp_23);   
        JCuda.cudaFree(part_fp_31);  
        JCuda.cudaFree(part_fp_32);    
        JCuda.cudaFree(part_fp_33);       
 
        JCuda.cudaFree(grid_mass); 
        JCuda.cudaFree(grid_pos_x);     
        JCuda.cudaFree(grid_pos_y);   
        JCuda.cudaFree(grid_pos_z);        
  
        JCuda.cudaFree(grid_vel_x);
        JCuda.cudaFree(grid_vel_y);     
        JCuda.cudaFree(grid_vel_z);        
 
        JCuda.cudaFree(grid_newvel_x);
        JCuda.cudaFree(grid_newvel_y);    
        JCuda.cudaFree(grid_newvel_z);        
 
        JCuda.cudaFree(grid_force_x);   
        JCuda.cudaFree(grid_force_y);
        JCuda.cudaFree(grid_force_z);
       
		//clearActiveNodes();
		//sim start time - time from when sim object was first instanced
		simStartTime = getCurTime();	
		System.out.println("Start rebuilding particles of sim : " + simStartTime);
		//partThdMgr.setSimStep(0);
		setSimFlags(CUDADevInit,false);
		setSimFlags(simIsBuiltIDX, false);
		try {
			this.cudaSetup();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		setSimFlags(simIsBuiltIDX, true);
	}//initOnce
	
	//called by partThdMgr runnable after particles are all built
	public void resetSimEnd() {
		System.out.println("Finished rebuilding particles of sim : " + getCurTime());
		//initial "sim" run - sets initial particle volume and density
		resetSimPriv();		
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

	//execute 1 step of simulation
	//modAmtMillis is in milliseconds, counting # of ms since last sim
	//returns true when simulation run is complete - true turns run sim flag off	
	public abstract boolean simMe(float modAmtMillis);	
	//sim method to show execution time for each step
	public abstract boolean simMeDebug(float modAmtMillis);	//simMeDebug	
	//sim method to run multi-threaded implementation
	public abstract boolean simMeMThd(float modAmtMillis);	//simMeDebug
	
	public abstract boolean simMeCuda(float modAmtMillis);	//simMeCuda
	
	public void gridCollisions() {
	}

	public static float getDeltaT() {return MPM_ABS_Sim.deltaT;}
	public void setDeltaT(float _delT) {MPM_ABS_Sim.deltaT = _delT;}
	
	//cube wall normals
	protected final myVectorf[] wallNorms = new myVectorf[] {
			new myVectorf( 1.0, 0.0, 0.0),new myVectorf(-1.0, 0.0, 0.0),
			new myVectorf( 0.0, 1.0, 0.0),new myVectorf( 0.0,-1.0, 0.0),
			new myVectorf( 0.0, 0.0, 1.0),new myVectorf( 0.0, 0.0,-1.0)
	};
	
	
	//collision detection for wall collisions, returns idx in wallNorms of collision
	public int checkWallCollision(myVectorf pos) {
		if(pos.x<=minSimBnds) {		return 0;	} else if(pos.x>=maxSimBnds) {	return 1;	}			
		if(pos.y<=minSimBnds) {		return 2;	} else if(pos.y>=maxSimBnds) {	return 3;	}
		if(pos.z<=minSimBnds) {		return 4;	} else if(pos.z>=maxSimBnds) {	return 5;	}
		return -1;
	}//checkWallCollision
	
	//check sim-specific central floating collider - return 0 vec if no collision, otherwise return normal of collision
	public abstract myVectorf checkColliderCollision(myVectorf pos);
	
	//looking ahead at position, to return new velocity that responds to collision if position has collided
	public myVectorf applyCollisions(myVectorf pos, myVectorf velocity) {
		//check collider collisions
		myVectorf sphereColNorm = checkColliderCollision(pos);
		if(sphereColNorm.magn != 0) {//colliding with central collider			
			return calcCollVel(collFric, velocity, sphereColNorm);
			
		} else {//not colliding with sphere collider, check walls	
			//check wall collisions
			int colType = checkWallCollision(pos);
			if(-1==colType) {return velocity;}
			//calc collision velocity if collision is going to occur
			return calcCollVel(wallFric, velocity, wallNorms[colType]);
		}
	}//applyCollisions
	
	//calc collision velocity if collision is going to occur - norm is collider normal
	public myVectorf calcCollVel(float fricCoeff, myVectorf velocity, myVectorf norm) {		
		float velNormDir=myVectorf._dot(velocity, norm);
		if(velNormDir >=0){return velocity;}
		//velocity in opposite direction of normal if dot prod is <0	
		myVectorf velTanDir = myVectorf._sub(velocity, myVectorf._mult(norm,velNormDir));
		float fricNormVel = fricCoeff*velNormDir;
		if(velTanDir.magn <= -fricNormVel) {//coulomb friction law
			//no bounce, so just stops
			return new myVectorf(0.0,0.0,0.0);
		}		
		return myVectorf._add(velTanDir,(myVectorf._normalize(velTanDir))._mult(fricNormVel));		 
	}//calcCollVel
	
	
	
//////////////////////////////////////////////////////////////////
//multi-threaded funcs
	//first clear operations
	protected void launchCalc_MT() {
		clearActiveNodes_MT();
	}
	
	protected void clearActiveNodes_MT() {
//        gridThdMgr.setSimStep(1);//set sim step in grid manager to "clear nodes" step
//        gridThdMgr.synchRun();
	}//clearActiveNodes_MT
	
	
//multi-threaded funcs end
//////////////////////////////////////////////////////////////////
	
	
	
	public void showTimeMsgSimStart(String _str) {System.out.println(_str+" Time Now : "+(getCurTime() - simStartTime));}
	//display message and time now
	public void showTimeMsgNow(String _str, long stTime) {		
		System.out.println(_str+" Time Now : "+(getCurTime() - stTime)+" ms");
	}
	
	//get time from "start time"
	//1518700615691 is epoch instant @ 8:17 am 2/15/18
	//public long getCurTime() {return getCurTime(1518700615691L);}
	public long getCurTime() {			
		Instant instant = Instant.now();
		long millis = instant.toEpochMilli() - execBuiltTime;//milliseconds since 1/1/1970, subtracting when this sim exec was built to keep millis low			
		return millis;
	}//getCurTime() 	
	//draw 1 frame of results
	//pa is the reference to the graphics engine (a processing applet).  use this to render results
	//animTimeMod is in seconds, counting # of seconds since last draw
	public abstract void drawMe(MPM_SimMain pa, float animTimeMod);	
	//reinitialize instanced Simulation- set up all resources - re call this every time new sim is being set up
	protected abstract void resetSimPriv();
	
	//some utility functions and constants
	protected myPointf Pf(myPointf O, float x, myVectorf I, double y, myVectorf J, double z, myVectorf K) {return new myPointf(O.x+x*I.x+y*J.x+z*K.x,O.y+x*I.y+y*J.y+z*K.y,O.z+x*I.z+y*J.z+z*K.z);}  // O+xI+yJ+zK
	//returns arrays of start and end points for cylinder of passed radius r going from point a to point b
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

//instance of sim world.  
class MPM_BaseSim extends MPM_ABS_Sim {
	//collider to demonstrate behavior
	//location of cylidrical collider (along z axis)
	protected myVectorf colLocation;
	//radius of cylinder
	protected float colRad, colSqRad;
	//points of cylinder, for display if using cylinder collider
	protected ArrayList<myPointf>[] cylPts;	
	
	public MPM_BaseSim(int _gridCount, float _h) {
		super(_gridCount, _h);
	}
		
	@Override
	protected void resetSimPriv() {
		//(re)build collider location and size
		//put collider location at 1/2 the distance between snowball and ground
		colLocation = new myVectorf(snoBallCtr.x, snoBallCtr.y, .25f*(snoBallCtr.z - minSimBnds)+ minSimBnds);
		//make colider radius 2x the radius of snowball
		colRad = snoBallRad * 1.5f; colSqRad = colRad * colRad;
		//make points to hold cylinder geometry, for display - along y axis
		//cylPts = buildCylinder(new myPointf(colLocation.x, maxSimBnds, colLocation.z), new myPointf(colLocation.x, minSimBnds, colLocation.z), colRad);
		
// 		//called once upon simulation start/reset
//		projectToGrid();
//		computePartVolumeAndDensity(); 
	}//resetSimPriv
	
	@Override
	//check central floating collider - return 0 vec if no collision, otherwise return normal of collision
	public synchronized myVectorf checkColliderCollision(myVectorf pos) {
		if(myPointf._SqrDist(pos, colLocation) <= colSqRad) {
			myVectorf tmp = new myVectorf(pos);
			tmp._sub(colLocation);
			return tmp._normalize();
		}
		return new myVectorf(0,0,0);
	}//checkColliderCollision
	
	
	//execute simStepsPerFrame step of simulation
	//modAmtMillis is in milliseconds, counting # of ms since last sim call
	//returns true when simulation run is complete - true turns run sim flag off
	public boolean simMeMThd(float modAmtMillis) {
		//launch either grid builder or part builder thread managers with appropriate sim step codes
//		for (int i=0;i<simStepsPerFrame;++i) {
//			//reset all active nodes from last iteration			
//			launchCalc_MT();
//
//		}
		return false;
	}//simMe}
	
	@Override
	//execute simStepsPerFrame step of simulation
	//modAmtMillis is in milliseconds, counting # of ms since last sim call
	//returns true when simulation run is complete - true turns run sim flag off
	public boolean simMe(float modAmtMillis) {

		return false;
	}//simMe}
	
	public boolean simMeCuda(float modAmtMillis) {
		/*cuDeviceGet(dev, 0);
	    cuCtxCreate(pctx, 0, dev);*/
	 	JCudaDriver.cuCtxSetCurrent(pctx);
	 	int updPartVelMemSz =  numCUDAThreads*6*Sizeof.FLOAT;
	 	
	 	int[] simGridSz = new int[] {numBlocksParticles*4,numBlocksGrid,numBlocksGrid,numBlocksParticles, numBlocksParticles*4,	numBlocksParticles,numBlocksParticles};
	 	int[] simSharedMemSz = new int[] {0,0,0,0,updPartVelMemSz,0,0};
 		for (int i=0;i<simStepsPerFrame;++i) {
			//reset all active nodes from last iteration	
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
	        
			for (int j=0;j<simStepFuncKeys.length;++j) {
				String key = simStepFuncKeys[j];
				cuLaunchKernel(cuFuncs.get(key), 
						simGridSz[j], 1, 1,           // Grid dimension 
		                numCUDAThreads, 1, 1,  // Block dimension
		                simSharedMemSz[j], null,           // Shared memory size and stream 
		                kernelParams.get(key), null // Kernel- and extra parameters
		            );
			}
		}

 		//JCudaDriver.cuCtxSetCurrent(pctx);
		//copy from device data to host particle position arrays
		cuMemcpyDtoH(Pointer.to(h_part_pos_x),part_pos_x, numPartsFloatSz);
		cuMemcpyDtoH(Pointer.to(h_part_pos_y),part_pos_y, numPartsFloatSz);
		cuMemcpyDtoH(Pointer.to(h_part_pos_z),part_pos_z, numPartsFloatSz);
		return false;
	}//simMe}
	
	//sim method to show execution time for each step
	public boolean simMeDebug(float modAmtMillis) {
		return false;
	}//simMeDebug
	
	
	//draw this sim's collider object
	private void drawCollider(MPM_SimMain pa) {
		pa.pushMatrix();pa.pushStyle();
		pa.sphereDetail(20);
		pa.translate(colLocation);
		pa.sphere(this.colRad);
		//not using cylinder
		//pa.drawCylinder(cylPts, pa.gui_DarkGreen, pa.gui_DarkCyan);

		pa.popStyle();pa.popMatrix();		
	}//drawCollider

	@Override
	public void drawMe(MPM_SimMain pa, float animTimeMod) {
		if(!getSimFlags(this.simIsBuiltIDX)) {return;}//if not built yet, don't try to draw anything

		pa.pushMatrix();pa.pushStyle();
		pa.sphereDetail(4);
		pa.scale(sclAmt);
		pa.strokeWeight(1.0f/sclAmt);
		
		pa.lights();
		//System.out.println("Draw");
		pa.ambientLight(102, 102, 102);
		pa.lightSpecular(204, 204, 204);
		pa.directionalLight(102, 102, 102, 0, 0, -1);
		pa.specular(255, 255, 255);

		pa.shininess(5.0f);
		pa.pushMatrix();pa.pushStyle();
		for(int i=0;i<numParts;++i) {
			pa.strokeWeight(3.0f);
			pa.stroke(h_part_clr[i][0], h_part_clr[i][1], h_part_clr[i][2]);
			pa.point(h_part_pos_x[i], h_part_pos_y[i], h_part_pos_z[i]);
			
		}
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

	
}//class MPM_BaseSim




