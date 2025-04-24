package MPM.CudaSim.sim;

import java.util.ArrayList;
import java.util.TreeMap;

import MPM.BaseSim.sim.Base_MPMSimFlags;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import MPM.CudaSim.sim.base.Base_MPMCudaSim;
import MPM.CudaSim.sim.base.MPM_CudaSimFlags;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.windowUI.base.Base_DispWindow;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;

/**
 * Instance of sim world with big snow boulders slamming into each other 
 * @author John Turner
 *
 */
public class MPM_CudaSnowBalls extends Base_MPMCudaSim {
	/**
	 * Centers of currently built spheres
	 */
	private myVectorf[] sphere_Ctrs;
	
	/**
	 * Initial Velocities of currently built spheres
	 */
	private myVectorf[] sphere_Vels;

	/**
	 * radius of sphere to make based on number of points and scale factor
	 */
	private float sphereRad;
	
	/**
	 * Abs Value of max value for any dim for sphere center (is max bnd - radius)
	 */
	private float maxCenterDim;
	
	/**
	 * Per sphere array of 2-element arrays of idxs (start, end) within partValues map for each sphere's values
	 */
	private int[][] idxsForSpheres;
	
	public MPM_CudaSnowBalls(IRenderInterface _pa, Base_MPMSimWindow _win, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa,_win,"Snowball Slam!", _currUIVals);
	}	
	
	@Override
	protected void updateCudaSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd) {
		//radius will be function of how many particles are built and grid dimensions
		sphereRad = (float)(3.0*Math.cbrt(numParts)/sclAmt);
		//derive half-length of single dimension of grid based on count of grid boxes and size of each box edge		
        //max valid location for center dof without breaching collider on start
        maxCenterDim = maxSimBnds-sphereRad;
        
        //array for start and end idx for each snowball in partVals array
        idxsForSpheres = new int[numSnowballs][];
	}//updateSimVals_FromUI_Indiv
	
	/**
	 * First instance of sim should be set to specified layout; subsequent instances are randomly synthesized
	 */
	//static public boolean doRand = false;	
	/**
	 * build particle layout for cuda sim - use multiples of h as radius
	 * @param partVals [OUT] map of particle locs, initial velocities and min/max vals being constructed
	 */
	@Override
	protected void buildPartLayoutMap(TreeMap<String, ArrayList<float[]>> partVals) {	
		win.getMsgObj().dispDebugMessage("MPM_CudaBalls("+simName+")", "buildPartLayoutMap", "Sphere Rad : " + sphereRad+ " | maxSimBnds : "+maxSimBnds);
		
       	buildSphereCtrsAndVels();
        //build point clouds and assign initial velocities
        reinitSimObjects(partVals);
        
        win.getMsgObj().dispDebugMessage("MPM_CudaBalls("+simName+")", "buildPartLayoutMap","Total # of particles in simulation : " + numParts+ " # parts float sz :"+numPartsFloatSz);
        //end create particle layout	
	}//buildPartLayout
	
	/**
	 * For when sphere counts have, or locations should, change. Rebuilds sphere centers and sphere velocities
	 */
	protected void buildSphereCtrsAndVels() {
    	sphere_Ctrs = new myVectorf[numSnowballs];
    	for (int i=0;i<sphere_Ctrs.length;++i) {
			sphere_Ctrs[i] = getRandSphereCenter(maxCenterDim);
		}
        //Find centralized location between spheres, to use as mutual target
    	myVectorf ctrTarget = new myVectorf();
        for (myVectorf vec : sphere_Ctrs) {
        	ctrTarget._add(vec);
        }
        //absolute target - vary individual target by random radius amt to give glancing blows
        ctrTarget._div(sphere_Ctrs.length);		

    	sphere_Vels = new myVectorf[sphere_Ctrs.length];           
        //float[] ctrTargetAra = ctrTarget.asArray(); 
        float tarRad = sphereRad*(1.25f + .05f*sphere_Ctrs.length);
        for (int i=0;i<sphere_Ctrs.length;++i) {
        	//pick a custom target for the ball to point at for glancing blow potential
        	//float[] custCtrTarget = getRandPosInSphereAra(tarRad, ctrTarget);
        	myPointf custCtrTarget = Base_DispWindow.AppMgr.getRandPosInSphere(tarRad, ctrTarget);
        	sphere_Vels[i] = new myVectorf(sphere_Ctrs[i], custCtrTarget);
        }   
		
	}//buildSphereCtrsAndVels

	/**
	 * Derives sample points around each sphere center. Also finds idxs in partVals map elements for each sphere
	 * @param partVals
	 */
	protected void buildSpherePoints(TreeMap<String, ArrayList<float[]>> partVals) {
        //Handle whether the requested number of parts can be evenly distributed amongst the request count of spheres
		int numPartsPerSphere = numParts/sphere_Ctrs.length; 
		int numPartsLeftOver =  numParts % sphere_Ctrs.length; 
		idxsForSpheres = new int[sphere_Ctrs.length][];  
		
		//use up extra particles, 1 per ball
        for (int i=0;i<numPartsLeftOver;++i) {
        	//add 1 to each ball for the left-over particles
        	idxsForSpheres[i] = createSphere(partVals, sphereRad, numPartsPerSphere+1, sphere_Ctrs[i]);
        }
        for (int i=numPartsLeftOver;i<sphere_Ctrs.length;++i) {
        	//make the balls without the extra points
        	idxsForSpheres[i] = createSphere(partVals, sphereRad, numPartsPerSphere, sphere_Ctrs[i]);
        }
	}//buildSpherePoints
	
	/**
	 * Derives and sets scaled initial velocities for each sphere's point cloud
	 * @param partVals
	 */
	protected void setSphereInitVelocities(TreeMap<String, ArrayList<float[]>> partVals) {
		myVectorf[] scaledVels = new myVectorf[sphere_Vels.length];
		float largestVelMag = -1.0f;
		//find largest magnitude velocity vector
		for (int i=0;i<sphere_Vels.length;++i) {
        	if(largestVelMag < sphere_Vels[i].magn) {
        		largestVelMag = sphere_Vels[i].magn;
        	}			
		}		
        float scaleVal = initVel/largestVelMag;
        for (int i=0;i<sphere_Vels.length;++i) {
        	scaledVels[i] = myVectorf._mult(sphere_Vels[i],scaleVal);
        } 			
        for (int i=0;i<idxsForSpheres.length;++i) {
        	setPartInitVelocities(partVals, idxsForSpheres[i][0],idxsForSpheres[i][1], scaledVels[i].asArray());
        }
	}//setSphereInitVelocities
	
	/**
	 * Rebuild sphere point clouds without re-deriving point centers
	 * @param partVals
	 */
	protected void reinitSimObjects(TreeMap<String, ArrayList<float[]>> partVals) {
		//use existing sphere centers to build sphere point clouds at specified sphere locations
		buildSpherePoints(partVals);
		//set initial velocities for spheres
        setSphereInitVelocities(partVals);
	}//reinitSimObjects

	/**
	 * draw scene-specific collider, if it exists
	 */
	@Override
	protected void drawColliders_Indiv(float animTimeMod) {}
	
	/**
	 * sim method to show execution time and debug information for each sim step
	 * @param modAmtMillis
	 * @return
	 */
	@Override
	public boolean simMeDebug(float modAmtMillis) {		return false;	}//simMeDebug

	@Override
	protected Base_MPMSimFlags buildSimFlags() {		return new MPM_CudaSimFlags(this);	}

	@Override
	public void handleSimFlagsDebug_Indiv(boolean val) {	}


}//class MPM_Cuda2Balls
