package MPM_CudaSim.sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import MPM_CudaSim.sim.base.base_MPMCudaSim;
import MPM_CudaSim.ui.MPM_SimWindow;
import MPM_CudaSim.utils.MPM_SimUpdateFromUIData;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;
import base_Math_Objects.vectorObjs.floats.myVectorf;

//instance of sim world with 2 big snow boulders slamming into each other 
public class MPM_Cuda2Balls extends base_MPMCudaSim {
	//scale w/timestep
	private static float initVel = 30.0f;
	
	public MPM_Cuda2Balls(IRenderInterface _pa, MPM_SimWindow _win, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa,_win,"Snowball Slam!", _currUIVals);
	}	

	/**
	 * Specify simulation-specific IDXs of UI components to ignore changes of when determining 
	 * whether or not to rebuild simulation based on UI changes
	 * @param IntIdxsToIgnore [Out] IDXs to Integer UI components to ignore changes
	 * @param FloatIdxsToIgnore [Out] IDXs to  UI components to ignore changes
	 * @param BoolIdxsToIgnore [Out] IDXs to Integer UI components to ignore changes
	 */
	@Override
	protected void setUIIdxsToIgnorePerSim(HashMap<Integer, Integer> IntIdxsToIgnore,
			HashMap<Integer, Integer> FloatIdxsToIgnore, HashMap<Integer, Integer> BoolIdxsToIgnore) {
	}//setUIIdxsToIgnorePerSim
	
	/**
	 * build particle layout for cuda sim - use multiples of h as radius
	 * @param partVals [OUT] map of particle locs, initial velocities and min/max vals being constructed
	 * @param numPartsRequested desired # of particles. May be off a bit from final value - ALWAYS USE SIZE OF PARTVALS FOR COUNT
	 */
	static private boolean doRand = false;
	@Override
	protected void buildPartLayoutMap(TreeMap<String, ArrayList<float[]>> partVals, int numPartsRequested) {	
		
		//build sphere of particles - scale volume of sphere based on cuberoot of # of particles, with 1000 particles being baseline sphere
		//radius will be function of how many particles are built
		float sphereRad = (float)(3.0*Math.cbrt(numPartsRequested)/sclAmt);		
		float offScl = .5f * (gridCount * cellSize);
        //max valid location for center dof without breaching collider on start
        //function of # particles requested and grid dims
        float maxDim = offScl-sphereRad;

		win.getMsgObj().dispInfoMessage("MPM_Cuda2Balls : "+simName, "buildPartLayoutMap", "Sphere Rad : " + sphereRad+ " | offScl : "+offScl);
		myVectorf[] sphere_Ctrs;
		myVectorf[] sphere_Vels;
		int numSpheres = numSnowballs;
        if (!doRand){
        	//first run default setup
        	numSpheres = 2;
        	sphere_Ctrs = new myVectorf[numSpheres];
        	sphere_Vels = new myVectorf[numSpheres];        	
        	sphere_Ctrs[0] = new myVectorf(.5f*offScl, .5f*offScl, .3f*offScl);
        	sphere_Ctrs[1] = new myVectorf(-.3f*offScl, .5f*offScl, .6f* offScl);
			sphere_Vels[0] = myVectorf._sub(sphere_Ctrs[1],sphere_Ctrs[0]);
			sphere_Vels[0]._mult(initVel/sphere_Vels[0].magn);
			sphere_Vels[0].z *=.5;
			sphere_Vels[1] = myVectorf._mult(sphere_Vels[0], -1.0f);
			sphere_Vels[1].z *=.5f;
			doRand = true;
       
        } else {
        	sphere_Ctrs = new myVectorf[numSpheres];
        	sphere_Vels = new myVectorf[numSpheres];        	       	
	        for (int i=0;i<sphere_Ctrs.length;++i) {
				sphere_Ctrs[i] = getRandSphereCenter(maxDim);
			}
	        
	        myVectorf ctrTarget = new myVectorf();
	        for (myVectorf vec : sphere_Ctrs) {
	        	ctrTarget._add(vec);
	        }
	        //absolute target - vary individual target by radius amt to give glancing blows
	        ctrTarget._div(sphere_Ctrs.length);
	        
	        float[] ctrTargetAra = ctrTarget.asArray(); 
	        float largestVelMag = -1.0f;
	        float tarRad = sphereRad*1.25f;
	        for (int i=0;i<sphere_Ctrs.length;++i) {
	        	//pick a custom target for the ball to point at for glancing blow potential
	        	float[] custCtrTarget = getRandPosInSphereAra(tarRad, ctrTargetAra);
	        	sphere_Vels[i] = new myVectorf(sphere_Ctrs[i],new myVectorf(custCtrTarget[0],custCtrTarget[1],custCtrTarget[2]));
	        	if(largestVelMag < sphere_Vels[i].magn) {
	        		largestVelMag = sphere_Vels[i].magn;
	        	}
	        }
	        for (int i=0;i<sphere_Ctrs.length;++i) {
	        	sphere_Vels[i]._mult(initVel/largestVelMag);
	        	
	        }   
        }

		int numPartsPerSphere = numPartsRequested/numSpheres; 
		int numPartsLeftOver =  numPartsRequested % numSpheres; 
        for (int i=0;i<sphere_Ctrs.length-1;++i) {
        	createSphere(partVals, sphereRad, numPartsPerSphere, sphere_Vels[i], sphere_Ctrs[i].asArray());
        }
        int i = sphere_Ctrs.length-1;
        createSphere(partVals, sphereRad, numPartsPerSphere + numPartsLeftOver, sphere_Vels[i], sphere_Ctrs[i].asArray());
//        //lower ball
//        //createSphere(partVals, sphereRad, numPartsPerSphere, new float [] {xVel, yVel, zVel}, sphere1_Ctr);
//        createSphere(partVals, sphereRad, numPartsPerSphere, sphere1_Vel, sphere1_Ctr.asArray());
//		
//
//        //upper ball        
//		//createSphere(partVals, sphereRad, numPartsPerSphere, new float [] {xVel, yVel, zVel}, sphere2_Ctr);
//		createSphere(partVals, sphereRad, numPartsPerSphere, sphere2_Vel, sphere2_Ctr.asArray());

         //end create particle layout	
	}//buildPartLayout

	
	
	
	@Override
	//draw scene-specific collider, if it exists
	protected void drawCollider(float animTimeMod) {}
	
	//sim method to show execution time for each step
	public boolean simMeDebug(float modAmtMillis) {
		//TODO any debugging that might be supportable here
		return false;
	}//simMeDebug
}//class MPM_Cuda2Balls
