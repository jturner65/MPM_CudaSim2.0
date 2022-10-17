package MPM_CudaSim.sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import MPM_CudaSim.sim.base.base_MPMCudaSim;
import MPM_CudaSim.ui.MPM_SimWindow;
import MPM_CudaSim.utils.MPM_SimUpdateFromUIData;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;

//instance of sim world with 2 big snow boulders slamming into each other 
public class MPM_Cuda2Balls extends base_MPMCudaSim {
	//scale w/timestep
	private static float initVel = 30.0f;
	
	public MPM_Cuda2Balls(IRenderInterface _pa, MPM_SimWindow _win, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa,_win,"2 Big Snowballs", _currUIVals);
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
	@Override
	protected void buildPartLayoutMap(TreeMap<String, ArrayList<Float[]>> partVals, int numPartsRequested) {	
        int numPartsPerSphere = numPartsRequested/2;

		float offScl = 600.0f/this.sclAmt;

        float xVel = -1.0f*initVel, 
        	yVel = 0, 
    		zVel = .2f*-initVel;
		
        float[] sphere1_Ctr = new float[] {.4f  *offScl, .5f*offScl, .4f *offScl};
		
		//build sphere of particles - scale volume of sphere based on cuberoot of # of particles, with 1000 particles being baseline sphere - radius will be function of how many particles are built
        float ballRad = (float) (3.0*Math.cbrt(numPartsRequested)/sclAmt);		

		//int incr = 1;
		//lower ball
		float sphereRad = createSphere(partVals, ballRad, numPartsPerSphere, new float [] {xVel, yVel, zVel}, sphere1_Ctr);
        float[] sphere2_Ctr = new float[] {-.1f*sphereRad *offScl, .5f*offScl, .4f*sphereRad*offScl};
				
		xVel *= -1;
		zVel = .5f*-initVel;
        //upper ball        
		createSphere(partVals, ballRad, numPartsPerSphere, new float [] {xVel, yVel, zVel}, sphere2_Ctr);

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
