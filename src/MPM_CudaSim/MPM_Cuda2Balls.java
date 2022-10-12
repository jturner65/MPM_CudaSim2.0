package MPM_CudaSim;

import java.util.ArrayList;
import java.util.TreeMap;

import MPM_CudaSim.base.base_MPMCudaSim;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;

//instance of sim world with 2 big snow boulders slamming into each other 
public class MPM_Cuda2Balls extends base_MPMCudaSim {
	//scale w/timestep
	private static float initVel = 30.0f;
	
	public MPM_Cuda2Balls(IRenderInterface _pa, MPM_SimWindow _win, int _gridCount, float _h, int _numParts, float _density) {
		super(_pa,_win,"2 Big Snowballs",_gridCount, _h,_numParts, _density);
	}	
	
	//build particle layout for cuda sim - use multiples of h as radius
	@Override
	protected void buildPartLayoutMap(TreeMap<String, ArrayList<Float[]>> partVals) {	
        int numPartsPerSphere = numParts/2;

		float offScl = 600.0f/this.sclAmt;

        float xVel = -1.0f*initVel, 
        	yVel = 0, 
    		zVel = .2f*-initVel;
		
        float[] sphere1_Ctr = new float[] {.4f  *offScl, .5f*offScl, .4f *offScl};
		
		//build sphere of particles - scale volume of sphere based on cuberoot of # of particles, with 1000 particles being baseline sphere - radius will be function of how many particles are built
        float ballRad = (float) (3.0*Math.cbrt(numParts)/sclAmt);		

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
