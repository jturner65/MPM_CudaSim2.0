package MPM_CudaSim;

import java.util.ArrayList;
import java.util.TreeMap;

import MPM_CudaSim.base.base_MPMCudaSim;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;

//instance of sim world with 2 big snow boulders slamming into each other 
public class MPM_Cuda2Balls extends base_MPMCudaSim {
	//scale w/timestep
	private static float initVel = 30.0f;
	
	public MPM_Cuda2Balls(IRenderInterface _pa,int _gridCount, float _h, int _numParts) {
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
