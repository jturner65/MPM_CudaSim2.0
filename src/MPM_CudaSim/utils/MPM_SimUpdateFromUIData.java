package MPM_CudaSim.utils;

import java.util.HashMap;
import java.util.Map;

import MPM_CudaSim.ui.MPM_SimWindow;
import base_UI_Objects.windowUI.uiData.UIDataUpdater;

public class MPM_SimUpdateFromUIData extends UIDataUpdater {

	public MPM_SimUpdateFromUIData(MPM_SimWindow _win) {		super(_win);	}
	public MPM_SimUpdateFromUIData(MPM_SimWindow _win, Map<Integer, Integer> _iVals, Map<Integer, Float> _fVals,
			Map<Integer, Boolean> _bVals) {
		super(_win, _iVals, _fVals, _bVals);
	}
	public MPM_SimUpdateFromUIData(MPM_SimUpdateFromUIData _otr) {		super(_otr);	}

	/**
	 * access app-specific ints
	 */
	public int getSimStepsPerFrame() {return intValues.get(MPM_SimWindow.gIDX_SimStepsPerFrame);}
	public int getNumParticles() {return intValues.get(MPM_SimWindow.gIDX_NumParticles);}
	public int getNumSnowballs() {return intValues.get(MPM_SimWindow.gIDX_NumSnowballs);}
	public int getGridCellsPerSide() {return intValues.get(MPM_SimWindow.gIDX_GridCount);}
	public int getDrawPtIncr() {return intValues.get(MPM_SimWindow.gIDX_DrawPointIncr);}
	
	/**
	 * access app-specific floats
	 */
	//sim related
	public float getTimeStep() {return floatValues.get(MPM_SimWindow.gIDX_TimeStep);}		 	          
	public float getPartMass() {return floatValues.get(MPM_SimWindow.gIDX_PartMass);}          
	public float getInitVel() {return floatValues.get(MPM_SimWindow.gIDX_InitVel);}
	public float getGridCellSize() {return floatValues.get(MPM_SimWindow.gIDX_GridCellSize);}
	public float getWallFricCoeff() {return floatValues.get(MPM_SimWindow.gIDX_wallFricCoeff);}		     
	public float getCollFricCoeff() {return floatValues.get(MPM_SimWindow.gIDX_CollFricCoeff);} 
	public float getDrawnVecScale() {return floatValues.get(MPM_SimWindow.gIDX_DrawnValScale);} 
	

	//material related
	public float getInitYoungMod() {return floatValues.get(MPM_SimWindow.gIDX_InitYoungMod);} 			      
	public float getPoissonRatio() {return floatValues.get(MPM_SimWindow.gIDX_PoissonRatio);} 			      
	public float getHardeningCoeff() {return floatValues.get(MPM_SimWindow.gIDX_HardeningCoeff);}		     
	public float getCriticalCompression() {return floatValues.get(MPM_SimWindow.gIDX_CriticalCompression);} 	 
	public float getCriticalStretch() {return floatValues.get(MPM_SimWindow.gIDX_CriticalStretch);}    
	public float getAlphaPicFlip() {return floatValues.get(MPM_SimWindow.gIDX_AlphaPicFlip);}	      
		
	/**
	 * access app-specific booleans
	 */		
	//TODO
	
	
	/**
	 * Utilities
	 */
	
	public boolean haveMaterialValsChanged(MPM_SimUpdateFromUIData _otr) {
		HashMap<Integer,Integer> FloatIdxsToCheck = new HashMap<Integer,Integer>();		
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_InitYoungMod, MPM_SimWindow.gIDX_InitYoungMod);		 
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_PoissonRatio,MPM_SimWindow.gIDX_PoissonRatio); 			 
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_HardeningCoeff,MPM_SimWindow.gIDX_HardeningCoeff);		 
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_CriticalCompression,MPM_SimWindow.gIDX_CriticalCompression);
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_CriticalStretch,MPM_SimWindow.gIDX_CriticalStretch);    
		FloatIdxsToCheck.put(MPM_SimWindow.gIDX_AlphaPicFlip,MPM_SimWindow.gIDX_AlphaPicFlip);	   
		return havePassedValuesChanged(_otr,new HashMap<Integer,Integer>(), FloatIdxsToCheck, new HashMap<Integer,Integer>());
	}
	
	/**
	 * Whether or not modifications to UI values that the sim relies on occur, in which case the sim should be reset 
	 * @return
	 */
	public boolean shouldSimBeReset(MPM_SimUpdateFromUIData _otr, 
			HashMap<Integer,Integer> IntIdxsToIgnore, 
			HashMap<Integer,Integer> FloatIdxsToIgnore, 
			HashMap<Integer,Integer> BoolIdxsToCheck) {
		return haveValuesChangedExceptPassed(_otr,IntIdxsToIgnore, FloatIdxsToIgnore, BoolIdxsToCheck);
	}
	
	/**
	 * 
	 * @param _otr
	 * @return
	 */
	public boolean checkPassedValuesChanged(MPM_SimUpdateFromUIData _otr, 
			HashMap<Integer,Integer> IntIdxsToCheck, 
			HashMap<Integer,Integer> FloatIdxsToCheck, 
			HashMap<Integer,Integer> BoolIdxsToCheck){
		return havePassedValuesChanged(_otr,IntIdxsToCheck, FloatIdxsToCheck, BoolIdxsToCheck);
	}
	
	@Override
	public String toString() {
		String res = "Owning Window Name: "+win.name+" | Tracked values : "+intValues.size() +" Integers, " +floatValues.size() +" Floats, " +boolValues.size() + " Booleans\n";
		if (intValues.size() > 0) {
			res+="Int Values: (" +intValues.size() +")\n";
			res+="\tSim Steps per Frame : "+getSimStepsPerFrame() +"\n";
			res+="\tNum Particles : "+getNumParticles() +"\n";
			res+="\tNum Snowballs : "+getNumSnowballs() +"\n";
			res+="\tNum Grid Cells Per Side : "+getGridCellsPerSide() +"\n";	
		} else {		res+="No Integer values present/being tracked";}
		
		if (floatValues.size() > 0) {
			res+="Float Values: (" +floatValues.size() +")\n";
			res+= "\tTime Step : "+getTimeStep() +"\n"; 	 	          
			res+= "\tParticle Mass : "+getPartMass()+"\n";         
			res+= "\tGrid Cell Size per Dim :"+getGridCellSize() +"\n"; 
			res+= "\tInit Young's Modulus : "+getInitYoungMod() +"\n";  			      
			res+= "\tPoisson Ratio : "+getPoissonRatio() +"\n";  			      
			res+= "\tHardening Coefficient : "+getHardeningCoeff() +"\n"; 		     
			res+= "\tCritical Compression : "+getCriticalCompression() +"\n";  
			res+= "\tCritical Stretch : "+getCriticalStretch() +"\n";    
			res+= "\tAlpha PIC/FLIP : "+getAlphaPicFlip() +"\n";     
			res+= "\tWall Friction Coefficient : "+getWallFricCoeff()+"\n"; 		     
			res+= "\tCollider Friction Coefficient : "+getCollFricCoeff()+"\n";
		} else {		res+="No Float values present/being tracked";	}
		
		if(boolValues.size() > 0) {
			res+="Boolean Values: (" +boolValues.size() +")\n";
			for (Map.Entry<Integer, Boolean> entry : boolValues.entrySet()) {
				res += "\tKey : "+entry.getKey()+" | Value : "+entry.getValue()+"\n";
			}		
		} else {		res+="No Boolean values present/being tracked";	}
		return res;
	}//toString	
	
}//class MPM_SimUpdateFromUIData
