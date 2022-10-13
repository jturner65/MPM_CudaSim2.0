package MPM_CudaSim.utils;

import java.util.HashMap;
import java.util.Map;

import MPM_CudaSim.ui.MPM_SimWindow;
import base_UI_Objects.windowUI.baseUI.base_UpdateFromUIData;

public class MPM_SimUpdateFromUIData extends base_UpdateFromUIData {
	
	protected int[] IntIDXsToCompare = new int[] {MPM_SimWindow.gIDX_NumParticles, MPM_SimWindow.gIDX_GridCount};
	protected int[] FloatIDXsToCompare = new int[] {
			MPM_SimWindow.gIDX_TimeStep, MPM_SimWindow.gIDX_PartMass,
			MPM_SimWindow.gIDX_GridCellSize, MPM_SimWindow.gIDX_wallFricCoeff,
			MPM_SimWindow.gIDX_CollFricCoeff, MPM_SimWindow.gIDX_InitYoungMod,
			MPM_SimWindow.gIDX_PoissonRatio, MPM_SimWindow.gIDX_HardeningCoeff,
			MPM_SimWindow.gIDX_CriticalCompression, MPM_SimWindow.gIDX_CriticalStretch,
			MPM_SimWindow.gIDX_AlphaPicFlip
	};
	
	protected int[] BoolIDXsToCompare = new int[] {};

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
	public int getGridCellsPerSide() {return intValues.get(MPM_SimWindow.gIDX_GridCount);}
	
	/**
	 * access app-specific floats
	 */
	//sim related
	public float getTimeStep() {return floatValues.get(MPM_SimWindow.gIDX_TimeStep);}		 	          
	public float getPartMass() {return floatValues.get(MPM_SimWindow.gIDX_PartMass);}          
	public float getGridCellSize() {return floatValues.get(MPM_SimWindow.gIDX_GridCellSize);}
	public float getWallFricCoeff() {return floatValues.get(MPM_SimWindow.gIDX_wallFricCoeff);}		     
	public float getCollFricCoeff() {return floatValues.get(MPM_SimWindow.gIDX_CollFricCoeff);}       

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
	
	/**
	 * Whether or not modifications to UI values that the sim relies on occur, in which case the sim should be reset 
	 * @return
	 */
	public boolean shouldSimBeReset(MPM_SimUpdateFromUIData _otr) {
		HashMap<Integer,Integer> IntIdxsToIgnore = new HashMap<Integer,Integer>();
		IntIdxsToIgnore.put(MPM_SimWindow.gIDX_SimStepsPerFrame, MPM_SimWindow.gIDX_SimStepsPerFrame);
		return haveValuesChanged(_otr,IntIdxsToIgnore, new HashMap<Integer,Integer>(),new HashMap<Integer,Integer>());
	}
	
	/**
	 * Rebuild simulation if any simulator-dependent variables have changed. These are variables that are sent to the cuda kernel
	 * @param _otr
	 * @return
	 */
	public boolean shouldSimBeRebuilt(MPM_SimUpdateFromUIData _otr) {
		for(int idx : IntIDXsToCompare) {
			if(intValues.get(idx) != _otr.intValues.get(idx)) {return true;}
		}
		for(int idx : FloatIDXsToCompare) {
			if(floatValues.get(idx) != _otr.floatValues.get(idx)) {return true;}
		}
		for(int idx : BoolIDXsToCompare) {
			if(boolValues.get(idx) != _otr.boolValues.get(idx)) {return true;}
		}			
			
		return false;
	}
	
	@Override
	public String toString() {
		String res = "Owning Window Name: "+win.name+" | Tracked values : "+intValues.size() +" Integers, " +floatValues.size() +" Floats, " +boolValues.size() + " Booleans\n";
		if (intValues.size() > 0) {
			res+="Int Values: (" +intValues.size() +")\n";
			res+="\tSim Steps per Frame : "+getSimStepsPerFrame() +"\n";
			res+="\tNum Particles : "+getNumParticles() +"\n";
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
