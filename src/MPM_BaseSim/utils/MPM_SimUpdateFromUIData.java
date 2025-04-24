package MPM_BaseSim.utils;

import java.util.HashMap;
import java.util.Map;

import MPM_BaseSim.ui.Base_MPMSimWindow;
import base_UI_Objects.windowUI.uiData.UIDataUpdater;

public class MPM_SimUpdateFromUIData extends UIDataUpdater {
	private static final int 
		RemakeKernelIDX = 0,
		ResetSimIDX = 1,
		RebuildSimIDX = 2,
		MatHasChangedIDX = 3;
	private static final int numStates = 4;
	
	private static final int
		intIDX = 0,
		floatIDX = 1,
		boolIDX = 2;
	//int, float and bool
	private static final int numTypes = 3;
	
	@SuppressWarnings("unchecked")
	private static HashMap<Integer, Integer>[][] simStateIDXsToCheck = new HashMap[numStates][numTypes]; 
	static { 
		for(int i=0;i<numStates;++i) {
			@SuppressWarnings("unchecked")
			HashMap<Integer, Integer>[] tmpMapAra = new HashMap[numTypes];
			for (int j=0;j<numTypes;++j) {tmpMapAra[j] = new HashMap<Integer,Integer>();}
			simStateIDXsToCheck[i] = tmpMapAra;
		}
		//requires rebuilding of simulation environment
		simStateIDXsToCheck[RebuildSimIDX][intIDX].put(Base_MPMSimWindow.gIDX_GridCount, Base_MPMSimWindow.gIDX_GridCount);
		simStateIDXsToCheck[RebuildSimIDX][intIDX].put(Base_MPMSimWindow.gIDX_NumSnowballs, Base_MPMSimWindow.gIDX_NumSnowballs);
		simStateIDXsToCheck[RebuildSimIDX][floatIDX].put(Base_MPMSimWindow.gIDX_GridCellSize, Base_MPMSimWindow.gIDX_GridCellSize);
		//requires only resetting existing sim
		simStateIDXsToCheck[ResetSimIDX][intIDX].put(Base_MPMSimWindow.gIDX_NumParticles, Base_MPMSimWindow.gIDX_NumParticles);
		simStateIDXsToCheck[ResetSimIDX][floatIDX].put(Base_MPMSimWindow.gIDX_InitVel, Base_MPMSimWindow.gIDX_InitVel);
		//requires rebuilding kernel/solver
		simStateIDXsToCheck[RemakeKernelIDX][floatIDX].put(Base_MPMSimWindow.gIDX_TimeStep, Base_MPMSimWindow.gIDX_TimeStep);		
		simStateIDXsToCheck[RemakeKernelIDX][floatIDX].put(Base_MPMSimWindow.gIDX_PartMass, Base_MPMSimWindow.gIDX_PartMass);		
		simStateIDXsToCheck[RemakeKernelIDX][floatIDX].put(Base_MPMSimWindow.gIDX_wallFricCoeff, Base_MPMSimWindow.gIDX_wallFricCoeff);		
		simStateIDXsToCheck[RemakeKernelIDX][floatIDX].put(Base_MPMSimWindow.gIDX_CollFricCoeff, Base_MPMSimWindow.gIDX_CollFricCoeff);		
		//Materials have changed, also requires rebuild
		simStateIDXsToCheck[MatHasChangedIDX][floatIDX].put(Base_MPMSimWindow.gIDX_InitYoungMod, Base_MPMSimWindow.gIDX_InitYoungMod);		 
		simStateIDXsToCheck[MatHasChangedIDX][floatIDX].put(Base_MPMSimWindow.gIDX_PoissonRatio, Base_MPMSimWindow.gIDX_PoissonRatio); 			 
		simStateIDXsToCheck[MatHasChangedIDX][floatIDX].put(Base_MPMSimWindow.gIDX_HardeningCoeff, Base_MPMSimWindow.gIDX_HardeningCoeff);		 
		simStateIDXsToCheck[MatHasChangedIDX][floatIDX].put(Base_MPMSimWindow.gIDX_CriticalCompression, Base_MPMSimWindow.gIDX_CriticalCompression);
		simStateIDXsToCheck[MatHasChangedIDX][floatIDX].put(Base_MPMSimWindow.gIDX_CriticalStretch, Base_MPMSimWindow.gIDX_CriticalStretch);    
		simStateIDXsToCheck[MatHasChangedIDX][floatIDX].put(Base_MPMSimWindow.gIDX_AlphaPicFlip, Base_MPMSimWindow.gIDX_AlphaPicFlip);	   
	}
	
	public MPM_SimUpdateFromUIData(Base_MPMSimWindow _win) {		super(_win);	}
	public MPM_SimUpdateFromUIData(Base_MPMSimWindow _win, Map<Integer, Integer> _iVals, Map<Integer, Float> _fVals,
			Map<Integer, Boolean> _bVals) {
		super(_win, _iVals, _fVals, _bVals);
	}
	public MPM_SimUpdateFromUIData(MPM_SimUpdateFromUIData _otr) {		super(_otr);	}

	/**
	 * access app-specific ints
	 */
	public int getSimStepsPerFrame() {return intValues.get(Base_MPMSimWindow.gIDX_SimStepsPerFrame);}
	public int getNumParticles() {return intValues.get(Base_MPMSimWindow.gIDX_NumParticles);}
	public int getNumSnowballs() {return intValues.get(Base_MPMSimWindow.gIDX_NumSnowballs);}
	public int getGridCellsPerSide() {return intValues.get(Base_MPMSimWindow.gIDX_GridCount);}
	public int getDrawPtIncr() {return intValues.get(Base_MPMSimWindow.gIDX_DrawPointIncr);}
	
	/**
	 * access app-specific floats
	 */
	//sim related
	public float getTimeStep() {return floatValues.get(Base_MPMSimWindow.gIDX_TimeStep);}		 	          
	public float getPartMass() {return floatValues.get(Base_MPMSimWindow.gIDX_PartMass);}          
	public float getInitVel() {return floatValues.get(Base_MPMSimWindow.gIDX_InitVel);}
	public float getGridCellSize() {return floatValues.get(Base_MPMSimWindow.gIDX_GridCellSize);}
	public float getWallFricCoeff() {return floatValues.get(Base_MPMSimWindow.gIDX_wallFricCoeff);}		     
	public float getCollFricCoeff() {return floatValues.get(Base_MPMSimWindow.gIDX_CollFricCoeff);} 
	public float getDrawnVecScale() {return floatValues.get(Base_MPMSimWindow.gIDX_DrawnValScale);} 
	

	//material related
	public float getInitYoungMod() {return floatValues.get(Base_MPMSimWindow.gIDX_InitYoungMod);} 			      
	public float getPoissonRatio() {return floatValues.get(Base_MPMSimWindow.gIDX_PoissonRatio);} 			      
	public float getHardeningCoeff() {return floatValues.get(Base_MPMSimWindow.gIDX_HardeningCoeff);}		     
	public float getCriticalCompression() {return floatValues.get(Base_MPMSimWindow.gIDX_CriticalCompression);} 	 
	public float getCriticalStretch() {return floatValues.get(Base_MPMSimWindow.gIDX_CriticalStretch);}    
	public float getAlphaPicFlip() {return floatValues.get(Base_MPMSimWindow.gIDX_AlphaPicFlip);}	      

	/**
	 * Utilities
	 */	
	/**
	 * Check whether simulation material values have changed via UI input
	 * @param _otr
	 * @return
	 */
	public boolean haveMaterialValsChanged(MPM_SimUpdateFromUIData _otr) {   
		return havePassedValuesChanged(_otr,
				simStateIDXsToCheck[MatHasChangedIDX][intIDX],
				simStateIDXsToCheck[MatHasChangedIDX][floatIDX],
				simStateIDXsToCheck[MatHasChangedIDX][boolIDX]);
	}

	/**
	 * Check whether the simulation should be rebuilt without remaking solver kernel
	 * @param _otr
	 * @return
	 */
	public boolean checkSimRebuild(MPM_SimUpdateFromUIData _otr) {
		return havePassedValuesChanged(_otr,
				simStateIDXsToCheck[RebuildSimIDX][intIDX],
				simStateIDXsToCheck[RebuildSimIDX][floatIDX],
				simStateIDXsToCheck[RebuildSimIDX][boolIDX]);
	}
	
	/**
	 * Check whether simulation should only be reset
	 * @param _otr
	 * @return 
	 */
	public boolean checkSimReset(MPM_SimUpdateFromUIData _otr) {
		return havePassedValuesChanged(_otr,
				simStateIDXsToCheck[ResetSimIDX][intIDX],
				simStateIDXsToCheck[ResetSimIDX][floatIDX],
				simStateIDXsToCheck[ResetSimIDX][boolIDX]);
	}
	
	/**
	 * Check whether solver kernel needs to be remade
	 * @param _otr
	 * @return
	 */
	public boolean checkSimKernelRebuilt(MPM_SimUpdateFromUIData _otr) {
		return havePassedValuesChanged(_otr,
				simStateIDXsToCheck[RemakeKernelIDX][intIDX],
				simStateIDXsToCheck[RemakeKernelIDX][floatIDX],
				simStateIDXsToCheck[RemakeKernelIDX][boolIDX]);
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
		String res = "Owning Window Name: "+win.getName()+" | Tracked values : "+intValues.size() +" Integers, " +floatValues.size() +" Floats, " +boolValues.size() + " Booleans\n";
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
