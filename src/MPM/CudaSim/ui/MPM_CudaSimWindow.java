/**
 * 
 */
package MPM.CudaSim.ui;

import java.util.TreeMap;

import MPM.BaseSim.sim.Base_MPMSim;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import MPM.CudaSim.sim.MPM_CudaSnowBalls;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.GUI_AppManager;

/**
 * @author John Turner
 *
 */
public class MPM_CudaSimWindow extends Base_MPMSimWindow {
	/**
	 * # of private boolean flags for this window - expands upon those determined in SOM_AnimWorldWin
	 */
	private final int numPrivFlags = numBaseMPMWinUIPrivFlags;	
	
	/**
	 * @param _p
	 * @param _AppMgr
	 * @param _winIdx
	 * @param _flagIdx
	 */
	public MPM_CudaSimWindow(IRenderInterface _p, GUI_AppManager _AppMgr, int _winIdx) {
		super(_p, _AppMgr, _winIdx);
		super.initThisWin(false);
	}
	
	@Override
	protected final Base_MPMSim buildSim() {
		return new MPM_CudaSnowBalls(ri, this, (MPM_SimUpdateFromUIData) getUIDataUpdater());
	}
	
	@Override
	protected double[] getMinMaxModParts() {return new double[]{1000, 1000000, 10000};}

	@Override
	protected double getInitNumParts() {return 1.0*initNumParts;}
	@Override
	protected void initDispFlags_Indiv() {}
	@Override
	protected final void initMe_Indiv() {
		//properly initialize gui objs
	}
	
	/**
	 * Build all UI objects to be shown in left side bar menu for this window.  This is the first child class function called by initThisWin
	 * @param tmpUIObjArray : map of object data, keyed by UI object idx, with array values being :                    
	 *           the first element double array of min/max/mod values                                                   
	 *           the 2nd element is starting value                                                                      
	 *           the 3rd elem is label for object                                                                       
	 *           the 4th element is object type (GUIObj_Type enum)
	 *           the 5th element is boolean array of : (unspecified values default to false)
	 *           	idx 0: value is sent to owning window,  
	 *           	idx 1: value is sent on any modifications (while being modified, not just on release), 
	 *           	idx 2: changes to value must be explicitly sent to consumer (are not automatically sent),
	 *           the 6th element is a boolean array of format values :(unspecified values default to false)
	 *           	idx 0: whether multi-line(stacked) or not                                                  
	 *              idx 1: if true, build prefix ornament                                                      
	 *              idx 2: if true and prefix ornament is built, make it the same color as the text fill color.
	 * @param tmpListObjVals : map of string arrays, keyed by UI object idx, with array values being each element in the list
	 * @param firstBtnIDX : first index to place button objects in @tmpBtnNamesArray 
	 * @param tmpBtnNamesArray : map of Object arrays to be built containing all button definitions, keyed by sequential value == objId
	 * 				the first element is true label
	 * 				the second element is false label
	 * 				the third element is integer flag idx 
	 */
	@Override
	protected void setupGUIObjsAras_Indiv(TreeMap<Integer, Object[]> tmpUIObjArray,	TreeMap<Integer, String[]> tmpListObjVals, int firstBtnIDX, TreeMap<Integer, Object[]> tmpBtnNamesArray) {}
	
	/**
	 * Retrieve the total number of defined privFlags booleans (application-specific state bools and interactive buttons)
	 */
	@Override
	public int getTotalNumOfPrivBools() { return numPrivFlags;}

	@Override
	protected final void setPrivFlags_Indiv(int idx, boolean val) {
		switch (idx) {//special actions for each flag
			default						: {return;}
		}
	}


	@Override
	protected final boolean setUI_IntValsCustom_Indiv(int UIidx, int ival, int oldVal) {
		return false;
	}

	@Override
	protected final boolean setUI_FloatValsCustom_Indiv(int UIidx, float ival, float oldVal) {
		return false;
	}

}//class MPM_CudaSimWindow
