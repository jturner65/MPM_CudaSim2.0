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
import base_UI_Objects.windowUI.uiObjs.base.GUIObj_Params;

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
	 * @param tmpUIObjMap : map of GUIObj_Params, keyed by unique string, with values describing the UI object
	 * 			- The object IDX                   
	 *          - A double array of min/max/mod values                                                   
	 *          - The starting value                                                                      
	 *          - The label for object                                                                       
	 *          - The object type (GUIObj_Type enum)
	 *          - A boolean array of behavior configuration values : (unspecified values default to false)
	 *           	idx 0: value is sent to owning window,  
	 *           	idx 1: value is sent on any modifications (while being modified, not just on release), 
	 *           	idx 2: changes to value must be explicitly sent to consumer (are not automatically sent),
	 *          - A boolean array of renderer format values :(unspecified values default to false)
	 *           	idx 0: whether multi-line(stacked) or not                                                  
	 *              idx 1: if true, build prefix ornament                                                      
	 *              idx 2: if true and prefix ornament is built, make it the same color as the text fill color.
	 */
	protected final void setupGUIObjsAras_Indiv(TreeMap<String, GUIObj_Params> tmpUIObjMap) {}

	/**
	 * Build all UI buttons to be shown in left side bar menu for this window. This is for instancing windows to add to button region
	 * USE tmpUIBtnObjMap.size() for start idx
	 * @param tmpUIBtnObjMap : map of GUIObj_Params to be built containing all button definitions, keyed by sequential value == objId
	 * 				the first element is the object index
	 * 				the second element is true label
	 * 				the third element is false label
	 * 				the final element is integer flag idx 
	 */
	protected final void setupGUIBtnAras_Indiv(TreeMap<String, GUIObj_Params> tmpUIBtnObjMap) {}
	
	
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
