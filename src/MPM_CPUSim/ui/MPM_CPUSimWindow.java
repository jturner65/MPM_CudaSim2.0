/**
 * 
 */
package MPM_CPUSim.ui;

import java.util.ArrayList;
import java.util.TreeMap;

import MPM_CPUSim.sim.MPM_CPUSnowBall;
import MPM_SimMain.sim.Base_MPMSim;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.GUI_AppManager;

/**
 * @author John Turner
 *
 */
public class MPM_CPUSimWindow extends Base_MPMSimWindow {
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
	public MPM_CPUSimWindow(IRenderInterface _p, GUI_AppManager _AppMgr, int _winIdx) {
		super(_p, _AppMgr, _winIdx);
		super.initThisWin(false);
	}

	@Override
	protected void initDispFlags_Indiv() {	}
	@Override
	protected final void initMe_Indiv() {
		privFlags.setFlag(showCollider, true);
	}
	
	@Override
	protected int initAllMPMPrivBtns_Indiv(ArrayList<Object[]> tmpBtnNamesArray) {
		return numPrivFlags;
	}
	
	@Override
	protected Base_MPMSim buildSim() {
		return new MPM_CPUSnowBall(ri, this, (MPM_SimUpdateFromUIData) getUIDataUpdater());
	}

	@Override
	protected void setPrivFlags_Indiv(int idx, boolean val) {
		switch (idx) {//special actions for each flag
			default						: {return;}
		}
	}

	/**
	 * Instancing class-specific (application driven) UI objects should be defined
	 * in this function.  Add an entry to tmpBtnNamesArray for each button, in the order 
	 * they are to be displayed
	 * @param tmpUIObjArray array list of Object arrays, where in each object array : 
	 * 			the first element double array of min/max/mod values
	 * 			the 2nd element is starting value
	 * 			the 3rd elem is label for object
	 * 			the 4th element is boolean array of {treat as int, has list values, value is sent to owning window}
	 * @param tmpListObjVals treemap keyed by object IDX and value is list of strings of values for all UI list select objects
	 */
	@Override
	protected void setupGUIObjsAras_Indiv(TreeMap<Integer, Object[]> tmpUIObjArray,
			TreeMap<Integer, String[]> tmpListObjVals) {

	}

	@Override
	protected final boolean setUI_IntValsCustom_Indiv(int UIidx, int ival, int oldVal) {
		return false;
	}

	@Override
	protected final boolean setUI_FloatValsCustom_Indiv(int UIidx, float ival, float oldVal) {
		return false;
	}

}//class MPM_CPUSimWindow
