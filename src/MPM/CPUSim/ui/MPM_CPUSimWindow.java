/**
 * 
 */
package MPM.CPUSim.ui;

import java.util.LinkedHashMap;

import MPM.BaseSim.sim.Base_MPMSim;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import MPM.CPUSim.sim.MPM_CPUSnowBall;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.uiObjs.base.GUIObj_Params;

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
    protected double[] getMinMaxModParts() {return new double[]{100, 100000, 1000};}

    @Override
    protected double getInitNumParts() {return 1000;}
    @Override
    protected void initDispFlags_Indiv() {    }
    @Override
    protected final void initMe_Indiv() {
        //Show sphere collider
        uiMgr.setPrivFlag(showCollider, true);
    }
   
    /**
     * Retrieve the total number of defined privFlags booleans (application-specific state bools and interactive buttons)
     */
    @Override
    public int getTotalNumOfPrivBools() { return numPrivFlags;}
    
    
    /**
     * Build all UI objects to be shown in left side bar menu for this window. This is the first child class function called by initThisWin
     * @param tmpUIObjMap : map of GUIObj_Params, keyed by unique string, with values describing the UI object
     *             - The object IDX                   
     *          - A double array of min/max/mod values                                                   
     *          - The starting value                                                                      
     *          - The label for object                                                                       
     *          - The object type (GUIObj_Type enum)
     *          - A boolean array of behavior configuration values : (unspecified values default to false)
     *               idx 0: value is sent to owning window,  
     *               idx 1: value is sent on any modifications (while being modified, not just on release), 
     *               idx 2: changes to value must be explicitly sent to consumer (are not automatically sent),
     *          - A boolean array of renderer format values :(unspecified values default to false) - Behavior Boolean array must also be provided!
     *                 - Should be multiline
     *                 - One object per row in UI space (i.e. default for multi-line and btn objects is false, single line non-buttons is true)
     *                 - Force this object to be on a new row/line (For side-by-side layouts)
     *                 - Text should be centered (default is false)
     *                 - Object should be rendered with outline (default for btns is true, for non-buttons is false)
     *                 - Should have ornament
     *                 - Ornament color should match label color  
     */
    @Override
    protected final void setupGUIObjsAras_Indiv(LinkedHashMap<String, GUIObj_Params> tmpUIObjMap) {}

    /**
     * Build all UI buttons to be shown in left side bar menu for this window. This is for instancing windows to add to button region
     * @param firstIdx : the first index to use in the map/as the objIdx
     * @param tmpUIBoolSwitchObjMap : map of GUIObj_Params to be built containing all flag-backed boolean switch definitions, keyed by sequential value == objId
     *                 the first element is the object index
     *                 the second element is true label
     *                 the third element is false label
     *                 the final element is integer flag idx 
     */
    @Override
    protected final void setupGUIBoolSwitchAras_Indiv(int firstIdx, LinkedHashMap<String, GUIObj_Params> tmpUIBoolSwitchObjMap)  {}
    
    
    @Override
    protected Base_MPMSim buildSim() {
        return new MPM_CPUSnowBall(ri, this, (MPM_SimUpdateFromUIData) getUIDataUpdater());
    }

    @Override
    protected void setPrivFlags_Indiv(int idx, boolean val) {
        switch (idx) {//special actions for each flag
            default                        : {return;}
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

}//class MPM_CPUSimWindow
