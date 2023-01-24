package MPM_CPUSim.sim.base;


import MPM_SimMain.sim.Base_MPMSim;
import MPM_SimMain.sim.SimResetProcess;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Render_Interface.IRenderInterface;

/**
 * Base class for CPU-based multi-threaded MPM solver/sim
 * @author John Turner
 */
public abstract class Base_MPMCPUSim extends Base_MPMSim {

	protected static final int numSimFlags = numBaseMPMSimFlags;
	
	/**
	 * @param _pa
	 * @param _win
	 * @param _simName
	 * @param _gravity
	 * @param _currUIVals
	 */
	public Base_MPMCPUSim(IRenderInterface _pa, Base_MPMSimWindow _win, String _simName, float[] _gravity,
			MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa, _win, _simName, _gravity, _currUIVals);
	}
	
	/**
	 * Instance-specific reset code
	 * @param rebuildSim
	 */
	@Override
	protected final void resetSim_Indiv(SimResetProcess rebuildSim) {
		
		initValues_Parts();
		
		
		initValues_Grids() ;
	}

	/**
	 * Initialize environmental layout particle holders/arrays
	 */
	protected final void initPartArrays() {		
	}
	
	/**
	 * build initial layout for particles for this simulation
	 * @param partVals [OUT] map of particle locs, initial velocities and min/max vals being constructed
	 */
	@Override
	protected final void buildPartLayouts() {
		// call instance class to build layout

	}
	/**
	 * Reinitialize existing sim - will resynthesize sample points but will not rederive locations of sim objs.
	 * @param partVals
	 */
	@Override
	protected final void reinitSimObjects() {
		// call instance class to re-init layout

	}
		
	@Override
	protected final void updateSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd) {
		//update instancing sim values
		updateCPUSimVals_FromUI_Indiv(upd);
	}//updateSimVals_FromUI_Indiv
	
	/**
	 * Update instancing class variables based on potential UI changes
	 * @param upd
	 */
	protected abstract void updateCPUSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd);

	@Override
	protected final void initValues_Parts() {
		// TODO Auto-generated method stub

	}

	@Override
	protected final void initValues_Grids() {
		// TODO Auto-generated method stub

	}

	@Override
	protected final boolean simMe_Indiv(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	protected final boolean simMePost_Indiv(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	protected final void _drawParts(float animTimeMod, boolean showLocColors) {
		// TODO Auto-generated method stub

	}

	@Override
	protected final void _drawPartVel(float animTimeMod, int pincr) {
		// TODO Auto-generated method stub

	}
	/**
	 * Draw any colliders if they exist
	 * @param animTimeMod
	 */
	@Override
	protected final void _drawColliders(float animTimeMod) {
		pa.pushMatState();	
		drawColliders_Indiv(animTimeMod);
		pa.popMatState();
	}
	/**
	 * draw internal-to-sim colliders, if they exist
	 * @param animTimeMod
	 */
	protected abstract void drawColliders_Indiv(float animTimeMod);

	@Override
	protected final void _drawGridVel(float animTimeMod) {
		// TODO Auto-generated method stub

	}

	@Override
	protected final void _drawGridAccel(float animTimeMod) {
		// TODO Auto-generated method stub

	}

	@Override
	protected final void _drawGridMass(float animTimeMod) {
		// TODO Auto-generated method stub

	}

	@Override
	protected final int getNumSimFlags() {return numSimFlags;	}

	/**
	 * set values for instancing class-specific boolean flags
	 * @param idx
	 * @param val
	 */
	@Override
	protected final void setPrivFlags_Indiv(int idx, boolean val) {}

}//class Base_MPMCPUSim
