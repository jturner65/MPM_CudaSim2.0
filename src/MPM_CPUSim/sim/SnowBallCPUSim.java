/**
 * 
 */
package MPM_CPUSim.sim;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Render_Interface.IRenderInterface;

/**
 * @author John Turner
 *
 */
public class SnowBallCPUSim extends Base_MPMCPUSim {

	/**
	 * @param _pa
	 * @param _win
	 * @param _simName
	 * @param _gravity
	 * @param _currUIVals
	 */
	public SnowBallCPUSim(IRenderInterface _pa, Base_MPMSimWindow _win, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa, _win, "Snowball Fall", new float[] {0.0f, -9.8f, 0.0f}, _currUIVals);
	}

	@Override
	protected void updateCPUSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd) {
	}//updateCPUSimVals_FromUI_Indiv

	@Override
	protected void drawColliders_Indiv(float animTimeMod) {}

	/**
	 * sim method to show execution time and debug information for each sim step
	 * @param modAmtMillis
	 * @return
	 */
	@Override
	public final boolean simMeDebug(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}
}//class SnowBallSim
