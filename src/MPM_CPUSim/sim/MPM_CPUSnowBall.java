/**
 * 
 */
package MPM_CPUSim.sim;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;
import MPM_CPUSim.sim.base.MPM_CPUSimFlags;
import MPM_SimMain.sim.Base_MPMSimFlags;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IRenderInterface;

/**
 * @author John Turner
 *
 */
public class MPM_CPUSnowBall extends Base_MPMCPUSim {

	/**
	 * @param _pa
	 * @param _win
	 * @param _simName
	 * @param _gravity
	 * @param _currUIVals
	 */
	public MPM_CPUSnowBall(IRenderInterface _pa, Base_MPMSimWindow _win, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa, _win, "Snowball Fall", _currUIVals);
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

	@Override
	protected Base_MPMSimFlags buildSimFlags() {
		// TODO Auto-generated method stub
		return new MPM_CPUSimFlags(this);
	}

	@Override
	public void handleSimFlagsDebug_Indiv(boolean val) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public myVectorf checkColliderCollision(myVectorf pos) {
		// TODO Auto-generated method stub
		return pos;
	}
}//class SnowBallSim
