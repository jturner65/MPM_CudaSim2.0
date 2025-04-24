package MPM.CPUSim.sim.base;

import MPM.BaseSim.sim.Base_MPMSimFlags;

/**
 * @author John Turner
 *
 */
public class MPM_CPUSimFlags extends Base_MPMSimFlags {

	/**
	 * @param _owner
	 * @param _numFlags
	 */
	public MPM_CPUSimFlags(Base_MPMCPUSim _owner) {
		super(_owner, Base_MPMSimFlags.numBaseMPMSimFlags);
	}

	/**
	 * @param _otr
	 */
	public MPM_CPUSimFlags(Base_MPMSimFlags _otr) {
		super(_otr);
	}

	@Override
	protected void handleMPMFlagSet_Indiv(int idx, boolean val, boolean oldval) {

	}

}//class MPM_CPUSimFlags
