package MPM.CudaSim.sim.base;

import MPM.BaseSim.sim.Base_MPMSimFlags;

public class MPM_CudaSimFlags extends Base_MPMSimFlags {
	
	public static final int
		CUDADevInitIDX				= numBaseMPMSimFlags + 0;			//if 1 time cuda device and kernel file init is complete
	protected static final int numCudaSimFlags = numBaseMPMSimFlags + 1;


	public MPM_CudaSimFlags(Base_MPMCudaSim _owner) {
		super(_owner, numCudaSimFlags);
	}

	public MPM_CudaSimFlags(MPM_CudaSimFlags _otr) {
		super(_otr);
	}
		
	/**
	 * Whether the 1 time cuda device and kernel file init is complete
	 * @param _val
	 */
	public final boolean getCudaDevInit() {return getFlag(CUDADevInitIDX);}
	/**
	 * Set whether the 1 time cuda device and kernel file init is complete
	 * @param _val
	 */
	public final void setCudaDevInit(boolean _val) {setFlag(CUDADevInitIDX, _val);}

	@Override
	protected final void handleMPMFlagSet_Indiv(int idx, boolean val, boolean oldval) {
		switch(idx) {
			case CUDADevInitIDX : {break;}
			default : {
				handleCudaMPMFlagSet_Indiv(idx, val, oldval);
			}
		}
	}//handleMPMFlagSet_Indiv
	
	/**
	 * Override this function if instancing class adds flags
	 * @param idx
	 * @param val
	 * @param oldval
	 */
	protected void handleCudaMPMFlagSet_Indiv(int idx, boolean val, boolean oldval) {}
	
}//class MPMCudaSimFlags
