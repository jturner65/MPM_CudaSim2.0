package MPM.CPUSim.sim.base;

import MPM.BaseSim.sim.Base_MPMSimFlags;

/**
 * MPM flags specific for multi-threaded CPU implementation
 * @author John Turner
 *
 */
public class MPM_CPUSimFlags extends Base_MPMSimFlags {
    
    public static final int
        gridIsBuiltIDX = numBaseMPMSimFlags + 0;
    protected static final int numCPUSimFlags = numBaseMPMSimFlags + 1;
    
    /**
     * @param _owner
     * @param _numFlags
     */
    public MPM_CPUSimFlags(Base_MPMCPUSim _owner) {
        super(_owner, numCPUSimFlags);
    }

    /**
     * @param _otr
     */
    public MPM_CPUSimFlags(MPM_CPUSimFlags _otr) {
        super(_otr);
    }

    /**
     * Whether the 1 time cuda device and kernel file init is complete
     * @param _val
     */
    public final boolean getGridIsBuilt() {return getFlag(gridIsBuiltIDX);}
    /**
     * Set whether the 1 time cuda device and kernel file init is complete
     * @param _val
     */
    public final void setGridIsBuilt(boolean _val) {setFlag(gridIsBuiltIDX, _val);}

    
    @Override
    protected void handleMPMFlagSet_Indiv(int idx, boolean val, boolean oldval) {
        switch(idx) {
        default : {
            handleCPUMPMFlagSet_Indiv(idx, val, oldval);
        }
    }
}//handleMPMFlagSet_Indiv
    
    /**
     * Override this function if instancing class adds flags
     * @param idx
     * @param val
     * @param oldval
     */
    protected void handleCPUMPMFlagSet_Indiv(int idx, boolean val, boolean oldval) {}

}//class MPM_CPUSimFlags
