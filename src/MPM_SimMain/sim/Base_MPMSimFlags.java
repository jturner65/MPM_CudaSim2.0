package MPM_SimMain.sim;

import base_Utils_Objects.tools.flags.Base_BoolFlags;

public abstract class Base_MPMSimFlags extends Base_BoolFlags {
	
	protected final Base_MPMSim owner;

	public static final int
		//debug is managed by Base_BoolFlags
		simIsBuiltIDX			= _numBaseFlags,				//simulation is finished being built
		showLocColorsIDX		= _numBaseFlags + 1,			//display particles by color of their initial location
		showColliderIDX			= _numBaseFlags + 2,			//display visual rep of collider object, if one exists
		showParticlesIDX		= _numBaseFlags + 3,			//display visual rep of material points
		showParticleVelVecsIDX  = _numBaseFlags + 4,			//plot velocity arrows for each particle	
		showGridIDX				= _numBaseFlags + 5,			//plot the computational grid
		showGridVelArrowsIDX 	= _numBaseFlags + 6,			//plot velocity arrows for each gridNode
		showGridAccelArrowsIDX 	= _numBaseFlags + 7,			//plot acceleration arrows for each gridNode
		showGridMassIDX  		= _numBaseFlags + 8,			//plot variable sized spheres proportional to gridnode mass
		showActiveNodesIDX		= _numBaseFlags + 9;			//show the grid nodes influenced by each particle
	protected static final int numBaseMPMSimFlags = _numBaseFlags + 10;
	
	public Base_MPMSimFlags(Base_MPMSim _owner, int _numFlags) {
		super(_numFlags);
		owner = _owner;
	}

	public Base_MPMSimFlags(Base_MPMSimFlags _otr) {
		super(_otr);
		owner = _otr.owner;
	}
		
	/**
	 * Whether the sim is built or not
	 * @param _val
	 */
	public final boolean getSimIsBuilt() {return getFlag(simIsBuiltIDX);}
	/**
	 * Set whether the sim is built or not
	 * @param _val
	 */
	public final void setSimIsBuilt(boolean _val) {setFlag(simIsBuiltIDX, _val);}

	/**
	 * Whether the material points should be shown with color based on location or grayscale
	 * @param _val
	 */
	public final boolean getShowLocColors() {return getFlag(showLocColorsIDX);}
	/**
	 * Set whether the material points should be shown with color based on location or grayscale
	 * @param _val
	 */
	public final void setShowLocColors(boolean _val) {setFlag(showLocColorsIDX, _val);}

	/**
	 * Whether the visual rep of collider object is displayed, if one exists
	 * @param _val
	 */
	public final boolean getShowCollider() {return getFlag(showColliderIDX);}
	/**
	 * Set whether the visual rep of collider object is displayed, if one exists
	 * @param _val
	 */
	public final void setShowCollider(boolean _val) {setFlag(showColliderIDX, _val);}

	/**
	 * Whether the material points should be shown
	 * @param _val
	 */
	public final boolean getShowParticles() {return getFlag(showParticlesIDX);}
	/**
	 * Set whether the material points should be shown
	 * @param _val
	 */
	public final void setShowParticles(boolean _val) {setFlag(showParticlesIDX, _val);}

	/**
	 * Whether the material point velocities should be shown
	 * @param _val
	 */
	public final boolean getShowPartVels() {return getFlag(showParticleVelVecsIDX);}
	/**
	 * Set whether the material point velocities should be shown
	 * @param _val
	 */
	public final void setShowPartVels(boolean _val) {setFlag(showParticleVelVecsIDX, _val);}	
	/**
	 * Whether the grid should be shown
	 * @param _val
	 */
	public final boolean getShowGrid() {return getFlag(showGridIDX);}
	/**
	 * Set whether the grid should be shown
	 * @param _val
	 */
	public final void setShowGrid(boolean _val) {setFlag(showGridIDX, _val);}

	/**
	 * Whether the grid velocities should be shown
	 * @param _val
	 */
	public final boolean getShowGridVel() {return getFlag(showGridVelArrowsIDX);}
	/**
	 * Set whether the grid velocities should be shown
	 * @param _val
	 */
	public final void setShowGridVel(boolean _val) {setFlag(showGridVelArrowsIDX, _val);}

	/**
	 * Whether the grid accelerations should be shown
	 * @param _val
	 */
	public final boolean getShowGridAccel() {return getFlag(showGridAccelArrowsIDX);}
	/**
	 * Set whether the grid accelerations should be shown
	 * @param _val
	 */
	public final void setShowGridAccel(boolean _val) {setFlag(showGridAccelArrowsIDX, _val);}

	/**
	 * Whether variable sized spheres proportional to gridnode mass should be shown
	 * @param _val
	 */
	public final boolean getShowGridMass() {return getFlag(showGridMassIDX);}
	/**
	 * Set whether variable sized spheres proportional to gridnode mass should be shown
	 * @param _val
	 */
	public final void setShowGridMass(boolean _val) {setFlag(showGridMassIDX, _val);}

	/**
	 * Whether to show the grid nodes influenced by each particle
	 * @param _val
	 */
	public final boolean getShowActiveNodes() {return getFlag(showActiveNodesIDX);}
	/**
	 * Set whether to show the grid nodes influenced by each particle
	 * @param _val
	 */
	public final void setShowActiveNodes(boolean _val) {setFlag(showActiveNodesIDX, _val);}

	
	
	@Override
	protected final void handleSettingDebug(boolean val) {owner.handleSimFlagsDebug(val);	}

	/**
	 * Custom handling of flag setting
	 */
	@Override
	protected final void handleFlagSet_Indiv(int idx, boolean val, boolean oldval) {
		switch(idx){
			case simIsBuiltIDX 				: {break;}			
			case showLocColorsIDX 			: {break;}
			case showColliderIDX 			: {break;}
			case showParticlesIDX			: {break;}
			case showParticleVelVecsIDX 	: {break;}
			case showGridIDX				: {break;}				
			case showGridVelArrowsIDX 		: {break;}		
			case showGridAccelArrowsIDX 	: {break;}
			case showGridMassIDX  			: {break;}		
			case showActiveNodesIDX  		: {break;}
			default :{
				handleMPMFlagSet_Indiv(idx, val, oldval);
			}					
		}
	}
	
	/**
	 * set values for instancing class-specific boolean flags
	 * @param idx
	 * @param val
	 */
	protected abstract void handleMPMFlagSet_Indiv(int idx, boolean val, boolean oldval);
}//class Base_MPMSimFlags
