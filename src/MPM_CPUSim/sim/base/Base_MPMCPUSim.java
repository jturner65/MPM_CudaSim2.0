package MPM_CPUSim.sim.base;


import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.jblas.FloatMatrix;

import MPM_CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM_CPUSim.sim.grid.MPM_CPUGridNode;
import MPM_CPUSim.sim.particles.MPM_CPURndrdPart;
import MPM_CPUSim.sim.threads.MPM_CPUGridBuilder;
import MPM_CPUSim.sim.threads.MPM_CPUPartBuilder;
import MPM_SimMain.sim.Base_MPMSim;
import MPM_SimMain.sim.SimResetProcess;
import MPM_SimMain.ui.Base_MPMSimWindow;
import MPM_SimMain.utils.MPM_SimUpdateFromUIData;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IRenderInterface;

/**
 * Base class for CPU-based multi-threaded MPM solver/sim
 * @author John Turner
 */
public abstract class Base_MPMCPUSim extends Base_MPMSim {
	
	//for multithreading 
	public ExecutorService th_exec;
	public int numThreadsAvail;	
	
	//runnable to launch threads to manage particle minipulations
	public MPM_CPUPartBuilder partThdMgr;
	//runnable to build and manage grid manipulations
	public MPM_CPUGridBuilder gridThdMgr;
	
	//all particles representing the material in this simulation
	public MPM_CPURndrdPart[] parts;
	//all grid nodes
	public MPM_CPUGridNode[][][] grid;	
	
	//friction coefficients of colliders
	public static float wallFric = 1.0f, collFric = 1.0f;

	//gravity 
	public final float gGravity = -9.8f;
	//const matrix for calculations
	public final FloatMatrix gravityMat = new FloatMatrix(new float[] {0, 0, gGravity});
	/**
	 * all grid nodes needing updates - using concurrent hashmap
	 * because there's no concurrent hash set in java, but the keys of CHM will work
	 */
	public ConcurrentHashMap<MPM_CPUGridNode,MPM_CPUActiveNodeAgg>[] activeNodes;
	/**
	 * @param _pa
	 * @param _win
	 * @param _simName
	 * @param _gravity
	 * @param _currUIVals
	 */
	public Base_MPMCPUSim(IRenderInterface _pa, Base_MPMSimWindow _win, String _simName, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa, _win, _simName, new float[] {0, 0, -9.8f}, _currUIVals);
		//for multithreading - do not use instanced version in PApplet - we may not use processing-based build to run simulation
		th_exec = Executors.newCachedThreadPool();	
		numThreadsAvail = win.getNumThreadsAvailable() - 2;

		//set up threads to be used to build the grid
		gridThdMgr = new MPM_CPUGridBuilder(this, numThreadsAvail);
		//set up thread builders to be used to build particles
		partThdMgr = new MPM_CPUPartBuilder(this, numThreadsAvail);		
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
	
	
	//add node to active set of nodes
	public MPM_CPUActiveNodeAgg addNodeToSet(MPM_CPUGridNode n) {	
		//build node aggregator whenever addded to set
		int nIdx = n.ID%numThreadsAvail;
		MPM_CPUActiveNodeAgg nodeAgg = activeNodes[nIdx].get(n);
		if(null == nodeAgg) {//not in map, put in map
			nodeAgg = new MPM_CPUActiveNodeAgg(n, numThreadsAvail);
			activeNodes[nIdx].put(n, nodeAgg);
		}
		return nodeAgg;
	}//addNodeToSet

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
	protected final void initValues_Parts() {}
	@Override
	protected final void initValues_Grids() {}
	@Override
	protected final boolean simMe_Indiv(float modAmtMillis) {		return false;	}
	@Override
	protected final boolean simMePost_Indiv(float modAmtMillis) {	return false;	}
	@Override
	protected final void _drawParts(float animTimeMod, boolean showLocColors) {}
	@Override
	protected final void _drawPartVel(float animTimeMod, int pincr) {}
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
	protected final void _drawGridVel(float animTimeMod) {}

	@Override
	protected final void _drawGridAccel(float animTimeMod) {}

	@Override
	protected final void _drawGridMass(float animTimeMod) {}

	/**
	 * check sim-specific central floating collider - 
	 * return 0 vec if no collision, otherwise return normal of collision
	 * @param pos
	 * @return
	 */
	public abstract myVectorf checkColliderCollision(myVectorf pos);
}//class Base_MPMCPUSim
