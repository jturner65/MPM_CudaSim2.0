package MPM.CPUSim.sim.base;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import MPM.BaseSim.sim.Base_MPMSim;
import MPM.BaseSim.sim.SimResetProcess;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import MPM.CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM.CPUSim.sim.grid.MPM_CPUGridNode;
import MPM.CPUSim.sim.particles.MPM_CPURndrdPart;
import MPM.CPUSim.sim.threads.MPM_CPUGridBuilder;
import MPM.CPUSim.sim.threads.MPM_CPUPartBuilder;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IRenderInterface;

/**
 * Base class for CPU-based multi-threaded MPM solver/sim
 * @author John Turner
 */
public abstract class Base_MPMCPUSim extends Base_MPMSim {
	
	/**
	 * for multithreading 
	 */
	public ExecutorService th_exec;
	public int numThreadsAvail;	
	
	/**
	 * runnable to launch threads to manage particle manipulations
	 */
	public MPM_CPUPartBuilder partThdMgr;
	/**
	 * runnable to build and manage grid manipulations
	 */
	public MPM_CPUGridBuilder gridThdMgr;
	
	/**
	 * all particles representing the material in this simulation
	 */
	public MPM_CPURndrdPart[] parts;
	/**
	 * Radius for drawn MPM points
	 */
	protected float partRad;
	/**
	 * all grid nodes
	 */
	public MPM_CPUGridNode[][][] grid;	
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
	@SuppressWarnings("unchecked")
	public Base_MPMCPUSim(IRenderInterface _pa, Base_MPMSimWindow _win, String _simName, MPM_SimUpdateFromUIData _currUIVals) {
		super(_pa, _win, _simName, new float[] {0, 0, -9.8f}, _currUIVals);
		
		th_exec = Executors.newCachedThreadPool();	
		numThreadsAvail = win.getNumThreadsAvailable() - 2;
		//initialize active nodes set - array of sets, array membership is node ID % numThreadsAvail
		//might not be balanced, but worst case all nodes in single thread.  not many nodes, and minimal calculation is node-specific ConcurrentHashMap<myGridNode,Boolean>[]
		activeNodes = new ConcurrentHashMap[numThreadsAvail];
		// a map of active nodes, keyed by grid node, for each thread
		for(int i=0;i<this.numThreadsAvail;++i) {activeNodes[i] = new ConcurrentHashMap<MPM_CPUGridNode, MPM_CPUActiveNodeAgg>();}		
		//set up threads to be used to build the grid
		gridThdMgr = new MPM_CPUGridBuilder(this, numThreadsAvail);
		//set up thread builders to be used to build particles
		partThdMgr = new MPM_CPUPartBuilder(this, numThreadsAvail);	
		
		//set up grid and initialize sim with UI values and reset sim
		updateSimVals_FromUI(_currUIVals);
	}
	
	/**
	 * Instance-specific reset code
	 * @param rebuildSim
	 */
	@Override
	protected final void resetSim_Indiv(SimResetProcess rebuildSim) {
		msgObj.dispDebugMessage("Base_MPMCPUSim("+simName+")", "resetSim_Indiv","Start resetting sim");
	
		initValues_Parts();		
		initValues_Grid();
		
		//notify gridbuilder that we have a new grid
		gridThdMgr.setNewSimVals();
		partThdMgr.setNewSimVals();
		msgObj.dispDebugMessage("Base_MPMCPUSim("+simName+")", "resetSim_Indiv","Finished resetting sim");
	}//resetSim_Indiv
	

	@Override
	protected final void initPartArrays() {
		//build array of particles. UI has been updated by here
		parts = new MPM_CPURndrdPart[numParts];		
		//radius of drawn particles
		partRad = 10.0f/this.minSclAmt;
	}

	@Override
	protected final void buildPartLayouts() {
		buildPartLayoutArray(parts);
	}

	@Override
	protected final void reinitSimObjects() {
		reinitSimObjects(parts);
	}
	
	/**
	 * build initial layout for particles for this simulation
	 * @param partVals [OUT] array of particles
	 */
	protected abstract void buildPartLayoutArray(MPM_CPURndrdPart[] parts);	
	
	/**
	 * Reinitialize existing sim - will resynthesize sample points but will not rederive locations of sim objs.
	 * @param partVals
	 */
	protected abstract void reinitSimObjects(MPM_CPURndrdPart[] parts);	
	
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
	
	/**
	 * Initialize values based on particle positions, if any
	 */
	@Override
	protected final void initValues_Parts() {
		
		
	}
	
	
	@Override
	protected final void initValues_Grid() {}
	
	
	@Override
	protected final boolean simMe_Indiv(float modAmtMillis) {		return false;	}
	@Override
	protected final boolean simMePost_Indiv(float modAmtMillis) {	return false;	}
	@Override
	protected final void _drawParts(float animTimeMod, boolean showLocColors) {
		
		
		
	}
	/**
	 * Draw instance class particle velocities
	 * @param animTimeMod
	 * @param minMag minimum magnitude per axis to draw vector
	 * @param pincr
	 */
	@Override
	protected final void _drawPartVel(float animTimeMod, float vecScale, int pincr) {
		
		
	}
	/**
	 * Draw any colliders if they exist
	 * @param animTimeMod
	 */
	@Override
	protected final void _drawColliders(float animTimeMod) {
		ri.pushMatState();	
		drawColliders_Indiv(animTimeMod);
		ri.popMatState();
	}
	/**
	 * draw internal-to-sim colliders, if they exist
	 * @param animTimeMod
	 */
	protected abstract void drawColliders_Indiv(float animTimeMod);

	/**
	 * Draw instance class grid velocities - use _drawGridVec method
	 * @param animTimeMod
	 * @param minMag minimum magnitude per axis to draw vector
	 */
	@Override
	protected final void _drawGridVel(float animTimeMod, float minMag) {}

	/**
	 * Draw instance class grid accelerations - use _drawGridVec method
	 * @param animTimeMod
	 * @param minMag minimum magnitude per axis to draw vector
	 */
	@Override
	protected final void _drawGridAccel(float animTimeMod, float minMag) {}

	/**
	 * Draw instance class grid masses - use _drawGridScalar method
	 * @param animTimeMod
	 * @param minMag minimum magnitude to draw scalar mass
	 */
	@Override
	protected final void _drawGridMass(float animTimeMod, float minMag) {}

	/**
	 * check sim-specific central floating collider - 
	 * return 0 vec if no collision, otherwise return normal of collision
	 * @param pos
	 * @return
	 */
	public abstract myVectorf checkColliderCollision(myPointf pos);
}//class Base_MPMCPUSim
