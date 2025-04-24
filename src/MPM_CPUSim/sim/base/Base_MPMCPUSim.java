package MPM_CPUSim.sim.base;


import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import MPM_BaseSim.sim.Base_MPMSim;
import MPM_BaseSim.sim.SimResetProcess;
import MPM_BaseSim.ui.Base_MPMSimWindow;
import MPM_BaseSim.utils.MPM_SimUpdateFromUIData;
import MPM_CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM_CPUSim.sim.grid.MPM_CPUGridNode;
import MPM_CPUSim.sim.particles.MPM_CPURndrdPart;
import MPM_CPUSim.sim.threads.MPM_CPUGridBuilder;
import MPM_CPUSim.sim.threads.MPM_CPUPartBuilder;
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
	 * all grid nodes
	 */
	public MPM_CPUGridNode[][][] grid;	
	/**
	 * all grid nodes needing updates - using concurrent hashmap
	 * because there's no concurrent hash set in java, but the keys of CHM will work
	 */
	public ConcurrentHashMap<MPM_CPUGridNode,MPM_CPUActiveNodeAgg>[] activeNodes;
	
	/**
	 * 
	 */
	
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
		for(int i=0;i<this.numThreadsAvail;++i) {activeNodes[i] = new ConcurrentHashMap<MPM_CPUGridNode, MPM_CPUActiveNodeAgg>();}		
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
		win.getMsgObj().dispDebugMessage("Base_MPMCPUSim("+simName+")", "resetSim_Indiv","Start resetting sim");
	
		initValues_Parts();
		
		
		initValues_Grids() ;
		win.getMsgObj().dispDebugMessage("Base_MPMCPUSim("+simName+")", "resetSim_Indiv","Finished resetting sim");
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
		
		
	}
	
	
	@Override
	protected final void initValues_Grids() {
		
		
	}
	
	
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
