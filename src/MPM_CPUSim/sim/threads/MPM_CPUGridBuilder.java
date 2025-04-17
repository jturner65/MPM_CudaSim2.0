package MPM_CPUSim.sim.threads;

import java.util.concurrent.ConcurrentHashMap;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;
import MPM_CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM_CPUSim.sim.grid.MPM_CPUGridNode;
import MPM_CPUSim.sim.threads.base.Base_MPMCPUSimThdBldr;

public class MPM_CPUGridBuilder extends Base_MPMCPUSimThdBldr {
	//array of hashsets of active nodes - reference to construct in sim
	protected ConcurrentHashMap<MPM_CPUGridNode, MPM_CPUActiveNodeAgg>[] activeNodes;
	
	public MPM_CPUGridBuilder(Base_MPMCPUSim _sim, int _numThds) {super(_sim, _numThds);}
	
	@Override
	public void initSimBuilder() {
		//rebuild sim grid
		int gridCount = sim.getGridSideCount();
		sim.grid = new MPM_CPUGridNode[gridCount][][];
		//initializes all worker-related variables like idxs for each worker to start/end in within grid array
		initWrkrIDXs(gridCount, "GridBuilder");
	}//initSimBuilder
	
	//build arraylist of thread worker callables
	@Override
	protected void initWrkrs() {
		activeNodes = sim.activeNodes;
		//build thread list structure - this won't change unless simulation structure will change, so on reinitialization should reuse these
		int gridCount = sim.getGridSideCount();
		for(int i =0;i<numThds;++i) {
			//System.out.println("i:"+i+" | start idx : " + stIDXs[i] + " | end : " + endIDXs[i]);
			callExecs.add(new MPM_CPUGridBuildWrkr(sim, i,stIDXs[i],endIDXs[i], sim.grid, gridCount, sim.getMinSimBnds(), sim.getMaxSimBnds(), sim.getCellSize()));
		}		
	}//initWrkrs
	
	//after all threads are finished this method is called, with specific code to execute depending on which sim step was run
	@Override
	protected void endRun() {
		switch(simStep) {//execute multiple different steps within this environment
			case 0 : {//initialization of all grid cells	
				//let simulation know that grid is built
				//sim.setSimGridIsBuilt();
				//build particles in sim - need to call this builder to rebuild all particles once grid is built
				//sim.resetSim();		
				break;}			
			case 1: {//after resetting all nodes, clear out nodeset arrays
				//System.out.println("All threads start nodeset clear : activenodes.size = " + activeNodes.length);
				for (int i=0; i< activeNodes.length; ++i) {	activeNodes[i].clear();	}
				//run project code
		      	sim.partThdMgr.setSimStep(1);
		      	sim.partThdMgr.synchRun();	
				break;}
			case 2 : {//aggregating all added masses and velocities				
				//now call sim routine here to compute grid forces
				sim.partThdMgr.setSimStep(2);
				sim.partThdMgr.synchRun();
				break;}			
			case 3 : {//aggregating all added forces and calc grid velocity and check for collisions
				//now we need to update the deformation gradient
				sim.partThdMgr.setSimStep(3);
				sim.partThdMgr.synchRun();				
				break;}			
			case 4 : {
				
				break;}
			
			default : {}
		}
	}//endRun

	
}//class myGridBuilder



