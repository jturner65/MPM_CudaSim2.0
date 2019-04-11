package MPM_CudaSim;


import java.util.*;
import java.util.concurrent.*;

////class to hold all necessary shared functionality used by sim thread manager/executors
//public abstract class mySimThdBldr implements Runnable {
//	//ref to owning sim
//	protected MPM_ABS_Sim sim;
//	//# of threads to spawn off - do this based on threads supported on machine
//	protected int numThds;
//	//structs to hold lists of callables to execute
//	protected List<mySimThreadExec> callExecs;
//	protected List<Future<Boolean>> callExecsFtrs;	
//	//values changed since reinitialization
//	protected boolean newSimVals;	
//	//what step in sim we are currently in - 0 is init for all types of mySimThdBldr, 
//	protected int simStep;	
//	//array of thread start and end indexes for objects being distributed to all threads
//	protected int[] stIDXs, endIDXs;
//	
//	public mySimThdBldr(MPM_ABS_Sim _sim, int _numThds) {
//		sim=_sim;
//		numThds = _numThds;
//		//array idxs to be 
//		stIDXs = new int[numThds];
//		endIDXs = new int[numThds];
//		newSimVals = true;
//		callExecs = new ArrayList<mySimThreadExec>();
//		callExecsFtrs = new ArrayList<Future<Boolean>>();
//	}//ctor	
//	
//	//individual thread builder/exec's init code - only called from run.  externally, setNewSimVals needs to be set if wanting to re-init mySimThdBldr objects
//	protected abstract void initSimBuilder();
//	
//	//build arraylist of thread worker callables
//	protected abstract void initWrkrs();
//
//	//call from initSimBuilder whenever workers need to be rebuilt due to changes in simulation topology
//	protected void initWrkrIDXs(int numObjs, String ownr) {
//		callExecs.clear();
//		callExecsFtrs.clear();
//		int numPartsPerThd= (numObjs > numThds ? (1 + (numObjs - 1)/numThds) : 1);
//		endIDXs[0]=numPartsPerThd;
//		for(int i=1;i<numThds;++i) {
//			stIDXs[i]=endIDXs[i-1];
//			endIDXs[i]=stIDXs[i]+numPartsPerThd;
//		}
//		if(endIDXs[numThds-1] >= numObjs) {endIDXs[numThds-1]=numObjs;}
//		else {System.out.println("mySimBuilder::initSimBuilder : Error initializing end idxs : not all objects are getting mapped to a thread in "+ ownr +" mySimBuilder : endIDXs[numThds-1] == "+endIDXs[numThds-1]);}
//		//for(int i=0;i<stIDXs.length;++i) {System.out.println(ownr + " : IDX : " + i + " StIDX : " + stIDXs[i]+ " | EndIDX : " + endIDXs[i]);}
//		//build workers specific to the instancing mySimThdBldr obj
//		initWrkrs();
//		//next call is initialization step for this sim builder
//		setSimStep(0);
//	}//initSimBuilder
//		
//	//set what step in simulation we are executing, and propagate into workers - only call externally for steps > 0 (initialization sets step 0 code)
//	public void setSimStep(int _s) {
//		simStep=_s;
//		//set step for all workers
//		for(mySimThreadExec sthd : callExecs) {sthd.setSimStep(_s);	}		
//	}
//	
//	//call whenever simulation values that impact building the simulation change (# of particles, grid dimensions, etc) - this will force sim step to be 0
//	public void setNewSimVals() {newSimVals = true;	}
//	
//	//after all threads are finished this method is called, with specific code to execute depending on which sim step was run
//	protected abstract void endRun();	
//	
//	
//	//call instead so waits on finishing
//	public void synchRun() {
//		//ref to particle array - handles if array is remade or other values impacting simulation architecture have changed, requiring a reset of the thread data
//		if(newSimVals) {initSimBuilder();newSimVals = false;}
//		try {
//			callExecsFtrs = sim.th_exec.invokeAll(callExecs);
//			for(Future<Boolean> f: callExecsFtrs) { f.get(); }
//		} catch (Exception e) { 
//			e.printStackTrace(); 
//		}	
//		endRun();	
//		
//		
//	}
//	
//	//re-initialize if new values have been set; launch all thread workers, block on futures from callables
//	@Override
//	public void run() {
//		//ref to particle array - handles if array is remade or other values impacting simulation architecture have changed, requiring a reset of the thread data
//		if(newSimVals) {initSimBuilder();newSimVals = false;}
//		try {
//			callExecsFtrs = sim.th_exec.invokeAll(callExecs);
//			for(Future<Boolean> f: callExecsFtrs) { f.get(); }
//		} catch (Exception e) { 
//			e.printStackTrace(); 
//		}	
//		endRun();
//	}//run
//
//}//abstract class mySimBuilder
//
//
////class to build a sphere of particles of specified size, dispersed amongst specified # of threads
//class myPartBuilder extends mySimThdBldr {
//	
//	public myPartBuilder(MPM_ABS_Sim _sim, int _numThds) {super(_sim, _numThds);}
//	
//	@Override
//	public void initSimBuilder() {	
//		initWrkrIDXs(sim.parts.length, "PartBuilder");
//	}//initSimBuilder	
//	
//	//build arraylist of thread worker callables
//	@Override
//	protected void initWrkrs() {
//		//build thread list structure - this won't change unless simulation structure will change, so on reinitialization should reuse these
//		for(int i =0;i<numThds;++i) {
//			//System.out.println("i:"+i+" | start idx : " + stIDXs[i] + " | end : " + endIDXs[i]);
//			callExecs.add(new partThreadExec(sim, i,stIDXs[i],endIDXs[i], sim.snoBallRad, sim.snoBallCtr, sim.parts));
//		}	
//	}//initWrkrs
//
//	//after all threads are finished this method is called, with specific code to execute depending on which sim step was run
//	@Override
//	protected void endRun() {
////		switch(simStep) {//execute multiple different steps within this environment
////			case 0 : {//initialization of all particles		
////				sim.setSimPartsAreBuilt();
////				sim.resetSimEnd();
////				break;}
////			case 1 : {//project to grid				
////				//once done, need to aggregate mass and vels in every grid node 
////		      	sim.gridThdMgr.setSimStep(2);//aggregate mass and vels
////		      	sim.gridThdMgr.synchRun();					
////				break;}
////			case 2 : {//compute grid forces from particles
////				//once done, need to aggregate forces in every grid node 
////				sim.gridThdMgr.setSimStep(3);//aggregate forces sim step
////				sim.gridThdMgr.synchRun();				      	
////				break;}
////			case 3 :{
////				break;}
////			case 4 : {//NOP				
////				break;}
////			default : {}
////		}
//	}//endRun
//
//}//class mySimBuilder
//
//
//class myGridBuilder extends mySimThdBldr {
//	//array of hashsets of active nodes - reference to construct in sim
//	protected ConcurrentHashMap<myGridNode, activeNodeAgg>[] activeNodes;
//	
//	public myGridBuilder(MPM_ABS_Sim _sim, int _numThds) {super(_sim, _numThds);}
//	
//	@Override
//	public void initSimBuilder() {
//		//rebuild sim grid
//		sim.grid = new myGridNode[sim.gridCount][][];
//		//initializes all worker-related variables like idxs for each worker to start/end in within grid array
//		initWrkrIDXs(sim.gridCount, "GridBuilder");
//	}//initSimBuilder
//	
//	//build arraylist of thread worker callables
//	@Override
//	protected void initWrkrs() {
//		activeNodes = sim.activeNodes;
//		//build thread list structure - this won't change unless simulation structure will change, so on reinitialization should reuse these
//		for(int i =0;i<numThds;++i) {
//			//System.out.println("i:"+i+" | start idx : " + stIDXs[i] + " | end : " + endIDXs[i]);
//			callExecs.add(new gridBuildWrkr(sim, i,stIDXs[i],endIDXs[i], sim.grid, sim.gridCount, sim.minSimBnds, sim.maxSimBnds, sim.h));
//		}		
//	}//initWrkrs
//	
//	//after all threads are finished this method is called, with specific code to execute depending on which sim step was run
//	@Override
//	protected void endRun() {
////		switch(simStep) {//execute multiple different steps within this environment
////			case 0 : {//initialization of all grid cells	
////				//let simulation know that grid is built
////				sim.setSimGridIsBuilt();
////				//build particles in sim - need to call this builder to rebuild all particles once grid is built
////				sim.resetSim();		
////				break;}			
////			case 1: {//after resetting all nodes, clear out nodeset arrays
////				//System.out.println("All threads start nodeset clear : activenodes.size = " + activeNodes.length);
////				for (int i=0; i< activeNodes.length; ++i) {	activeNodes[i].clear();	}
////				//run project code
////		      	sim.partThdMgr.setSimStep(1);
////		      	sim.partThdMgr.synchRun();	
////				break;}
////			case 2 : {//aggregating all added masses and velocities				
////				//now call sim routine here to compute grid forces
////				sim.partThdMgr.setSimStep(2);
////				sim.partThdMgr.synchRun();
////				break;}			
////			case 3 : {//aggregating all added forces and calc grid velocity and check for collisions
////				//now we need to update the deformation gradient
////				sim.partThdMgr.setSimStep(3);
////				sim.partThdMgr.synchRun();				
////				break;}			
////			case 4 : {
////				
////				break;}
////			
////			default : {}
////		}
//	}//endRun
//
//	
//}//class myGridBuilder
//
//
//
