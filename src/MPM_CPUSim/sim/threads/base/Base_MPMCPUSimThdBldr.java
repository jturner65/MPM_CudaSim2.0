package MPM_CPUSim.sim.threads.base;


import java.util.*;
import java.util.concurrent.*;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;

//class to hold all necessary shared functionality used by sim thread manager/executors
public abstract class Base_MPMCPUSimThdBldr implements Runnable {
	//ref to owning sim
	protected Base_MPMCPUSim sim;
	//# of threads to spawn off - do this based on threads supported on machine
	protected int numThds;
	//structs to hold lists of callables to execute
	protected List<Base_MPMCPUSimThreadExec> callExecs;
	protected List<Future<Boolean>> callExecsFtrs;	
	//values changed since reinitialization
	protected boolean newSimVals;	
	//what step in sim we are currently in - 0 is init for all types of mySimThdBldr, 
	protected int simStep;	
	//array of thread start and end indexes for objects being distributed to all threads
	protected int[] stIDXs, endIDXs;
	
	public Base_MPMCPUSimThdBldr(Base_MPMCPUSim _sim, int _numThds) {
		sim=_sim;
		numThds = _numThds;
		//array idxs to be 
		stIDXs = new int[numThds];
		endIDXs = new int[numThds];
		newSimVals = true;
		callExecs = new ArrayList<Base_MPMCPUSimThreadExec>();
		callExecsFtrs = new ArrayList<Future<Boolean>>();
	}//ctor	
	
	//individual thread builder/exec's init code - only called from run.  externally, setNewSimVals needs to be set if wanting to re-init mySimThdBldr objects
	protected abstract void initSimBuilder();
	
	//build arraylist of thread worker callables
	protected abstract void initWrkrs();

	//call from initSimBuilder whenever workers need to be rebuilt due to changes in simulation topology
	protected void initWrkrIDXs(int numObjs, String ownr) {
		callExecs.clear();
		callExecsFtrs.clear();
		int numPartsPerThd= (numObjs > numThds ? (1 + (numObjs - 1)/numThds) : 1);
		endIDXs[0]=numPartsPerThd;
		for(int i=1;i<numThds;++i) {
			stIDXs[i]=endIDXs[i-1];
			endIDXs[i]=stIDXs[i]+numPartsPerThd;
		}
		if(endIDXs[numThds-1] >= numObjs) {endIDXs[numThds-1]=numObjs;}
		else {System.out.println("mySimBuilder::initSimBuilder : Error initializing end idxs : not all objects are getting mapped to a thread in "+ ownr +" mySimBuilder : endIDXs[numThds-1] == "+endIDXs[numThds-1]);}
		//for(int i=0;i<stIDXs.length;++i) {System.out.println(ownr + " : IDX : " + i + " StIDX : " + stIDXs[i]+ " | EndIDX : " + endIDXs[i]);}
		//build workers specific to the instancing mySimThdBldr obj
		initWrkrs();
		//next call is initialization step for this sim builder
		setSimStep(0);
	}//initSimBuilder
		
	//set what step in simulation we are executing, and propagate into workers - only call externally for steps > 0 (initialization sets step 0 code)
	public void setSimStep(int _s) {
		simStep=_s;
		//set step for all workers
		for(Base_MPMCPUSimThreadExec sthd : callExecs) {sthd.setSimStep(_s);	}		
	}
	
	//call whenever simulation values that impact building the simulation change (# of particles, grid dimensions, etc) - this will force sim step to be 0
	public void setNewSimVals() {newSimVals = true;	}
	
	//after all threads are finished this method is called, with specific code to execute depending on which sim step was run
	protected abstract void endRun();	
	
	
	//call instead so waits on finishing
	public void synchRun() {
		//ref to particle array - handles if array is remade or other values impacting simulation architecture have changed, requiring a reset of the thread data
		if(newSimVals) {initSimBuilder();newSimVals = false;}
		try {
			callExecsFtrs = sim.th_exec.invokeAll(callExecs);
			for(Future<Boolean> f: callExecsFtrs) { f.get(); }
		} catch (Exception e) { 
			e.printStackTrace(); 
		}	
		endRun();	
		
		
	}
	
	//re-initialize if new values have been set; launch all thread workers, block on futures from callables
	@Override
	public void run() {
		//ref to particle array - handles if array is remade or other values impacting simulation architecture have changed, requiring a reset of the thread data
		if(newSimVals) {initSimBuilder();newSimVals = false;}
		try {
			callExecsFtrs = sim.th_exec.invokeAll(callExecs);
			for(Future<Boolean> f: callExecsFtrs) { f.get(); }
		} catch (Exception e) { 
			e.printStackTrace(); 
		}	
		endRun();
	}//run

}//abstract class mySimBuilder



