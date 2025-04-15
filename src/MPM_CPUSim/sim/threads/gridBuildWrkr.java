package MPM_CPUSim.sim.threads;

import java.util.concurrent.ConcurrentHashMap;

import org.jblas.FloatMatrix;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;
import MPM_CPUSim.sim.grid.activeNodeAgg;
import MPM_CPUSim.sim.grid.myGridNode;
import MPM_CPUSim.sim.threads.base.mySimThreadExec;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;

public class gridBuildWrkr extends mySimThreadExec {
	myGridNode[][][] grid;
	float minSimBnds, maxSimBnds, cellSize; 
	int gridCount;
	ConcurrentHashMap<myGridNode, activeNodeAgg> nodeSet;
	FloatMatrix gravBase;
	
	//will build the grid based on the passed span of idxs
	public gridBuildWrkr(Base_MPMCPUSim _sim, int _thIDX, int _stIDX, int _endIDX, myGridNode[][][] _grid, int _gridSize, float _minBnd, float _maxBnd, float _cellSize) {
		super(_sim,_thIDX, _stIDX, _endIDX);
		grid=_grid;
		gridCount = _gridSize;
		minSimBnds = _minBnd;
		maxSimBnds = _maxBnd;
		cellSize = _cellSize;
		//System.out.println("thdIdx : "+ thIDX + " | size : " + sim.activeNodes.length);
		nodeSet = sim.activeNodes[thIDX];
		gravBase = sim.gravityMat.dup();
	}
	
	//initialize new grid
	@Override
	protected void execSimStep0() {
		myGridNode[][] g1tmp;
		myGridNode[] g2tmp;
		//do a slice of the grid for each thread - slice along i axis of grid
		for(int i=stIDX; i<endIDX; ++i) {
			float gridPosX = i*cellSize + minSimBnds;
			g1tmp = new myGridNode[gridCount][];
			for(int j=0;j<gridCount;++j) {
				float gridPosY = j*cellSize + minSimBnds;
				g2tmp = new myGridNode[gridCount];
				for(int k=0;k<gridCount;++k) {
					g2tmp[k] = new myGridNode(new int[] {i,j,k}, new myPointf(gridPosX, gridPosY, k*cellSize + minSimBnds));
				}
				g1tmp[j]=g2tmp;
			}
			grid[i]=g1tmp;
		}
	}//execSimStep0

	//reset grid nodes
	@Override
	protected void execSimStep1() {
		for (myGridNode n : nodeSet.keySet()) {n.reset();}
	}
	//aggregate all particle forces for every active grid node
	@Override
	protected void execSimStep2() {
		for(activeNodeAgg agg : nodeSet.values()) {agg.aggregateMassAndVel();}
	}
	//aggregate all particle forces for every active grid node
	@Override
	protected void execSimStep3() {
		float deltaT = sim.getDeltaT();
		for(activeNodeAgg agg : nodeSet.values()) {agg.aggregateForces();}
		for(myGridNode n : nodeSet.keySet()) {
			if (n.hasMass()) {			n.calcVel(deltaT, gravBase.dup());}
		}		
		for(myGridNode n : nodeSet.keySet()) {
			myVectorf newPos = myVectorf._add(n.pos, myVectorf._mult(n.newvelocity,deltaT));
			n.newvelocity=applyCollisions(newPos,n.newvelocity);
		}	
	}
	
	//compute grid velocities and calc collisions
	@Override
	protected void execSimStep4() {
	
	}
	//calc collisions
	@Override
	protected void execSimStep5() {

	}	

}//
