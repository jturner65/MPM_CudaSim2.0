package MPM_CPUSim.sim.threads;

import java.util.concurrent.ConcurrentHashMap;

import org.jblas.FloatMatrix;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;
import MPM_CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM_CPUSim.sim.grid.MPM_CPUGridNode;
import MPM_CPUSim.sim.threads.base.Base_MPMCPUSimThreadExec;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;

public class MPM_CPUGridBuildWrkr extends Base_MPMCPUSimThreadExec {
	MPM_CPUGridNode[][][] grid;
	float minSimBnds, maxSimBnds, cellSize; 
	int gridCount;
	ConcurrentHashMap<MPM_CPUGridNode, MPM_CPUActiveNodeAgg> nodeSet;
	FloatMatrix gravBase;
	
	//will build the grid based on the passed span of idxs
	public MPM_CPUGridBuildWrkr(Base_MPMCPUSim _sim, int _thIDX, int _stIDX, int _endIDX, MPM_CPUGridNode[][][] _grid, int _gridSize, float _minBnd, float _maxBnd, float _cellSize) {
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
		MPM_CPUGridNode[][] g1tmp;
		MPM_CPUGridNode[] g2tmp;
		//do a slice of the grid for each thread - slice along i axis of grid
		for(int i=stIDX; i<endIDX; ++i) {
			float gridPosX = i*cellSize + minSimBnds;
			g1tmp = new MPM_CPUGridNode[gridCount][];
			for(int j=0;j<gridCount;++j) {
				float gridPosY = j*cellSize + minSimBnds;
				g2tmp = new MPM_CPUGridNode[gridCount];
				for(int k=0;k<gridCount;++k) {
					g2tmp[k] = new MPM_CPUGridNode(new int[] {i,j,k}, new myPointf(gridPosX, gridPosY, k*cellSize + minSimBnds));
				}
				g1tmp[j]=g2tmp;
			}
			grid[i]=g1tmp;
		}
	}//execSimStep0

	//reset grid nodes
	@Override
	protected void execSimStep1() {
		for (MPM_CPUGridNode n : nodeSet.keySet()) {n.reset();}
	}
	//aggregate all particle forces for every active grid node
	@Override
	protected void execSimStep2() {
		for(MPM_CPUActiveNodeAgg agg : nodeSet.values()) {agg.aggregateMassAndVel();}
	}
	//aggregate all particle forces for every active grid node
	@Override
	protected void execSimStep3() {
		float deltaT = sim.getDeltaT();
		for(MPM_CPUActiveNodeAgg agg : nodeSet.values()) {agg.aggregateForces();}
		for(MPM_CPUGridNode n : nodeSet.keySet()) {
			if (n.hasMass()) {			n.calcVel(deltaT, gravBase.dup());}
		}		
		for(MPM_CPUGridNode n : nodeSet.keySet()) {
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
