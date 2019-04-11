package MPM_CudaSim;

import java.util.HashSet;
import java.util.concurrent.*;

import org.jblas.FloatMatrix;

////1 thread's worth of execution for executing work
//public abstract class mySimThreadExec implements Callable<Boolean> {
//	//ref to owning sim
//	protected MPM_ABS_Sim sim;
//
//	//sim step now running
//	protected int simStep;
//	//idxs in component array that this thread will map to
//	protected int stIDX, endIDX, thIDX;
//	
//	public mySimThreadExec(MPM_ABS_Sim _sim, int _thIDX, int _stIDX, int _endIDX) {
//		sim=_sim;stIDX = _stIDX; endIDX = _endIDX;thIDX = _thIDX;		
//		simStep=0;
//	}
//	//set what step of sim is being executed
//	public void setSimStep(int _s) {simStep = _s;}	
//			
//	//methods to execute different simulations steps - not all thread exec processes will use all step functs
//	protected abstract void execSimStep0();
//	protected abstract void execSimStep1();
//	protected abstract void execSimStep2();
//	protected abstract void execSimStep3();
//	protected abstract void execSimStep4();
//	protected abstract void execSimStep5();
//	
//	@Override
//	public Boolean call() throws Exception { 
//		switch(simStep) {
//			case 0 : {	 execSimStep0();return true;		}
//			case 1 : {	 execSimStep1();return true;		}
//			case 2 : {	 execSimStep2();return true;		}
//			case 3 : {	 execSimStep3();return true;		}
//			case 4 : {	 execSimStep4();return true;		}
//			case 5 : {	 execSimStep5();return true;		}
//			default : { return true;}
//		}
//	}
//	
//	
//	//cube wall normals
//	protected final myVectorf[] wallNorms = new myVectorf[] {
//			new myVectorf( 1.0, 0.0, 0.0),new myVectorf(-1.0, 0.0, 0.0),
//			new myVectorf( 0.0, 1.0, 0.0),new myVectorf( 0.0,-1.0, 0.0),
//			new myVectorf( 0.0, 0.0, 1.0),new myVectorf( 0.0, 0.0,-1.0)
//	};
//	
//	//calc collision velocity if collision is going to occur - norm is collider normal
//	protected myVectorf calcCollVel(float fricCoeff, myVectorf velocity, myVectorf norm) {		
//		float velNormDir=myVectorf._dot(velocity, norm);
//		if(velNormDir >=0){return velocity;}
//		//velocity in opposite direction of normal if dot prod is <0	
//		myVectorf velTanDir = myVectorf._sub(velocity, myVectorf._mult(norm,velNormDir));
//		float fricNormVel = fricCoeff*velNormDir;
//		if(velTanDir.magn <= -fricNormVel) {//coulomb friction law
//			//no bounce, so just stops
//			return new myVectorf(0.0,0.0,0.0);
//		}		
//		return myVectorf._add(velTanDir,(myVectorf._normalize(velTanDir))._mult(fricNormVel));
//		 
//	}//calcCollVel
//	
//	//collision detection for wall collisions, returns idx in wallNorms of collision
//	protected int checkWallCollision(myVectorf pos) {
//		if(pos.x<=sim.minSimBnds) {		return 0;	} else if(pos.x>=sim.maxSimBnds) {	return 1;	}			
//		if(pos.y<=sim.minSimBnds) {		return 2;	} else if(pos.y>=sim.maxSimBnds) {	return 3;	}
//		if(pos.z<=sim.minSimBnds) {		return 4;	} else if(pos.z>=sim.maxSimBnds) {	return 5;	}
//		return -1;
//	}//checkWallCollision
//	
//	//looking ahead at position, to return new velocity that responds to collision if position has collided
//	protected myVectorf applyCollisions(myVectorf pos, myVectorf velocity) {
//		//check collider collisions
//		myVectorf sphereColNorm = sim.checkColliderCollision(pos);
//		if(sphereColNorm.magn != 0) {//colliding with central collider			
//			return calcCollVel(sim.collFric, velocity, sphereColNorm);
//			
//		} else {//not colliding with sphere collider, check walls	
//			//check wall collisions
//			int colType = checkWallCollision(pos);
//			if(-1==colType) {return velocity;}
//			//calc collision velocity if collision is going to occur
//			return calcCollVel(sim.wallFric, velocity, wallNorms[colType]);
//		}
//	}//applyCollisions
//
//	
//}//class simThreadExec
//
////1 thread's worth of execution for building particles
//class partThreadExec extends mySimThreadExec {
//	//ref to particle array to be filled
//	myRndrdPart[] parts;
//	//center of sphere being built
//	myVectorf ctr;
//	//radius of sphere to be built
//	float rad;
//	
//	public partThreadExec(MPM_ABS_Sim _sim, int _thIDX, int _stIDX, int _endIDX, float _rad, myVectorf _ctr, myRndrdPart[] _parts) {
//		super(_sim,_thIDX, _stIDX, _endIDX);
//		parts=_parts;	
//		ctr = new myVectorf(_ctr);
//		rad = _rad;
//	}	
//	
//	private static final double lcl_third = 1.0/3.0;
//	//return a random position within a sphere of radius rad centered at ctr
//	private myVectorf getRandPosInSphere(double rad, myVectorf ctr){
//		myVectorf pos = new myVectorf();
//		double u = ThreadLocalRandom.current().nextDouble(0,1), r = rad * Math.pow(u, lcl_third);
//		do{
//			pos.set(ThreadLocalRandom.current().nextDouble(-1,1), ThreadLocalRandom.current().nextDouble(-1,1),ThreadLocalRandom.current().nextDouble(-1,1));
//		} while (pos.sqMagn > 1.0f);
//		pos._mult(r);
//		pos._add(ctr);
//		return pos;
//	}	
//	//initialization - build particles
//	@Override
//	protected void execSimStep0() {
//		for(int i=stIDX; i<endIDX;++i) {
//			parts[i]= new myRndrdPart(getRandPosInSphere(rad, ctr));
//			parts[i].thdIDX = thIDX;
//		}
//	}
//	
//	
//	private static final float twoThirds = 2.0f/3.0f, fourThirds = 2.0f*twoThirds;
//
//	private float[] calcWtDerivWt(MPM_ABS_Sim sim, float xVal) {
//		float xxAbs = Math.abs(xVal), xxSq = xVal * xVal, xxCu = xxAbs*xxSq;
//		float xxSign=Math.signum(xVal);
//		if(xxAbs < 1) {
//			return new float[] {(xxCu/2.0f - xxSq + twoThirds), (1.5f * xxSq*xxSign - 2.0f * xVal) / sim.h};			
//		} else if(xxAbs <2) {
//			//should never be >=2 based on how idxs are constructed so this will always be only if xAbs >=1 and <2
//			return new float[] {(-xxCu/6.0f + xxSq - 2*xxAbs + fourThirds), (-0.5f * xxSq *xxSign + 2.0f * xVal- 2.0f * xxSign)/ sim.h};
//		} else {
//			return new float[] {0.0f, 0.0f};
//		}
//	}//calcWtDerivWt	
//	
//	
//	private void findNeighborNodes(myRndrdPart p) {
//		p.neighbors.clear();
//		//precalculate to prevent repeated calcs in loops
//		float[] posRelPreCalc = {(p.pos.x - sim.minSimBnds)/ sim.h, (p.pos.y - sim.minSimBnds )/ sim.h,(p.pos.z - sim.minSimBnds)/ sim.h};
//		//no need to calculate end index - will always be 3 higher than start index
//		int stIdxI  = (int) (posRelPreCalc[0]) - 1,endIdxI = stIdxI + 3,
//			stIdxJ  = (int) (posRelPreCalc[1]) - 1,endIdxJ = stIdxJ + 3,
//			stIdxK  = (int) (posRelPreCalc[2]) - 1,endIdxK = stIdxK + 3;
//		//System.out.println("st idx : " + stIdxI + " | end idx : " + endIdxI + " Pos x (- min) : " + (pos.x - sim.minSimBnds) + " | posPreCalc : ["+posRelPreCalc[0]+","+posRelPreCalc[1]+","+posRelPreCalc[2]+"]");
//		for (int i = stIdxI; i <= endIdxI; ++i) {
//			if (0 <= i && i<sim.gridCount) {
//				//calc weight deriv has weight in first idx, and derive w/respect to the current index in the 2nd pos.
//				float[] x_wtVals = calcWtDerivWt(sim,posRelPreCalc[0] - i);
//				for (int j = stIdxJ; j <= endIdxJ; ++j) {
//					if (0 <= j && j<sim.gridCount) {
//						float[] y_wtVals = calcWtDerivWt(sim,posRelPreCalc[1] - j);
//						//precalcs
//						float xyWt00 = x_wtVals[0] * y_wtVals[0], xyWt01 = x_wtVals[0] * y_wtVals[1], xyWt10 = x_wtVals[1] * y_wtVals[0];
//						for (int k = stIdxK; k <= endIdxK; ++k) {
//							if (0 <= k && k<sim.gridCount) {
//								myGridNode node = sim.grid[i][j][k];
//								float[] z_wtVals = calcWtDerivWt(sim,posRelPreCalc[2] - k);	
//								//weights and differential weights/wrespect to each dir calculations
//								float weight = xyWt00 * z_wtVals[0]; 
//								float dweightx = xyWt10 * z_wtVals[0]; 
//								float dweighty = xyWt01 * z_wtVals[0];
//								float dweightz = xyWt00 * z_wtVals[1];								
//								activeNodeAgg agg = sim.addNodeToSet(node);	
//								nghbrNodeInfo tmp = new nghbrNodeInfo(p, node, weight, dweightx, dweighty, dweightz, agg);
//								p.neighbors.put(tmp, agg);
//							}//if k in bound
//						}//for k
//					}//if j in bound
//				}//for j
//			}//if i in bound
//		}//for i		
//	}//findNeighborNodes
//	
//	
//	//find neighborhood and add mass/vel to grid nodes
//	@Override
//	protected void execSimStep1() {
//		for(int i=stIDX; i<endIDX;++i) {findNeighborNodes(parts[i]);}//parts[i].findNeighborNodes(sim, sim.grid);}			
//		//KEEP SEPARATED
//		for(int i=stIDX; i<endIDX;++i) {for (nghbrNodeInfo ndInfo : parts[i].neighbors.keySet()) {ndInfo.addMassVelToNode(parts[i].mass);}}
//		//NEED TO CALL AGGREGATOR ONCE ALL THREADS ARE DONE
//	}//execSimStep1
//	
//	private float calcDet(FloatMatrix m) {
//		return m.get(0, 0) * m.get(1, 1) * m.get(2, 2)
//			   + m.get(0, 1) * m.get(1, 2) * m.get(2, 0)
//			   + m.get(1, 0) * m.get(2, 1) * m.get(0, 2)
//				   
//			   - m.get(2, 0) * m.get(1, 1) * m.get(0, 2)
//			   - m.get(1, 0) * m.get(0, 1) * m.get(2, 2)
//			   - m.get(2, 1) * m.get(1, 2) * m.get(0, 0);
//	}
//	
//	
//	//compute grid forces
//	@Override
//	protected void execSimStep2() {
//		float lambda0 = sim.mat.getLambda0(), mu0 = sim.mat.getMu0();
//		for(int i=stIDX; i<endIDX;++i) {
//			//Polar decomposition
//			FloatMatrix[] svd = org.jblas.Singular.fullSVD(parts[i].elasticDeformationGrad);
//			FloatMatrix Re = svd[0].mmul(svd[2].transpose());
//
//			//Compute the determinent of Fe and Fp
//			float Je = calcDet(parts[i].elasticDeformationGrad),  Jp = calcDet(parts[i].plasticDeformationGrad);	
//			
//			float hc1mJp = (float) Math.exp(sim.mat.hardeningCoeff * (1-Jp));
//			float scal1 = 2.0f * mu0 * (hc1mJp);
//			float scal2 = lambda0 * (hc1mJp) * (Je - 1);
//			FloatMatrix mat1 = (parts[i].elasticDeformationGrad.sub(Re)).mmul(parts[i].elasticDeformationGrad.transpose());
//			FloatMatrix mat2 = FloatMatrix.eye(3).mul((Je * scal2));
//
//			//FloatMatrix cauchyStressWoJpn = (mat1.mul(scal1)).add(mat2.mul(scal2));
//			FloatMatrix cauchyStressWoJpnMVol = ((mat1.mul(scal1)).add(mat2)).mul(parts[i].vol);
//			for (nghbrNodeInfo ndInfo : parts[i].neighbors.keySet()) {
//				ndInfo.addPartForce(cauchyStressWoJpnMVol);
//			}			//requires force aggregation in gridnode call
//			//NEED TO AGGREGATE ALL FORCES FOR ALL ACTIVE NODES AND SUBTRACT FROM EXISTING FORCES ONCE ALL THREADS ARE DONE
//		}		
//	}//execSimStep2
//	
//	//update deformation gradient
//	@Override
//	protected void execSimStep3() {
//		float delT = sim.getDeltaT();
//		//update the deformation gradient for each particle
//		for(int i=stIDX; i<endIDX;++i) {		parts[i].updDeformationGradient(sim.mat, delT);		}
//		for(int i=stIDX; i<endIDX;++i) {		parts[i].calcPartVel(sim.mat);	}
//		
//		for(int i=stIDX; i<endIDX;++i) {		
//			myVectorf pos=new myVectorf(parts[i].pos);
//			pos._add(myVectorf._mult(parts[i].vel,delT));
//			parts[i].vel=applyCollisions(pos,parts[i].vel);
//		}
//		
//		for(int i=stIDX; i<endIDX;++i) {		parts[i].advanceParticle(delT);	}		
//	}//execSimStep3
//	
//	@Override
//	protected void execSimStep4() {
//	}
//	
//	@Override
//	protected void execSimStep5() {}
//
//}//class simThreadExec
//
//class gridBuildWrkr extends mySimThreadExec {
//	myGridNode[][][] grid;
//	float minSimBnds, maxSimBnds, cellSize;
//	int gridCount;
//	ConcurrentHashMap<myGridNode, activeNodeAgg> nodeSet;
//	FloatMatrix gravBase;
//	
//	//will build the grid based on the passed span of idxs
//	public gridBuildWrkr(MPM_ABS_Sim _sim, int _thIDX, int _stIDX, int _endIDX, myGridNode[][][] _grid, int _gridSize, float _minBnd, float _maxBnd, float _cellSize) {
//		super(_sim,_thIDX, _stIDX, _endIDX);
//		grid=_grid;
//		gridCount = _gridSize;
//		minSimBnds = _minBnd;
//		maxSimBnds = _maxBnd;
//		cellSize = _cellSize;
//		//System.out.println("thdIdx : "+ thIDX + " | size : " + sim.activeNodes.length);
//		nodeSet = sim.activeNodes[thIDX];
//		gravBase = sim.gravity.dup();
//	}
//	
//	//initialize new grid
//	@Override
//	protected void execSimStep0() {
//		myGridNode[][] g1tmp;
//		myGridNode[] g2tmp;
//		//do a slice of the grid for each thread - slice along i axis of grid
//		for(int i=stIDX; i<endIDX; ++i) {
//			float gridPosX = i*cellSize + minSimBnds;
//			g1tmp = new myGridNode[gridCount][];
//			for(int j=0;j<gridCount;++j) {
//				float gridPosY = j*cellSize + minSimBnds;
//				g2tmp = new myGridNode[gridCount];
//				for(int k=0;k<gridCount;++k) {
//					g2tmp[k] = new myGridNode(new int[] {i,j,k}, new myVectorf(gridPosX, gridPosY, k*cellSize + minSimBnds));
//				}
//				g1tmp[j]=g2tmp;
//			}
//			grid[i]=g1tmp;
//		}
//	}//execSimStep0
//
//	//reset grid nodes
//	@Override
//	protected void execSimStep1() {
//		for (myGridNode n : nodeSet.keySet()) {n.reset();}
//	}
//	//aggregate all particle forces for every active grid node
//	@Override
//	protected void execSimStep2() {
//		for(activeNodeAgg agg : sim.activeNodes[thIDX].values()) {agg.aggregateMassAndVel();}
//	}
//	//aggregate all particle forces for every active grid node
//	@Override
//	protected void execSimStep3() {
//		float delT = sim.getDeltaT();
//		for(activeNodeAgg agg : sim.activeNodes[thIDX].values()) {agg.aggregateForces();}
//		for(myGridNode n : sim.activeNodes[thIDX].keySet()) {
//			if (n.hasMass()) {			n.calcVel(delT, gravBase.dup());}
//		}		
//		for(myGridNode n : sim.activeNodes[thIDX].keySet()) {
//			myVectorf newPos = myVectorf._add(n.pos, myVectorf._mult(n.newvelocity,delT));
//			n.newvelocity=applyCollisions(newPos,n.newvelocity);
//		}	
//	}
//	
//	//compute grid velocities and calc collisions
//	@Override
//	protected void execSimStep4() {
//	
//	}
//	//calc collisions
//	@Override
//	protected void execSimStep5() {
//
//	}	
//		
//}//class gridBuildWrkr
