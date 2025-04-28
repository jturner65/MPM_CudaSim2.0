package MPM.CPUSim.sim.threads;

import java.util.concurrent.ThreadLocalRandom;

import org.jblas.FloatMatrix;

import MPM.BaseSim.material.Base_MPMMaterial;
import MPM.CPUSim.sim.base.Base_MPMCPUSim;
import MPM.CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM.CPUSim.sim.grid.MPM_CPUGridNode;
import MPM.CPUSim.sim.particles.MPM_CPUNeighborNodeInfo;
import MPM.CPUSim.sim.particles.MPM_CPURndrdPart;
import MPM.CPUSim.sim.particles.base.Base_MPMCPUParticle;
import MPM.CPUSim.sim.threads.base.Base_MPMCPUSimThreadExec;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;


//1 thread's worth of execution for building particles
public class MPM_CPUPartBuildWrkr extends Base_MPMCPUSimThreadExec {
	/**
	 * ref to particle array to be filled
	 */
	protected MPM_CPURndrdPart[] parts;
	/**
	 * center of snowball sphere being built
	 */
	protected myVectorf snowBallCtr;
	/**
	 * radius of sphere of particles to be built
	 */
	protected float snowBallRad;
	/**
	 * Material quantities of particles
	 */
	protected Base_MPMMaterial mat;
		
	public MPM_CPUPartBuildWrkr(Base_MPMCPUSim _sim, int _stIDX, int _endIDX, int _thIDX) {
		super(_sim, _stIDX, _endIDX, _thIDX);
		parts=sim.parts;
		mat = sim.mat;
//		snowBallCtr = new myVectorf(_ctr);
//		snowBallRad = _rad;
	}	
	
	private static final double lcl_third = 1.0/3.0;
	//return a random position within a sphere of radius rad centered at ctr
	private myPointf getRandPosInSphere(double rad, myVectorf ctr){
		myVectorf pos = new myVectorf();
		double u = ThreadLocalRandom.current().nextDouble(0,1), r = rad * Math.pow(u, lcl_third);
		do{
			pos.set(ThreadLocalRandom.current().nextDouble(-1,1), ThreadLocalRandom.current().nextDouble(-1,1),ThreadLocalRandom.current().nextDouble(-1,1));
		} while (pos.sqMagn > 1.0f);
		pos._mult(r);
		pos._add(ctr);
		return pos;
	}	
	/**
	 * initialization - build particles
	 */
	@Override
	protected void execSimStep0() {
		for(int i=stIDX; i<endIDX;++i) {
			parts[i]= new MPM_CPURndrdPart(getRandPosInSphere(snowBallRad, snowBallCtr));
			parts[i].thdIDX = thIDX;
		}
	}
	
	
	private static final float twoThirds = 2.0f/3.0f, fourThirds = 2.0f*twoThirds;

	private float[] calcWtDerivWt(Base_MPMCPUSim sim, float xVal) {
		float xxAbs = Math.abs(xVal), xxSq = xVal * xVal, xxCu = xxAbs*xxSq;
		float xxSign=Math.signum(xVal);
		float cellSize = sim.getCellSize();
		if(xxAbs < 1) {
			return new float[] {(xxCu/2.0f - xxSq + twoThirds), (1.5f * xxSq*xxSign - 2.0f * xVal) / cellSize};			
		} else if(xxAbs <2) {
			//should never be >=2 based on how idxs are constructed so this will always be only if xAbs >=1 and <2
			return new float[] {(-xxCu/6.0f + xxSq - 2*xxAbs + fourThirds), (-0.5f * xxSq *xxSign + 2.0f * xVal- 2.0f * xxSign)/ cellSize};
		} else {
			return new float[] {0.0f, 0.0f};
		}
	}//calcWtDerivWt	
	
	
	private void findNeighborNodes(MPM_CPURndrdPart p) {
		p.neighbors.clear();
		float cellSize = sim.getCellSize();
		int gridCount = sim.getGridSideCount();
		//precalculate to prevent repeated calcs in loops
		float[] posRelPreCalc = {(p.pos.x - minSimBnds.x)/ cellSize, (p.pos.y - minSimBnds.y)/ cellSize,(p.pos.z - minSimBnds.z)/ cellSize};
		//no need to calculate end index - will always be 3 higher than start index
		int stIdxI  = (int) (posRelPreCalc[0]) - 1,endIdxI = stIdxI + 3,
			stIdxJ  = (int) (posRelPreCalc[1]) - 1,endIdxJ = stIdxJ + 3,
			stIdxK  = (int) (posRelPreCalc[2]) - 1,endIdxK = stIdxK + 3;
		//System.out.println("st idx : " + stIdxI + " | end idx : " + endIdxI + " Pos x (- min) : " + (pos.x - sim.getMinSimBnds()) + " | posPreCalc : ["+posRelPreCalc[0]+","+posRelPreCalc[1]+","+posRelPreCalc[2]+"]");
		for (int i = stIdxI; i <= endIdxI; ++i) {
			if (0 <= i && i<gridCount) {
				//calc weight deriv has weight in first idx, and derive w/respect to the current index in the 2nd pos.
				float[] x_wtVals = calcWtDerivWt(sim,posRelPreCalc[0] - i);
				for (int j = stIdxJ; j <= endIdxJ; ++j) {
					if (0 <= j && j<gridCount) {
						float[] y_wtVals = calcWtDerivWt(sim,posRelPreCalc[1] - j);
						//precalcs
						float xyWt00 = x_wtVals[0] * y_wtVals[0], xyWt01 = x_wtVals[0] * y_wtVals[1], xyWt10 = x_wtVals[1] * y_wtVals[0];
						for (int k = stIdxK; k <= endIdxK; ++k) {
							if (0 <= k && k<gridCount) {
								MPM_CPUGridNode node = sim.grid[i][j][k];
								float[] z_wtVals = calcWtDerivWt(sim,posRelPreCalc[2] - k);	
								//weights and differential weights/wrespect to each dir calculations
								float weight = xyWt00 * z_wtVals[0]; 
								float dweightx = xyWt10 * z_wtVals[0]; 
								float dweighty = xyWt01 * z_wtVals[0];
								float dweightz = xyWt00 * z_wtVals[1];								
								MPM_CPUActiveNodeAgg agg = sim.addNodeToSet(node);	
								MPM_CPUNeighborNodeInfo tmp = new MPM_CPUNeighborNodeInfo(p, node, weight, dweightx, dweighty, dweightz, agg);
								p.neighbors.put(tmp, agg);
							}//if k in bound
						}//for k
					}//if j in bound
				}//for j
			}//if i in bound
		}//for i		
	}//findNeighborNodes
	
	
	//find neighborhood and add mass/vel to grid nodes
	@Override
	protected void execSimStep1() {
		for(int i=stIDX; i<endIDX;++i) {findNeighborNodes(parts[i]);}//parts[i].findNeighborNodes(sim, sim.grid);}			
		//KEEP SEPARATED
		for(int i=stIDX; i<endIDX;++i) {for (MPM_CPUNeighborNodeInfo ndInfo : parts[i].neighbors.keySet()) {ndInfo.addMassVelToNode(Base_MPMCPUParticle.mass);}}
		//NEED TO CALL AGGREGATOR ONCE ALL THREADS ARE DONE
	}//execSimStep1
	
	private float calcDet(FloatMatrix m) {
		return m.get(0, 0) * m.get(1, 1) * m.get(2, 2)
			   + m.get(0, 1) * m.get(1, 2) * m.get(2, 0)
			   + m.get(1, 0) * m.get(2, 1) * m.get(0, 2)
				   
			   - m.get(2, 0) * m.get(1, 1) * m.get(0, 2)
			   - m.get(1, 0) * m.get(0, 1) * m.get(2, 2)
			   - m.get(2, 1) * m.get(1, 2) * m.get(0, 0);
	}
	
	
	//compute grid forces
	@Override
	protected void execSimStep2() {
		float lambda0 = sim.mat.getLambda0Ptr()[0], mu0 = sim.mat.getMu0Ptr()[0];
		for(int i=stIDX; i<endIDX;++i) {
			//Polar decomposition
			FloatMatrix[] svd = org.jblas.Singular.fullSVD(parts[i].elasticDeformationGrad);
			FloatMatrix Re = svd[0].mmul(svd[2].transpose());

			//Compute the determinant of Fe and Fp
			float Je = calcDet(parts[i].elasticDeformationGrad),  Jp = calcDet(parts[i].plasticDeformationGrad);	
			
			float hc1mJp = (float) Math.exp(sim.mat.getHardeningCoeffPtr()[0] * (1-Jp));
			float scal1 = 2.0f * mu0 * (hc1mJp);
			float scal2 = lambda0 * (hc1mJp) * (Je - 1);
			var elastDefGrad = parts[i].elasticDeformationGrad.dup();
			FloatMatrix mat1 = (elastDefGrad.sub(Re)).mmul(parts[i].elasticDeformationGrad.transpose());
			FloatMatrix mat2 = FloatMatrix.eye(3).mul((Je * scal2));

			//FloatMatrix cauchyStressWoJpn = (mat1.mul(scal1)).add(mat2.mul(scal2));
			parts[i].cauchyStressWoJpnMVol = ((mat1.mul(scal1)).add(mat2)).mul(parts[i].vol);
			for (MPM_CPUNeighborNodeInfo ndInfo : parts[i].neighbors.keySet()) {
				ndInfo.addPartForce(parts[i].cauchyStressWoJpnMVol);
			}			//requires force aggregation in gridnode call
			//NEED TO AGGREGATE ALL FORCES FOR ALL ACTIVE NODES AND SUBTRACT FROM EXISTING FORCES ONCE ALL THREADS ARE DONE
		}		
	}//execSimStep2
	
	//update deformation gradient
	@Override
	protected void execSimStep3() {
		//update the deformation gradient for each particle
		float deltaT = sim.getDeltaT();
		for(int i=stIDX; i<endIDX;++i) {		parts[i].updDeformationGradient(sim.mat,deltaT);		}
		for(int i=stIDX; i<endIDX;++i) {		parts[i].calcPartVel(sim.mat);	}
		
		for(int i=stIDX; i<endIDX;++i) {		
			myPointf pos=new myPointf(parts[i].pos);
			pos._add(myVectorf._mult(parts[i].vel,deltaT));
			parts[i].vel=applyCollisions(pos,parts[i].vel);
		}
		
		for(int i=stIDX; i<endIDX;++i) {		parts[i].advanceParticle(deltaT);	}		
	}//execSimStep3
	
	@Override
	protected void execSimStep4() {
	}
	
	@Override
	protected void execSimStep5() {}

}//class partThreadExec