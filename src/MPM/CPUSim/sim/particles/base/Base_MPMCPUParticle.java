package MPM.CPUSim.sim.particles.base;

/**
 * this file holds all physically-based sim components : classes to implement forces, particles, colliders, springs.
 */

import java.util.concurrent.ConcurrentHashMap;

import org.jblas.FloatMatrix;

import MPM.BaseSim.material.MPM_Material;
import MPM.CPUSim.sim.base.Base_MPMCPUSim;
import MPM.CPUSim.sim.grid.MPM_CPUActiveNodeAgg;
import MPM.CPUSim.sim.grid.MPM_CPUGridNode;
import MPM.CPUSim.sim.particles.MPM_CPUNeighborNodeInfo;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;


public class Base_MPMCPUParticle {
	public final int ID;
	public static int IDgen = 0;
	
	//thread index owning this particle
	public int thdIDX;
	
	public myPointf pos;
	public myVectorf vel, frc;
	//3 incex values of next lowest node IDX == discretized version of location
	public int[] nodeIDX;	

	//particle mass and volume
	public static float mass = .010f;
	public float vol; 
	
	//for outputting results
	protected static final String strfmt = "%.5f";

	
	public ConcurrentHashMap<MPM_CPUNeighborNodeInfo,MPM_CPUActiveNodeAgg> neighbors;
	
	public FloatMatrix plasticDeformationGrad = new FloatMatrix(new float[][] {{1, 0, 0},
		   																	   {0, 1, 0},
		   																	   {0, 0, 1}});
	public FloatMatrix elasticDeformationGrad = new FloatMatrix(new float[][] {{1, 0, 0},
			   																   {0, 1, 0},
			   																   {0, 0, 1}});	
	
	//deformation gradient calc structures - pre-init to speed up calculation
	private FloatMatrix clampedDiagM = new FloatMatrix(new float[][] {{1, 0, 0},
 		  																{0, 1, 0},
 		  																{0, 0, 1}});
	private FloatMatrix invClampedDiagM = new FloatMatrix(new float[][] {{1, 0, 0},
																		 {0, 1, 0},
																		 {0, 0, 1}});
	
	//values used for grid force and deformation gradient update calculations
	public FloatMatrix cauchyStressWoJpnMVol;
	private FloatMatrix dvp;
	
		
	public Base_MPMCPUParticle(myPointf _iPos) {
		ID = IDgen++;
		pos = new myPointf(_iPos);
		vel = new myVectorf(0.0,0.0,0.0);
		frc = new myVectorf();
		neighbors = new ConcurrentHashMap<MPM_CPUNeighborNodeInfo,MPM_CPUActiveNodeAgg>();
	}
	
	public void reset() {
		
	}
	
	public void applyForce(myVectorf _force) {frc._add(_force);}//applyforce
	
	//for each particle, update using flip and pic velocities
	//alpha commonly set to .95
	public void updateVelocity(float alpha, myVectorf velPic, myVectorf velFlip) {
		//perform like this so that intermediate magnitudes are not calculated over and over (as they would be using the myVectorf funcs)
		float x = velPic.x + alpha*(velFlip.x - velPic.x);
		float y = velPic.y + alpha*(velFlip.y - velPic.y);
		float z = velPic.z + alpha*(velFlip.z - velPic.z);
		vel.set(x,y,z);	
	}
	
	//move particle (fwd euler update)
	public void advanceParticle(float deltaT){
		//perform like this so that intermediate magnitudes are not calculated over and over (as they would be using the myVectorf funcs)		
		float x = pos.x + deltaT*vel.x;
		float y = pos.y + deltaT*vel.y;
		float z = pos.z + deltaT*vel.z;		
		pos.set(x,y,z);
	}
	
	private synchronized float[] calcWtDerivWt(float h, float xVal) {
		float xxAbs = Math.abs(xVal), xxSq = xVal * xVal, xxCu = xxAbs*xxSq;
		float xxSign=Math.signum(xVal);
		if(xxAbs < 1) {
			return new float[] {(xxCu/2.0f - xxSq + twoThirds), (1.5f * xxSq*xxSign - 2.0f * xVal) / h};			
		} else if(xxAbs <2) {
			//should never be >=2 based on how idxs are constructed so this will always be only if xAbs >=1 and <2
			return new float[] {(-xxCu/6.0f + xxSq - 2*xxAbs + fourThirds), (-0.5f * xxSq *xxSign + 2.0f * xVal- 2.0f * xxSign)/ h};
		} else {
			return new float[] {0.0f, 0.0f};
		}
	}//calcWtDerivWt
	
	private static final float twoThirds = 2.0f/3.0f, fourThirds = 2.0f*twoThirds;
	//find close gridnodes
	public void findNeighborNodes(Base_MPMCPUSim sim, MPM_CPUGridNode[][][] grid){		
		neighbors.clear();
		float cellSize = sim.getCellSize();
		int gridCount = sim.getGridSideCount();
		//precalculate to prevent repeated calcs in loops
		float[] posRelPreCalc = {(pos.x - sim.getMinSimBnds().x)/ cellSize, (pos.y - sim.getMinSimBnds().y)/ cellSize,(pos.z - sim.getMinSimBnds().z)/ cellSize};
		//no need to calculate end index - will always be 3 higher than start index
		int stIdxI  = (int) (posRelPreCalc[0]) - 1,endIdxI = stIdxI + 3,
			stIdxJ  = (int) (posRelPreCalc[1]) - 1,endIdxJ = stIdxJ + 3,
			stIdxK  = (int) (posRelPreCalc[2]) - 1,endIdxK = stIdxK + 3;
		//System.out.println("st idx : " + stIdxI + " | end idx : " + endIdxI + " Pos x (- min) : " + (pos.x - sim.minSimBnds) + " | posPreCalc : ["+posRelPreCalc[0]+","+posRelPreCalc[1]+","+posRelPreCalc[2]+"]");
		if(0 > stIdxI) {stIdxI = 0;}
		if(0 > stIdxJ) {stIdxJ = 0;}
		if(0 > stIdxK) {stIdxK = 0;}
		if(endIdxI >= gridCount) {endIdxI = gridCount -1;}
		if(endIdxJ >= gridCount) {endIdxJ = gridCount -1;}
		if(endIdxK >= gridCount) {endIdxK = gridCount -1;}
		//precalc weights
		float[][] x_wtVals = new float[(endIdxI - stIdxI + 1)][];
		float[][] y_wtVals = new float[(endIdxJ - stIdxJ + 1)][];
		float[][] z_wtVals = new float[(endIdxK - stIdxK + 1)][];
		int idxI = 0, idxJ = 0, idxK = 0;
		//calc weight deriv has weight in first idx, and derive w/respect to the current index in the 2nd pos.
		for (int i = stIdxI; i <= endIdxI; ++i) {x_wtVals[idxI++] = calcWtDerivWt(cellSize,posRelPreCalc[0] - i);}
		for (int j = stIdxJ; j <= endIdxJ; ++j) {y_wtVals[idxJ++] = calcWtDerivWt(cellSize,posRelPreCalc[1] - j);}
		for (int k = stIdxK; k <= endIdxK; ++k) {z_wtVals[idxK++] = calcWtDerivWt(cellSize,posRelPreCalc[2] - k);}
		
		idxI = 0;
		for (int i = stIdxI; i <= endIdxI; ++i) {
			idxJ = 0;
			for (int j = stIdxJ; j <= endIdxJ; ++j) {
				//precalcs
				float xyWt00 = x_wtVals[idxI][0] * y_wtVals[idxJ][0], xyWt01 = x_wtVals[idxI][0] * y_wtVals[idxJ][1], xyWt10 = x_wtVals[idxI][1] * y_wtVals[idxJ][0];
				idxK = 0;
				for (int k = stIdxK; k <= endIdxK; ++k) {
					MPM_CPUGridNode node = grid[i][j][k];
					//weights and differential weights w/respect to each dir calculations
					float weight = xyWt00 * z_wtVals[idxK][0]; 
					float dweightx = xyWt10 * z_wtVals[idxK][0]; 
					float dweighty = xyWt01 * z_wtVals[idxK][0];
					float dweightz = xyWt00 * z_wtVals[idxK][1];								
					MPM_CPUActiveNodeAgg agg = sim.addNodeToSet(node);	
					MPM_CPUNeighborNodeInfo tmp = new MPM_CPUNeighborNodeInfo(this, node, weight, dweightx, dweighty, dweightz, agg);
					neighbors.put(tmp, agg);
					++idxK;
				}//for k
				++idxJ;
			}//for j
			++idxI;
		}//for i
	}//findNeighborNodes
	
	
	public void addMassVelToGrid() {
		for (MPM_CPUNeighborNodeInfo ndInfo : neighbors.keySet()) {ndInfo.addMassVelToNode(mass);}	
	}
	
	public void computePartVolumeAndDensity(float h3) {
		//System.out.println("Compute particle volumes and densities during the first iteration");
		float partDensity = 0;
		for (MPM_CPUNeighborNodeInfo ndInfo : neighbors.keySet()) {
			partDensity += (ndInfo.getWeightedNodeMass()) / h3;
		}
		vol = mass / partDensity;
	}//computePartVolumeAndDensity	
	
	private float calcDet(FloatMatrix m) {
		return m.get(0, 0) * m.get(1, 1) * m.get(2, 2)
			   + m.get(0, 1) * m.get(1, 2) * m.get(2, 0)
			   + m.get(1, 0) * m.get(2, 1) * m.get(0, 2)
				   
			   - m.get(2, 0) * m.get(1, 1) * m.get(0, 2)
			   - m.get(1, 0) * m.get(0, 1) * m.get(2, 2)
			   - m.get(2, 1) * m.get(1, 2) * m.get(0, 0);
	}
		
	
	public void compGridForces(float lambda0, float mu0, float hardeningCoeff) {
		//Compute grid forces		
		//Polar decomposition
		FloatMatrix[] svd = org.jblas.Singular.fullSVD(elasticDeformationGrad);
		FloatMatrix Re = svd[0].mmul(svd[2].transpose());

		//Compute the determinant of Fe and Fp
		float Je = calcDet(elasticDeformationGrad),  Jp = calcDet(plasticDeformationGrad);	
		
		float hc1mJp = (float) Math.exp(hardeningCoeff * (1-Jp));
		float scal1 = 2.0f * mu0 * (hc1mJp);
		float scal2 = lambda0 * (hc1mJp) * (Je - 1);
		FloatMatrix mat1 = (elasticDeformationGrad.sub(Re)).mmul(elasticDeformationGrad.transpose());
		FloatMatrix mat2 = FloatMatrix.eye(3).mul((Je * scal2));

		cauchyStressWoJpnMVol = ((mat1.mul(scal1)).add(mat2)).mul(vol);
		//distribute this quantity to all neighbor nodes
		for (MPM_CPUNeighborNodeInfo ndInfo : neighbors.keySet()) {
			ndInfo.addPartForce(cauchyStressWoJpnMVol);
		}
					
	}//compGridForces

	
//	//FloatMatrix toString results are in ROW MAJOR order - across each row first, then to next row - need transpose so values will be in correct order for matlab 
////	e.g. this function used on tmpDbgMat = new FloatMatrix(new float[][] {{1, 2, 3},
////		   															      {4, 5, 6},
////																	      {7, 8, 9}});]
////	gives output 1,4,7,2,5,8,3,6,9
//	private String getMatStr(FloatMatrix m, String fmt) {
//		//jblas FloatMatrix toString args : (String fmt, String open, String close, String colSep, String rowSep)
//		return m.transpose().toString(fmt, "","",",",",");
//	}
//	
//	//returns a string holding the format of the matrices that is used for the csv strings - all mats are assumed to be 9 element
//	public static String getMatHdrStr(String matName) {
//		String res = "";
//		for(int c=0;c<3;++c) {for(int r=0;r<3;++r) {res += matName+"_r"+r+"_c"+c+",";}}
//		return res;		
//	}
//	
//	//this method will retrieve the values used by the grid force calculation, and the result generated.
//	//they will be formatted in a comma-separated string to be saved to a file to be used of the matlab verification process
//	//call this right after compGridForces is calculated
//	public String getGridFrcCalcStrCSV(MPM_Material mat) {
//		return ""+mat.toGridFrcCalcStrCSV(strfmt)+","+getMatStr(elasticDeformationGrad, strfmt)+","+getMatStr(plasticDeformationGrad,strfmt) + getMatStr(cauchyStressWoJpnMVol,strfmt);	
//	}	
//	
//	//this method will retrieve the values used by the deformation gradient calculation, and the result generated.
//	//they will be formatted in a comma-separated string to be saved to a file to be used of the matlab verification process
//	//call this right after updDeformationGradient is calculated
//	public String getUpdDefGradStrCSV(MPM_Material mat, float deltaT) {
//		return ""+String.format(strfmt,  deltaT)+","+mat.toDefGradUpdCalcStrCSV(strfmt)+"," + getMatStr(dvp, strfmt)+ "," +getMatStr(elasticDeformationGrad, strfmt)+","+getMatStr(plasticDeformationGrad,strfmt);		
//	}
	
	/**
	 * update this particle's deformation structures based on neighboring node state	
	 * @param mat
	 * @param deltaT
	 */
	public void updDeformationGradient(MPM_Material mat, float deltaT) {
		dvp = FloatMatrix.zeros(3, 3);
		for (MPM_CPUNeighborNodeInfo ndInfo : neighbors.keySet()) {
			FloatMatrix toAdd = ndInfo.buildDVP();
			dvp.addi(toAdd);
		}			
		FloatMatrix tempNewFe = (FloatMatrix.eye(3).add(dvp.mul(deltaT))).mmul(elasticDeformationGrad);
		FloatMatrix[] svdDVP_ED = org.jblas.Singular.fullSVD(tempNewFe);

		float diag;
		for (int i=0; i<3; ++i) {
			diag = Math.min(Math.max(svdDVP_ED[1].get(i), 1-mat.getCriticalCompressionPtr()[0]), 1+mat.getCriticalStretchPtr()[0]);
			clampedDiagM.put(i, i,diag);
			invClampedDiagM.put(i,i,1.0f/diag);		
		}				
		elasticDeformationGrad = svdDVP_ED[0].mmul(clampedDiagM).mmul(svdDVP_ED[2].transpose());		
		plasticDeformationGrad = svdDVP_ED[2].mmul(invClampedDiagM).mmul(svdDVP_ED[0].transpose()).mmul((tempNewFe).mmul(plasticDeformationGrad));					
	}//updDeformationGradient	
	
	public myVectorf vPIC = new myVectorf();
	public void calcPartVel(MPM_Material mat) {	
		vPIC.clear();
		myVectorf vFLIP = new myVectorf(vel);
		for (MPM_CPUNeighborNodeInfo ndInfo : neighbors.keySet()) {		
			//the "calcPartVelPicFlip" function also updates this particle's vPIC
			vFLIP._add(ndInfo.calcPartVelPicFlip());
		}
		updateVelocity(mat.getAlphaPicFlipPtr()[00], vPIC, vFLIP);
	}//calcPartVel
		
	@Override
	public String toString(){
		String res = "ID : " + ID + "\tMass:"+mass+"\n";
		res +="\tPosition:"+pos.toStrBrf()+"\n";
		res +="\tVelocity:"+vel.toStrBrf()+"\n";
		res +="\tCurrentForces:"+frc.toStrBrf()+"\n";
		
		return res;		
	}
}//myParticle

