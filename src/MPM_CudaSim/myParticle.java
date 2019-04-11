package MPM_CudaSim;

/**
 * this file holds all physically-based sim components : classes to implement forces, particles, colliders, springs.
 */

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import org.jblas.FloatMatrix;


//public class myParticle {
//	public final int ID;
//	public static int IDgen = 0;
//	
//	//thread index owning this particle
//	public int thdIDX;
//	
//	public myVectorf pos, vel, frc;
//	//3 incex values of next lowest node IDX == discretized version of location
//	public int[] nodeIDX;	
//
//	//particle mass and volume
//	public static float mass = .01f*.01f*.01f*100.0f;
//	public static float density = 100.0f;
//	public float vol; 
//	
//	ConcurrentHashMap<nghbrNodeInfo,activeNodeAgg> neighbors;
//	
//	public FloatMatrix plasticDeformationGrad = new FloatMatrix(new float[][] {{1, 0, 0},
//		   																	   {0, 1, 0},
//		   																	   {0, 0, 1}});
//	public FloatMatrix elasticDeformationGrad = new FloatMatrix(new float[][] {{1, 0, 0},
//			   																   {0, 1, 0},
//			   																   {0, 0, 1}});	
//	
//	//deformation gradient calc structures - pre-init to speed up calculation
//	private FloatMatrix clampedDiagM = new FloatMatrix(new float[][] {{1, 0, 0},
// 		  																{0, 1, 0},
// 		  																{0, 0, 1}});
//	private FloatMatrix invClampedDiagM = new FloatMatrix(new float[][] {{1, 0, 0},
//																		 {0, 1, 0},
//																		 {0, 0, 1}});
//	
//		
//	public myParticle(myVectorf _iPos) {
//		ID = IDgen++;
//		pos = new myVectorf(_iPos);
//		vel = new myVectorf(0.0,0.0,0.0);
//		frc = new myVectorf();
//		neighbors = new ConcurrentHashMap<nghbrNodeInfo,activeNodeAgg>();
//	}
//	
//	
//	public void applyForce(myVectorf _force) {frc._add(_force);}//applyforce
//	
//	//for each particle, update using flip and pic velocities
//	//alpha commonly set to .95
//	public void updateVelocity(float alpha, myVectorf velPic, myVectorf velFlip) {
//		//perform like this so that intermediate magnitudes are not calculated over and over (as they would be using the myVectorf funcs)
//		float x = velPic.x + alpha*(velFlip.x - velPic.x);
//		float y = velPic.y + alpha*(velFlip.y - velPic.y);
//		float z = velPic.z + alpha*(velFlip.z - velPic.z);
//		vel.set(x,y,z);	
//	}
//	
//	//move particle (fwd euler update)
//	public void advanceParticle(float deltaT){
//		//perform like this so that intermediate magnitudes are not calculated over and over (as they would be using the myVectorf funcs)		
//		float x = pos.x + deltaT*vel.x;
//		float y = pos.y + deltaT*vel.y;
//		float z = pos.z + deltaT*vel.z;		
//		pos.set(x,y,z);
//	}
//	
//	private synchronized float[] calcWtDerivWt(MPM_ABS_Sim sim, float xVal) {
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
//	private static final float twoThirds = 2.0f/3.0f, fourThirds = 2.0f*twoThirds;
//	//find close gridnodes
//	public void findNeighborNodes(MPM_ABS_Sim sim, myGridNode[][][] grid){		
//		neighbors.clear();
//		//precalculate to prevent repeated calcs in loops
//		float[] posRelPreCalc = {(pos.x - sim.minSimBnds)/ sim.h, (pos.y - sim.minSimBnds )/ sim.h,(pos.z - sim.minSimBnds)/ sim.h};
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
//								myGridNode node = grid[i][j][k];
//								float[] z_wtVals = calcWtDerivWt(sim,posRelPreCalc[2] - k);	
//								//weights and differential weights/wrespect to each dir calculations
//								float weight = xyWt00 * z_wtVals[0]; 
//								float dweightx = xyWt10 * z_wtVals[0]; 
//								float dweighty = xyWt01 * z_wtVals[0];
//								float dweightz = xyWt00 * z_wtVals[1];								
//								activeNodeAgg agg = sim.addNodeToSet(node);	
//								nghbrNodeInfo tmp = new nghbrNodeInfo(this, node, weight, dweightx, dweighty, dweightz, agg);
//								neighbors.put(tmp, agg);
//							}//if k in bound
//						}//for k
//					}//if j in bound
//				}//for j
//			}//if i in bound
//		}//for i
//	}//findNeighborNodes
//	
//	
//	public void addMassVelToGrid() {
//		for (nghbrNodeInfo ndInfo : neighbors.keySet()) {ndInfo.addMassVelToNode(mass);}	
//	}
//	
//	public void computePartVolumeAndDensity(float h3) {
//		//System.out.println("Compute particle volumes and densities during the first iteration");
//		float partDensity = 0;
//		for (nghbrNodeInfo ndInfo : neighbors.keySet()) {
//			partDensity += (ndInfo.getWeightedNodeMass()) / h3;
//		}
//		vol = mass / partDensity;
//	}//computePartVolumeAndDensity	
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
//	public void compGridForces(float lambda0, float mu0, float hardeningCoeff) {
//		//Compute grid forces
//
//		//FloatMatrix Fe = elasticDeformationGrad;
//		//FloatMatrix Fp = plasticDeformationGrad;
//		
//		//Polar decomposition
//		FloatMatrix[] svd = org.jblas.Singular.fullSVD(elasticDeformationGrad);
//		//FloatMatrix Re = svdED[0].mmul(svdED[2].transpose());
//		FloatMatrix Re = svd[0].mmul(svd[2].transpose());
//
//		//Compute the determinent of Fe and Fp
//		float Je = calcDet(elasticDeformationGrad),  Jp = calcDet(plasticDeformationGrad);	
//		
//		float hc1mJp = (float) Math.exp(hardeningCoeff * (1-Jp));
//		float scal1 = 2.0f * mu0 * (hc1mJp);
//		float scal2 = lambda0 * (hc1mJp) * (Je - 1);
//		FloatMatrix mat1 = (elasticDeformationGrad.sub(Re)).mmul(elasticDeformationGrad.transpose());
//		FloatMatrix mat2 = FloatMatrix.eye(3).mul((Je * scal2));
//
//		//FloatMatrix cauchyStressWoJpn = (mat1.mul(scal1)).add(mat2.mul(scal2));
//		FloatMatrix cauchyStressWoJpnMVol = ((mat1.mul(scal1)).add(mat2)).mul(vol);
//		for (nghbrNodeInfo ndInfo : neighbors.keySet()) {
//			ndInfo.addPartForce(cauchyStressWoJpnMVol);
//		}
//					
//	}//compGridForces
//	
//		
//	//update this particle's deformation structures based on neighboring node state	
//	public void updDeformationGradient(myMaterial mat, float deltaT) {
//		FloatMatrix dvp = FloatMatrix.zeros(3, 3);
//		for (nghbrNodeInfo ndInfo : neighbors.keySet()) {
//			FloatMatrix toAdd = ndInfo.buildDVP();
//			dvp.addi(toAdd);
//		}			
//		//if(ID==100) {System.out.println("");}
//		FloatMatrix tempNewFe = (FloatMatrix.eye(3).add(dvp.mul(deltaT))).mmul(elasticDeformationGrad);
//		//if(ID==100) {System.out.println(tempNewFe+"\n");}
//		FloatMatrix[] svdDVP_ED = org.jblas.Singular.fullSVD(tempNewFe);
//
//		float diag;
//		for (int i=0; i<3; ++i) {
//			diag = Math.min(Math.max(svdDVP_ED[1].get(i), 1-mat.criticalCompression), 1+mat.criticalStretch);
//			clampedDiagM.put(i, i,diag);
//			invClampedDiagM.put(i,i,1.0f/diag);		
//		}				
//		elasticDeformationGrad = svdDVP_ED[0].mmul(clampedDiagM).mmul(svdDVP_ED[2].transpose());		
//		plasticDeformationGrad = svdDVP_ED[2].mmul(invClampedDiagM).mmul(svdDVP_ED[0].transpose()).mmul((tempNewFe).mmul(plasticDeformationGrad));					
//	}//updDeformationGradient	
//	
//	public myVectorf vPIC = new myVectorf();
//	public void calcPartVel(myMaterial mat) {	
//		vPIC.clear();
//		myVectorf vFLIP = new myVectorf(vel);
//		for (nghbrNodeInfo ndInfo : neighbors.keySet()) {		
//			//the "calcPartVelPicFlip" function also updates this particle's vPIC
//			vFLIP._add(ndInfo.calcPartVelPicFlip());
//		}
//		updateVelocity(mat.alphaPicFlip, vPIC, vFLIP);
//	}//calcPartVel
//		
//	@Override
//	public String toString(){
//		String res = "ID : " + ID + "\tMass:"+mass+"\n";
//		res +="\tPosition:"+pos.toStrBrf()+"\n";
//		res +="\tVelocity:"+vel.toStrBrf()+"\n";
//		res +="\tCurrentForces:"+frc.toStrBrf()+"\n";
//		
//		return res;		
//	}
//}//myParticle
//
////individual particle of snow
//class myRndrdPart extends myParticle{
//	public int[] color;
//	//so we can reset color if we change it
//	public static final int[] origColor = new int[] {255,255,255,255};
//	
//	public myRndrdPart(myVectorf _iPos) {
//		super(_iPos);
//		color = new int[] {255,255,255,255};
//	}
//	
//	public void drawMe(MPM_SimMain pa, float partRad) {
//		pa.pushMatrix();pa.pushStyle();
//		pa.translate(pos);
//		pa.sphere(partRad);
//		pa.popStyle();pa.popMatrix();
//
//	}//drawMe
//	
//	public void drawMeVel(MPM_SimMain pa, float partRad) {
//		pa.pushMatrix();pa.pushStyle();
//		pa.translate(pos);
//		pa.sphere(partRad);
//		
//		pa.stroke(190,20,70,155);
//		pa.line(0,0,0,vel.x,vel.y,vel.z);
//		//if(ID==100) {pa.outStr2Scr(""+vel);}
//		pa.popStyle();pa.popMatrix();
//	}//drawMeDebug
//}//class myRndrdPart
//
////class to hold information regarding neighboring nodes to a particle
////put in single struct-like object so only have single list of nghbrNodeInfos
//class nghbrNodeInfo {
//	//particle this node is a neighbor to
//	private myParticle p;
//	//neighbor node. howdy neighbor!
//	private myGridNode node;
//	//this node's aggregator;
//	private activeNodeAgg agg;
//	private float weightToNode;
//	private FloatMatrix dWeightToNode;
//	
//	public nghbrNodeInfo(myParticle _p, myGridNode _node, float _weightToNode, float dweightx, float dweighty, float dweightz, activeNodeAgg _agg) {
//		p=_p;node=_node;agg=_agg;
//		weightToNode=_weightToNode; 		
//		dWeightToNode= new FloatMatrix(new float[] {dweightx, dweighty, dweightz});
//	}
//	//add owning particle's mass to this neighbor node
//	public void addMassVelToNode(float mass) {
//		float wtdMass = weightToNode * mass;
//		//momentum quantityt
//		myVectorf wtdVel = myVectorf._mult(p.vel,wtdMass);
//		agg.addPartMassAndWtdVel(wtdMass, wtdVel, p.thdIDX);
//	}
//	//get this neighbor node's mass weighted by particle's weighting on node
//	public float getWeightedNodeMass() {return node.getWtdMass(weightToNode);}
//	
//	//add forces from stress/strain calc to node
//	public void addPartForce(FloatMatrix cauchyStressWoJpnMVol) {
//		//add weighted particle force to aggregator
//		//FloatMatrix res = (cauchyStressWoJpnMVol).mmul(dWeightToNode);
//		agg.addPartFrc((cauchyStressWoJpnMVol).mmul(dWeightToNode), p.thdIDX);
//	}
//	
//	public FloatMatrix buildDVP() {
//		FloatMatrix vMat = new FloatMatrix(new float[] {node.newvelocity.x, node.newvelocity.y,node.newvelocity.z});
//		return vMat.mmul(dWeightToNode.transpose());
//	}
//	
//	public myVectorf calcPartVelPicFlip() {
//		p.vPIC._add(myVectorf._mult(node.newvelocity,weightToNode));
//		return node.getWtDiffVel(weightToNode);
//	}
//	
//	
//}//class nghbrNodeInfo
