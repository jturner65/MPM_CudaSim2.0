package MPM_CPUSim.sim;

import org.jblas.FloatMatrix;

import MPM_CPUSim.sim.grid.activeNodeAgg;
import MPM_CPUSim.sim.grid.myGridNode;
import MPM_CPUSim.sim.particles.base.myParticle;
import base_Math_Objects.vectorObjs.floats.myVectorf;

/**
 * class to hold information regarding neighboring nodes to a particle, 
 * put in single struct-like object so only have single list of nghbrNodeInfos
 */
public class nghbrNodeInfo {
	//particle this node is a neighbor to
	private myParticle p;
	//neighbor node. howdy neighbor!
	private myGridNode node;
	//this node's aggregator;
	private activeNodeAgg agg;
	private float weightToNode;
	private FloatMatrix dWeightToNode;
	
	public nghbrNodeInfo(myParticle _p, myGridNode _node, float _weightToNode, float dweightx, float dweighty, float dweightz, activeNodeAgg _agg) {
		p=_p;node=_node;agg=_agg;
		weightToNode=_weightToNode; 		
		dWeightToNode= new FloatMatrix(new float[] {dweightx, dweighty, dweightz});
	}
	//add owning particle's mass to this neighbor node
	public void addMassVelToNode(float mass) {
		float wtdMass = weightToNode * mass;
		//momentum quantityt
		myVectorf wtdVel = myVectorf._mult(p.vel,wtdMass);
		agg.addPartMassAndWtdVel(wtdMass, wtdVel, p.thdIDX);
	}
	//get this neighbor node's mass weighted by particle's weighting on node
	public float getWeightedNodeMass() {return node.getWtdMass(weightToNode);}
	
	//add forces from stress/strain calc to node
	public void addPartForce(FloatMatrix cauchyStressWoJpnMVol) {
		//add weighted particle force to aggregator
		//FloatMatrix res = (cauchyStressWoJpnMVol).mmul(dWeightToNode);
		agg.addPartFrc((cauchyStressWoJpnMVol).mmul(dWeightToNode), p.thdIDX);
	}
	
	public FloatMatrix buildDVP() {
		FloatMatrix vMat = new FloatMatrix(new float[] {node.newvelocity.x, node.newvelocity.y,node.newvelocity.z});
		return vMat.mmul(dWeightToNode.transpose());
	}
	
	//the "calcPartVelPicFlip" function also updates this particle's vPIC
	public myVectorf calcPartVelPicFlip() {
		p.vPIC._add(myVectorf._mult(node.newvelocity,weightToNode));
		return node.getWtDiffVel(weightToNode);
	}
	
	
}//class nghbrNodeInfo
