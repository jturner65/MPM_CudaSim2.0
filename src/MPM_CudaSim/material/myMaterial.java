package MPM_CudaSim.material;

import MPM_CudaSim.utils.MPM_SimUpdateFromUIData;

//class for the material that the simulation consists of
public class myMaterial {
	private float initYoungMod, poissonRatio, lambda0, mu0;
	private float hardeningCoeff, criticalCompression, criticalStretch, alphaPicFlip;
		
	public myMaterial(MPM_SimUpdateFromUIData _currUIVals) {
		updateMatVals_FromUI(_currUIVals);
	}//ctor
		

	//recalc lambda0 and mu0 when young mod or poisson ratio changes
	private void recalcParams() {	
		lambda0 = (initYoungMod * poissonRatio) / ((1.0f+poissonRatio) * (1.0f-2.0f*poissonRatio));
		mu0 = initYoungMod / (2.0f * (1.0f+poissonRatio));
	}
	
	public void updateMatVals_FromUI(MPM_SimUpdateFromUIData upd) {		
		initYoungMod = upd.getInitYoungMod();
		poissonRatio = upd.getPoissonRatio();		
		hardeningCoeff = upd.getHardeningCoeff();     
		criticalCompression = upd.getCriticalCompression();
		criticalStretch = upd.getCriticalStretch();   
		alphaPicFlip = upd.getAlphaPicFlip();
		//recalc mu0 when ym or poisson ratio changes
		recalcParams();
	}//updateMatVals_FromUI
	
	
	//private to govern recalculation when youngmod and poisson ratio change via UI input
	
	// Return ptrs to values for kernel params
	public float[] getLambda0Ptr() {return new float[] {lambda0};}
	public float[] getMu0Ptr() {return new float[] {mu0};}	
	public float[] getHardeningCoeffPtr() {return new float[] {hardeningCoeff};}     
	public float[] getCriticalCompressionPtr() {return new float[] {criticalCompression};}  
	public float[] getCriticalStretchPtr() {return new float[] {criticalStretch};} 
	public float[] getAlphaPicFlipPtr() {return new float[] {alphaPicFlip};}      
	
	
}//myMaterial
