package MPM_SimMain.material;

import MPM_SimMain.utils.MPM_SimUpdateFromUIData;

/**
 * class for the material that the simulation consists of
 * @author John Turner
 *
 */
public class MPM_Material {
	private float initYoungMod, poissonRatio, lambda0, mu0;
	private float hardeningCoeff, criticalCompression, criticalStretch, alphaPicFlip;
		
	public MPM_Material(MPM_SimUpdateFromUIData _currUIVals) {
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
	
	/**
	 * Return ptr to lambda0 value (precalced from Young's modulus and Poisson ration) for kernel params
	 * @return
	 */
	public float[] getLambda0Ptr() {return new float[] {lambda0};}
	/**
	 * Return ptr to Mu0 value (precalced from Young's modulus and Poisson ration) for kernel params
	 * @return
	 */
	public float[] getMu0Ptr() {return new float[] {mu0};}	
	/**
	 * Return ptr to Hardening Coefficient value for kernel params
	 * @return
	 */
	public float[] getHardeningCoeffPtr() {return new float[] {hardeningCoeff};}     
	/**
	 * Return ptr to Critical Compression value for kernel params
	 * @return
	 */
	public float[] getCriticalCompressionPtr() {return new float[] {criticalCompression};}  
	/**
	 * Return ptr to Critical Stretch value for kernel params
	 * @return
	 */
	public float[] getCriticalStretchPtr() {return new float[] {criticalStretch};} 
	/**
	 * Return ptr to Alpha PICFLIP value for kernel params
	 * @return
	 */
	public float[] getAlphaPicFlipPtr() {return new float[] {alphaPicFlip};}      
	
	
}//myMaterial
