package MPM.BaseSim.material;

import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;

/**
 * class for the material that the simulation consists of
 * @author John Turner
 *
 */
public class MPM_Material {
	//Set as ptrs(aras) to facilitate consumption
	private float[] initYoungMod,
					poissonRatio,
					hardeningCoeff, 
					criticalCompression, 
					criticalStretch, 
					alphaPicFlip, 
					lambda0, 
					mu0;
		
	public MPM_Material(MPM_SimUpdateFromUIData _currUIVals) {
		initYoungMod = new float[1];
		poissonRatio = new float[1];
		hardeningCoeff = new float[1];
		criticalCompression = new float[1];
		criticalStretch = new float[1];
		alphaPicFlip = new float[1];
		lambda0 = new float[1];
		mu0 = new float[1];
		updateMatVals_FromUI(_currUIVals);
	}//ctor
		

	//recalc lambda0 and mu0 when young mod or poisson ratio changes
	private void recalcParams() {
		float poissonRatioP1 = (1.0f+poissonRatio[0]);
		lambda0[0] = (initYoungMod[0] * poissonRatio[0]) / (poissonRatioP1 * (1.0f-2.0f*poissonRatio[0]));
		mu0[0] = initYoungMod[0] / (2.0f * poissonRatioP1);
	}
	
	public void updateMatVals_FromUI(MPM_SimUpdateFromUIData upd) {		
		initYoungMod[0] = upd.getInitYoungMod();
		poissonRatio[0] = upd.getPoissonRatio();		
		hardeningCoeff[0]= upd.getHardeningCoeff();     
		criticalCompression[0] = upd.getCriticalCompression();
		criticalStretch[0] = upd.getCriticalStretch();   
		alphaPicFlip[0] = upd.getAlphaPicFlip();
		//recalc mu0 when ym or poisson ratio changes
		recalcParams();
	}//updateMatVals_FromUI
	
	/**
	 * Return ptr to lambda0 value (precalced from Young's modulus and Poisson ration) for kernel params
	 * @return
	 */
	public float[] getLambda0Ptr() {return lambda0;}
	/**
	 * Return ptr to Mu0 value (precalced from Young's modulus and Poisson ration) for kernel params
	 * @return
	 */
	public float[] getMu0Ptr() {return mu0;}	
	/**
	 * Return ptr to Hardening Coefficient value for kernel params
	 * @return
	 */
	public float[] getHardeningCoeffPtr() {return hardeningCoeff;}     
	/**
	 * Return ptr to Critical Compression value for kernel params
	 * @return
	 */
	public float[] getCriticalCompressionPtr() {return criticalCompression;}  
	/**
	 * Return ptr to Critical Stretch value for kernel params
	 * @return
	 */
	public float[] getCriticalStretchPtr() {return criticalStretch;} 
	/**
	 * Return ptr to Alpha PICFLIP value for kernel params
	 * @return
	 */
	public float[] getAlphaPicFlipPtr() {return alphaPicFlip;}      
	
	//gets csv string of values used by particles in grid force calculation
	public String toGridFrcCalcStrCSV(String fmt) {
		return String.format(fmt, lambda0[0])+","+String.format(fmt, mu0[0])+","+String.format(fmt, hardeningCoeff[0]);
		
	}
	public String toDefGradUpdCalcStrCSV(String fmt) {
		return String.format(fmt, criticalCompression[0])+","+String.format(fmt, criticalStretch[0]);
		
	}
}//MPM_Material
