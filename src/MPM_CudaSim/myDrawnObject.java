package MPM_CudaSim;

import java.util.ArrayList;
import java.util.Arrays;

import processing.core.PApplet;
import processing.core.PConstants;

public abstract class myDrawnObject {
	public static MPM_SimMain pa;

	public int fillClr, strkClr;
	public float len;								//length of object
	protected static final int numReps = 4;				//default number of repetitions of subdivid/tuck/untuck
	
	public static final int trajPtRad = 2;			//radius of points
	
	public myVector canvasNorm;							//normal to drawing canvas == normal to plane of poly
	protected myPoint[] origpts;							//originally drawn points making up this curve
	protected myPoint[] pts;								//points making up this curve
	protected float[] dpts;						//array holding distance to each point from beginning
	
	//beautiful pts
	public cntlPt[] cntlPts;						//control points describing object, if used	
	
	protected float[] d_cntlPts;
	protected float cntl_len;
	
	public myPoint COV,									//center of verts
			COM;
	//boolean flags about object
	public boolean[] flags;							//init'ed to all false, then set specifically when appropriate
	//flag idxs
	public final int 
		//construction
			isClosed 		= 0,					//object is a closed poly
			isMade 			= 1,					//whether or not the object is finished being drawn
			isFlipped		= 2,					//points being displayed are flipped (reflected)
			usesCntlPts 	= 3,					//if this curve is governed by control points (as opposed to just drawn freehand)
	    //calculation
			reCalcPoints	= 4,					//recalculate points from cntl point radii - use if radius is changed on stroke from user input
			cntlWInvRad		= 5,					//whether the weight impact on cntl radius is inverse or direct - inverse has slow drawing be wider, direct has slow drawing be narrower
			interpStroke	= 6,					//whether or not control-point based strokes are interpolating or not
		//display			 
			showCntlPnts 	= 7,					//show this object's cntl points
			vertNorms 		= 8,					//use vertex normals to shade curve
			drawNorms 		= 9,					//display normals for this object as small arrows
			drawCntlRad 	= 10,	
			useProcCurve	= 11;					//toggles whether we use straight lines in vertex building or processing's curve vertex			
	public final int numFlags = 12;					//always 1 more than last flag const
	
	//flags about type of operation that uses interpolation being done
	public int lnI_Typ,								//what interpolation type will this curve use for line operations (tuck, find myPoint at ratio of length, etc) 
				sbI_Typ,							//interp type for subdivision
				swI_Typ;							//what kind of interpolation will be used for this curve as it is swemyPoint around axis (if this is a closed sweep-poly)
	
	//flags about which interpolation type should be done
	public final int linear_int = 0,				//denotes linear interpolation 
					ball_int = 1;					//ball-morph interpolation, for sweep morph
	//add more here when we have more 
	
	public final int 	//point processing flags
		_subdivide		=0,
		_tuck			=1,
		_equaldist		=2,
		_resample		=3;

	//array of normals, tans for every control point,
	protected myVector[] c_nAra, c_tAra, c_bAra;
	
	
	protected float[][] distFromIntAxis;			//keeps the distance to internal axis of stroke			
	protected static final int 
	//idxs in distFromAxis array
			_d = 0,									//index for dist from axis
			_t = 1,									//interpolant along rot axis for proj
			_a = 2,									//idx1 for myPoint in ara for rotational axis - proj lies between these two
			_b = 3;									//idx2 for myPoint in ara for rotational axis
	protected final int nAxisVals = 4;
	
	public myDrawnObject(MPM_SimMain _pa, myVector _canvNorm) {
		pa = _pa;
		canvasNorm = _canvNorm;		//c.drawSNorm  draw surface normal
		initFlags();
		lnI_Typ = linear_int;
		swI_Typ = linear_int;
		sbI_Typ = linear_int;
		c_nAra = new myVector[0];c_tAra = new myVector[0];c_bAra = new myVector[0];
		COM = pa.P();COV = pa.P();
//		flags[firstDrawn] = true;
	}
	
	//initialize point referencing structs - using both for speed concerns.
	public void startDrawing(){	pts = new myPoint[0];len = 0; dpts = new float[0]; flags[drawCntlRad] = false;}	
	public void addPt(myPoint p){
		ArrayList<myPoint> tmp = new ArrayList<myPoint>(Arrays.asList(pts));
		tmp.add(p);
		setPts(tmp);
	}//setPt
	
	//subdivide, tuck, respace, etc, cntlpts of this curve
	public void processCntlPts(int numPts, int numReps){
		float origLen = cntl_len;
		setCPts(procCntlPt(_subdivide, cntlPts, 2, origLen));										//makes 1 extra vert  equilspaced between each vert, to increase resolution of curve
		for(int i = 0; i < numReps; ++i){
			//setCPts(procCntlPt(_subdivide, cntlPts, 2, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, .5f, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, -.5f, origLen));
		}		//smooth curve - J4
		setCPts(procCntlPt(_equaldist, cntlPts, .5f, origLen));
		for(int i = 0; i < numReps; ++i){
			//setCPts(procCntlPt(_subdivide, cntlPts, 2, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, .5f, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, -.5f, origLen));
		}		//smooth curve - J4
		setCPts(procCntlPt(_resample, cntlPts, numPts, origLen));
	}			

	//subdivide, tuck, respace, resample, etc. pts of this curve
	public void processPts(myPoint[] pts, int numPts, int numReps){
		setPts(procPts(_subdivide, pts, 2, len, flags[isClosed]));										//makes 1 extra vert  equilspaced between each vert, to increase resolution of curve
		for(int i = 0; i < numReps; ++i){
			setPts(procPts(_subdivide, pts, 2, len, flags[isClosed]));
			setPts(procPts(_tuck, pts, .5f, len, flags[isClosed]));
			setPts(procPts(_tuck, pts, -.5f, len, flags[isClosed]));
		}		//smooth curve - J4
		setPts(procPts(_equaldist, pts, .5f, len, flags[isClosed]));
		for(int i = 0; i < numReps; ++i){
			setPts(procPts(_subdivide, pts, 2, len, flags[isClosed]));
			setPts(procPts(_tuck, pts, .5f, len, flags[isClosed]));
			setPts(procPts(_tuck, pts, -.5f, len, flags[isClosed]));
		}		//smooth curve - J4
		setPts(procPts(_resample, pts, numPts, len, flags[isClosed]));		
	}	
	
	
	//sets required info for points array - points and dist between pts, length, etc
	protected void setPts(ArrayList<myPoint> tmp){
		pts = tmp.toArray(new myPoint[0]);
		dpts = getPtDist(pts, flags[isClosed]);	
		len=length(pts, flags[isClosed]);
	}//setPts	
	//make a new point interpolated between either 2 or 3 points in pts ara, described by # of idxs
	public myPoint makeNewPoint(myPoint[] pts, int[] idxs, float s){	return _Interp(pts[idxs[0]], s, (idxs.length == 2 ? pts[idxs[1]] : _Interp(pts[idxs[1]],.5f,pts[idxs[2]], lnI_Typ)),lnI_Typ );	}
	
	/**
	 * process all points using passed algorithm on passed array of points - not all args are used by all algs.
	 * @param _typ type of point processing
	 * @param pts array to be processed
	 * @param val quantity used by variou processing : subdivision-> # of new pts +1, tuck-> amt to be tucked,  resample-> # of new verts
	 * @param len length of segment described by points, including ends if closed
	 * @param wrap whether the point list wraps around or not
	 * @return arraylist of processed points
	 */
	public ArrayList<myPoint> procPts(int _typ, myPoint[] pts, float val, float _len, boolean wrap){
		ArrayList<myPoint> tmp = new ArrayList<myPoint>(); // temporary array
		switch(_typ){
			case _subdivide	:{
			    for(int i = 0; i < pts.length-1; ++i){tmp.add(pts[i]); for(int j=1;j<val;++j){tmp.add(makeNewPoint(pts,new int[]{i,i+1}, (j/(val))));}}
			    tmp.add(pts[pts.length-1]);				
			    return tmp;}
			case _tuck		:{
				if(wrap){tmp.add(makeNewPoint(pts,new int[]{0,pts.length-1,1}, val));} else {tmp.add(0,pts[0]);}
			    for(int i = 1; i < pts.length-1; ++i){	tmp.add(i,makeNewPoint(pts,new int[]{i,i-1,i+1}, val));   }
		    	if(wrap){tmp.add(makeNewPoint(pts,new int[]{pts.length-1,pts.length-2,0}, val));} else {tmp.add(pts[pts.length-1]);}			
		    	return tmp;}
			case _equaldist	:{
				float ratio = _len/(1.0f * pts.length),curDist = 0;					 //new distance between each vertex, iterative dist travelled so far			 
				for(int i =0; i<pts.length; ++i){tmp.add(at(curDist/_len));curDist+=ratio;}	
				tmp.add(pts[pts.length-1]);				
				return tmp;}	
			case _resample	:{
				float ratio = pts.length/(1.0f * (val-1)),f;					//distance between each vertex		 
				int idx, newIdx=0;		
				for(float i = 0; i<pts.length-1; i+=ratio){idx = (int)i;	f = i-idx;tmp.add(newIdx++,makeNewPoint(pts,new int[]{idx,idx+1},f));}
				if(wrap) {
					if(myPoint._dist(tmp.get(newIdx-1), tmp.get(0)) > ratio){	tmp.add(makeNewPoint(new myPoint[]{tmp.get(newIdx-1), tmp.get(0)},new int[]{0,1},.5f));}		//want to only add another point if last 2 points are further than ratio appart
				} else {		tmp.add(pts[pts.length-1]);}			//always add another point if open line/loop - want to preserve end point
				break;}	
			default :
		}
		
		return tmp;
	}
		
	//CONTROL POINT-RELATED FUNCTIONS
	//build essential orientation vectors for control points
	public void buildCntlFrameVecAras(){
		c_nAra = buildNormals(cntlPts);
		c_tAra = buildTangents(cntlPts, false);
		c_bAra = buildBinormals(cntlPts, c_nAra, c_tAra);		//use these with cntl point radius to build stroke pts
	}//buildCntlFrameVecAras
		
	//sets required info for points array - points and dist between pts, length, etc
	protected void setCPts(ArrayList<cntlPt> tmp){
		cntlPts = tmp.toArray(new cntlPt[0]);
		d_cntlPts = getPtDist(cntlPts, false);	
		cntl_len=length(cntlPts, false);
	}//setPts	
	//make a new point interpolated between either 2 or 3 points in pts ara, described by # of idxs
	public cntlPt makeNewPoint(cntlPt[] pts, int[] idxs, float s){	return _Interp(pts[idxs[0]], s, (idxs.length == 2 ? pts[idxs[1]] : _Interp(pts[idxs[1]],.5f,pts[idxs[2]], lnI_Typ)),lnI_Typ );	}
	/**
	 * process all points using passed algorithm on passed array of points - not all args are used by all algs.
	 * @param _typ type of point processing
	 * @param pts array to be processed
	 * @param val quantity used by variou processing : subdivision-> # of new pts +1, tuck-> amt to be tucked,  resample-> # of new verts
	 * @param len length of segment described by points, including ends if closed
	 * @param wrap whether the point list wraps around or not
	 * @return arraylist of processed points
	 */	
	public ArrayList<cntlPt> procCntlPt(int _typ, cntlPt[] pts, float val, float _len){
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>(); // temporary array
		switch(_typ){
			case _subdivide	:{
			    for(int i = 0; i < pts.length-1; ++i){tmp.add(pts[i]); for(int j=1;j<val;++j){tmp.add(makeNewPoint(pts,new int[]{i,i+1}, (j/(val))));}}
			    tmp.add(pts[pts.length-1]);				
			    return tmp;}
			case _tuck		:{
				tmp.add(0,pts[0]);//no wrap on control points, so no  need to check
			    for(int i = 1; i < pts.length-1; ++i){	tmp.add(i,makeNewPoint(pts,new int[]{i,i-1,i+1}, val));   }
		    	tmp.add(pts[pts.length-1]);			
		    	return tmp;}
			case _equaldist	:{
				float ratio = _len/(1.0f * pts.length),curDist = 0;					 //new distance between each vertex, iterative dist travelled so far			 
				for(int i =0; i<pts.length; ++i){tmp.add(at_C(curDist/_len, pts));curDist+=ratio;}	
				tmp.add(pts[pts.length-1]);				
				return tmp;}	
			case _resample	:{
				float ratio = pts.length/(1.0f * (val-1)),f;					//distance between each vertex		 
				int idx, newIdx=0;		
				for(float i = 0; i<pts.length-1; i+=ratio){idx = (int)i;	f = i-idx;tmp.add(newIdx++,makeNewPoint(pts,new int[]{idx,idx+1},f));}			
				tmp.add(pts[pts.length-1]);
				break;}	
			default :
		}
		return tmp;
	}	
	//end cntlmyPoint related
	//normals, tangents, binormals at each point
	public myVector[] buildNormals(myPoint[] _pts){
		ArrayList<myVector> tmp = new ArrayList<myVector>();				
		for(int i =0; i<_pts.length; ++i){tmp.add(pa.U(canvasNorm));	}		//make normal the canvas normal
		return tmp.toArray(new myVector[0]);
	}
	public myVector[] buildTangents(myPoint[] _pts, boolean close){
		ArrayList<myVector> tmp = new ArrayList<myVector>();
		for(int i=0; i<_pts.length-1; ++i){tmp.add(pa.U(_pts[i], _pts[i+1]));}
		if(close){tmp.add(pa.U(_pts[_pts.length-1], _pts[0]));}
		else {tmp.add(pa.U(_pts[_pts.length-2], _pts[_pts.length-1]));}
		return tmp.toArray(new myVector[0]);
	}	
	public myVector[] buildBinormals(myPoint[] _pts, myVector[] n_ara, myVector[] t_ara){//build last
		ArrayList<myVector> tmp = new ArrayList<myVector>();
		for(int i=0; i<_pts.length; ++i){tmp.add(pa.U(pa.N(n_ara[i],t_ara[i])));}
		return tmp.toArray(new myVector[0]);
	}

	//find location of center of verts
	public myPoint calcCOV(){myPoint C = pa.P();for(int i=0;i<pts.length;++i){C._add(pts[i]);}myPoint Ct = pa.P(1.0f/pts.length,C); COV=pa.P(Ct);return COV;}
	//find COV of passed verts
	public myPoint calcCOVOfAra(myPoint[] pts){myPoint C = pa.P();for(int i=0;i<pts.length;++i){C._add(pts[i]);}myPoint Ct = pa.P(1.0f/pts.length,C); return Ct;}

	
	/**
	 * return the interpolated myVectortor between two myPoint's myVectortors given the adjacent idx's of two points in pts and the array of myVectortors
	 * @param idxA, idxB : 2 idxs in myVector aras to be interped between
	 * @param s : interpolant
	 * @param vAra : array to find myVectors in
	 * @return interpolated myVectortor
	 */
	public myVector getInterpVec(int idxA, float s, int idxB, myVector[] vAra, int interpMech){return _Interp(vAra[idxA], s, vAra[idxB], interpMech);}
	public myVector getUInterpVec(int idxA, float s, int idxB, myVector[] vAra, int interpMech){return pa.U(_Interp(vAra[idxA], s, vAra[idxB], interpMech));}
	
	/**
	 * using arc length parameterisation this will return a point along the curve at a 
	 * particular fraction of the length of the curve (0,1 will return endpoints, .5 will return halfway along curve)
	 * @param t fraction of curve length we are interested in returning a point - should be 0-1
	 * @return point @t along curve
	 */
	public myPoint at(float t){return at(t,new float[1], len, pts, dpts);}//put interpolant between adjacent axis points in s ara if needed
	public myPoint at(float t, float[] s){return at(t,s, len, pts, dpts);}//put interpolant between adjacent axis points in s ara if needed
	public myPoint at(float t, float[] s, float _len, myPoint[] pts, float[] _dpts){//call directly if wanting interpolant between adj axis points too
		if(t<0){PApplet.println("In at : t="+t+" needs to be [0,1]");return pts[0];} else if (t>1){PApplet.println("In at : t="+t+" needs to be [0,1]");return pts[pts.length-1];}
		float dist = t * _len;
		for(int i=0; i<_dpts.length-1; ++i){										//built off dpts so that it will get wrap for closed curve
			if(_dpts[i+1] >= dist){
				s[0] = ((dist-_dpts[i])/(_dpts[i+1]-_dpts[i]));					//needs to stay between 0 and 1 (since interpolation functions between pts will be 0-1 based), so normalize by distance dpts[i]
				return makeNewPoint(pts,new int[]{i,((i+1)%pts.length)}, s[0]);		//put interpolant between adjacent axis points in s ara if needed		
			}					
		}		
		return pts[0];
	}//at	
	
	public cntlPt at_C(float t, cntlPt[] pts){float[] _dpts = this.getPtDist(pts, false);float _len = this.length(pts, false);return at_C(t,new float[1], _len, pts, _dpts);}//put interpolant between adjacent axis points in s ara if needed
	public cntlPt at_C(float t, float[] s, float _len, cntlPt[] pts, float[] _dpts){//call directly if wanting interpolant between adj axis points too
		if(t<0){PApplet.println("In at : t="+t+" needs to be [0,1]");return pts[0];} else if (t>1){PApplet.println("In at : t="+t+" needs to be [0,1]");return pts[pts.length-1];}
		float dist = t * _len;
		for(int i=0; i<_dpts.length-1; ++i){										//built off dpts so that it will get wrap for closed curve
			if(_dpts[i+1] >= dist){
				s[0] = ((dist-_dpts[i])/(_dpts[i+1]-_dpts[i]));					//needs to stay between 0 and 1 (since interpolation functions between pts will be 0-1 based), so normalize by distance dpts[i]
				return makeNewPoint(pts,new int[]{i,((i+1)%pts.length)}, s[0]);		//put interpolant between adjacent axis points in s ara if needed		
			}			
		}		
		return new cntlPt(pa);
	}//at_C	
	
	/**
	 * this will conduct the appropriate interpolation on the two passed points, based on what type is set for the interpolation requested
	 * as more interpolation schemes are implemented, add cases
	 * @param A, B : points to be interpolated
	 * @param s : interpolant
	 * @param _typ : what is being interpolated - smoothing a line, sweeping the curve around the axis, etc.  for each _typ, a value should be set for this curve (i.e.smInterpTyp)
	 * @return : resultant point
	 */
	protected myPoint _Interp(myPoint A, float s, myPoint B, int _typ){
		switch (_typ){
			case linear_int : {	return pa.L(A, s, B);}
			//add more cases for different interpolation		
			default : {	return pa.L(A, s, B);}			//defaults to linear
		}	
	}//_Interp
	//same as above, but with myVectortors
	protected myVector _Interp(myVector A, float s, myVector B, int _typ){
		switch (_typ){
			case linear_int : {	return pa.V(A, s, B);}
			//add more cases for different interpolation		
			default : {	return pa.V(A, s, B);}			//defaults to linear
		}	
	}//_Interp
	//same as above but with doubles
	protected float _Interp(float A, float s, float B, int _typ){
		switch (_typ){
			case linear_int : {	return (1-s)*A + (s*B);}
			//add more cases for different interpolation		
			default : {	return (1-s)*A + (s*B);}			//defaults to linear
		}	
	}//_Interp
	protected cntlPt _Interp(cntlPt A, float s, cntlPt B, int _typ){
		switch (_typ){
			case linear_int : {	return cntlPt.L(A, s, B);}
			//add more cases for different interpolation		
			default : {	return cntlPt.L(A, s, B);}			//defaults to linear
		}	
	}//_Interp	

	//draw currently selected point
	public void drawSelPoint(int i ){
		pa.pushMatrix();		pa.pushStyle();
		pa.stroke(255,255,0,255);
		if(flags[usesCntlPts]){pa.show(cntlPts[i], 3);} else {pa.show(pts[i], 3);}
		pa.popStyle();		pa.popMatrix();
	}
	

	public abstract void rebuildPolyPts();

	//makes a copy of the points in order
	public myPoint[] cpyPoints(myPoint[] pts){myPoint[] tmp = new myPoint[pts.length]; for(int i=0; i<pts.length; ++i){	tmp[i]=pa.P(pts[i]);}	return tmp;}//cpyPoints
	//makes a copy of the points in order
	public cntlPt[] cpyPoints(cntlPt[] pts){cntlPt[] tmp = new cntlPt[pts.length]; for(int i=0; i<pts.length; ++i){	tmp[i]=new cntlPt(pts[i]);}	return tmp;}//cpyPoints

	//move points by the passed myVectortor 
	public myPoint[] movePoints(myVector move, myPoint[] pts){for(int i =0; i<pts.length; ++i){	pts[i]._add(move);	}	return pts;}	
	//move points by the passed myVectortor 
	public cntlPt[] movePoints(myVector move, cntlPt[] pts){for(int i =0; i<pts.length; ++i){	pts[i]._add(move);	}	return pts;}	
	//flip the passed points and move them based on the displacement from the passed movement myVectortor
//	public myPoint[] flipPtsAndMove(myDrawnObject _obj, myPoint[] pts, myVector move,  myVector covAxis){
//		myPoint[] tmp = movePoints(move, cpyPoints(pts));
//		return tmp;
//	}//flipPtsAndMove
//	//flip the passed cntlpts and move them based on the displacement from the passed movement myVectortor
//	public cntlPt[] flipPtsAndMove(myDrawnObject _obj, cntlPt[] pts, myVector move,  myVector covAxis, boolean reverse){
//		cntlPt[] tmp = movePoints(move, (reverse ? rotPtsAroundCOV(pts, PApplet.PI, _obj.COV, pa.U(_obj.canvasNorm), covAxis) : cpyPoints(pts)));
//		return tmp;
//	}//flipPtsAndMove
	//set this object's points to be passed points, for copying
	public void setPtsToArrayList(cntlPt[] pts){ArrayList<cntlPt> tmp = new ArrayList<cntlPt>(Arrays.asList(pts));setCPts(tmp);}	
	//set this object's points to be passed points, for copying
	public void setPtsToArrayList(myPoint[] pts){ArrayList<myPoint> tmp = new ArrayList<myPoint>(Arrays.asList(pts));setPts(tmp);}	

	
	//rotate points around axis that is xprod of canvas norm and the lstroke cov to the rstroke cov myVector at stroke cov.
	//should make mirror image of pts
	public myPoint[] rotPtsAroundCOV(myPoint[] pts, float angle, myPoint cov, myVector _canvasNorm, myVector _covNorm){//need to pass canvas norm since it cannot be static
		ArrayList<myPoint> tmp = new ArrayList<myPoint>();//				res.get(sl)[i] = pa.R(old, sliceA, canvasNorm, bv, myPointOnAxis);
		for(int i=0; i<pts.length; ++i){tmp.add(pa.R(pts[i], angle, pa.V(_canvasNorm), _covNorm, cov));}
		return tmp.toArray(new myPoint[0]);			
	}//	
	//rotate points around axis that is xprod of canvas norm and the lstroke cov to the rstroke cov myVector at stroke cov.
	//should make mirror image of pts
	public cntlPt[] rotPtsAroundCOV(cntlPt[] pts, float angle, myPoint cov, myVector _canvasNorm, myVector _covNorm){//need to pass canvas norm since it cannot be static
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>();//				res.get(sl)[i] = pa.R(old, sliceA, canvasNorm, bv, myPointOnAxis);
		for(int i=0; i<pts.length; ++i){tmp.add(pa.R(pts[i], angle, pa.V(_canvasNorm), _covNorm, cov));}
		return tmp.toArray(new cntlPt[0]);			
	}//	
	//finds index of point with largest projection on passed myVectortor in passed myPoint ara
	protected int findLargestProjection(myVector v, myPoint c, myPoint[] pts){
		double prjLen = -1, d;
		int res = -1;
		for(int i=0; i<pts.length; ++i){d = myVector._dot(v,pa.V(c, pts[i]));if(d > prjLen){prjLen = d;res = i;}}	
		return res;
	}//findLargestProjection : largest projection on passed myVectortor in passed myPoint ara		
	
	public abstract int findClosestPt(myPoint p, double[] d);
	
	//finds closest point to p in sPts - put dist in d
	protected final int findClosestPt(myPoint p, double[] d, myPoint[] _pts){
		int res = -1;
		double mindist = 99999999, _d;
		for(int i=0; i<_pts.length; ++i){_d = myPoint._dist(p,_pts[i]);if(_d < mindist){mindist = _d; d[0]=_d;res = i;}}	
		return res;
	}
	//reorder verts so that they start at newStIdx - only done for closed  meshes, so should wrap around
	//made static so doesn't have to be used with pts array - should be called by implementing class
	protected ArrayList<myPoint> reorderPts(int newStIdx, myPoint[] pts){
		ArrayList<myPoint> tmp = new ArrayList<myPoint>(Arrays.asList(pts));
		for(int i = 0; i<pts.length; ++i){tmp.set(i, pts[(i+newStIdx)%pts.length]);	}
		return tmp;
	}		
	
	public abstract void remakeDrawnTraj(boolean useVels);

	//run at the end of drawing a curve - will set appropriate flags, execute smoothing, subdivision, resampling, etc
	public abstract void finalizeDrawing(boolean procPts);	

	public void drawMe() {
		pa.pushMatrix();
		pa.pushStyle();
		pa.setColorValFill(fillClr);
		pa.setColorValStroke(strkClr);
			pa.strokeWeight(1);
//			if(flags[useProcCurve]){pa.show(pts);} 
//			else {			
				pa.curve(pts);
				//}
		pa.popStyle();			
		pa.popMatrix();
//		if(flags[drawNorms] && (nAra != null)&& (tAra != null)&& (bAra != null)){drawNorms(pts, nAra,tAra,bAra,20);}
//		drawCOV();
	}
	public abstract void dragPicked(myVector disp, int idx);
	public abstract void dragAll(myVector disp);			//move COV to location pointed at by mouse, and move all points by same displacement
	//drag point represented by passed idx in passed array - either point or cntl point
	protected void dragPicked(myVector dispVec, int idx, myPoint[] pts) {if((-1 != idx) && (pts[idx] != null)){pts[idx]._add(dispVec);flags[reCalcPoints]=true;}}
	protected void dragAll(myVector dispVec, myPoint[] pts){if((pts == null)||(pts.length == 0)){return;}for(int i=0;i<pts.length;++i){pts[i]._add(dispVec);}flags[reCalcPoints]=true;}
	//rotate points around axis that is xprod of canvas norm and the lstroke cov to the rstroke cov vec at stroke cov.
	//should make mirror image of pts
	public myPoint[] rotPtsAroundCOV(myPoint[] _pts, double angle, myPoint cov, myVector _canvasNorm, myVector _covNorm){//need to pass canvas norm since it cannot be static
		ArrayList<myPoint> tmp = new ArrayList<myPoint>();//				res.get(sl)[i] = pa.R(old, sliceA, canvasNorm, bv, ptOnAxis);
		for(int i=0; i<_pts.length; ++i){tmp.add(pa.R(_pts[i], angle, pa.V(_canvasNorm), _covNorm, cov));}
		return tmp.toArray(new myPoint[0]);			
	}//	
	//rotate points around axis that is xprod of canvas norm and the lstroke cov to the rstroke cov vec at stroke cov.
	//should make mirror image of pts
	public cntlPt[] rotPtsAroundCOV(cntlPt[] _pts, double angle, myPoint cov, myVector _canvasNorm, myVector _covNorm){//need to pass canvas norm since it cannot be static
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>();//				res.get(sl)[i] = pa.R(old, sliceA, canvasNorm, bv, ptOnAxis);
		for(int i=0; i<_pts.length; ++i){tmp.add(pa.R(_pts[i], angle, pa.V(_canvasNorm), _covNorm, cov));}
		return tmp.toArray(new cntlPt[0]);			
	}//		
	
	//returns array of distances to each point from beginning - needs to retain dist from last vert to first if closed
	//public final float[] getPtDist(){float[] res = new float[pts.length+1];res[0]=0;for(int i=1; i<pts.length; ++i){res[i] = res[i-1] + pa.d(pts[i-1],pts[i]);}if(flags[isClosed]){res[pts.length] = res[pts.length-1] +  pa.d(pts[pts.length-1],pts[0]);}return res;}
	//public final float[] getPtDist(){return getPtDist(pts, flags[isClosed]);}
	//returns array of distances to each point from beginning - needs to retain dist from last vert to first if closed
	public final float[] getPtDist(myPoint[] pts, boolean wrap){
		float[] res = new float[pts.length+1];
		res[0]=0;
		for(int i=1; i<pts.length; ++i){
			//System.out.println("i : "+i);
			res[i] = res[i-1] + (float)myPoint._dist(pts[i-1],pts[i]);
			}
		if(wrap){
			//System.out.println("wrap");
			res[pts.length] = res[pts.length-1] +  (float)myPoint._dist(pts[pts.length-1],pts[0]);
		} else {
			//System.out.println("no wrap");
			
			res[pts.length] = 0;
		}
		
		return res;}
	//returns length of curve, including endpoint if closed
	//public final float length(){return length(pts, flags[isClosed]);}//{float res = 0;for(int i =0; i<pts.length-1; ++i){res += pa.d(pts[i],pts[i+1]);}if(flags[isClosed]){res+=pa.d(pts[pts.length-1],pts[0]);}return res;}
	public final float length(myPoint[] pts, boolean closed){float res = 0;for(int i =0; i<pts.length-1; ++i){res += (float)myPoint._dist(pts[i],pts[i+1]);}if(closed){res+=(float)myPoint._dist(pts[pts.length-1],pts[0]);}return res;}


	public void drawCOV(){		if(COV == null) {return;}		pa.pushMatrix();		pa.pushStyle();	pa.stroke(255,0,255,255);		pa.show(COV,3);		pa.popStyle();		pa.popMatrix();	}
	//drawCntlRad
	public myPoint getPt(int i){return pts[i];}
	
	public void initFlags(){flags = new boolean[numFlags]; for(int i=0;i<numFlags;++i){flags[i]=false;}}
	public String toString(){
		String res = "#pts : "+pts.length+" len : "+ len ;
		return res;
	}

}//myDrawnObject


class myVariStroke extends myDrawnObject {

	protected final int numVerts = 200;							

	public int offsetType;						//1:Q-bspline w/normal offset, 2:Q-bspline w/ball offset, 3:Q-bspline w/radial offset
	public myOffset _offset;					//offset used by this stroke to calculate poly loop
	public final int numIntCntlPts = 200, numCntlPts = 6;			//# of control points to use to draw line
	private boolean ptsDerived;					//whether or not the points making up the loop of this stroke have been derived yet

	
	public myPoint[] drawnCntlPts;						//control point locations - start point +  vel myVectortor (scaled tangent), repeat
	public float[] d_drwnCntlPts;					//distance between each drawnCntlPts points
	public float drwnCntlLen;						//len of arc of drawnCntlPts pts

	public float[] cntlPtIntrps;				//interpolants for each control point as part of total, based upon radius of controlpoint - larger will be bigger. 
												//for each cntl myPoint, this value will be cntl point rad / total control point rads.
	//interpolated drawn curve, weighted by drawn speed
	public myPoint[] interpCntlPts;					//interpolated control point locations - start point +  vel myVectortor (scaled tangent), repeat
	public float[] d_interpPts;					//distance between interpolated points
	public float interpLen;						//len of interp pts
	
	public myVariStroke(MPM_SimMain _pa, myVector _canvNorm, int fillClrCnst, int strkClrCnst) {
		super(_pa, _canvNorm);
		flags[isClosed] = false;	
		fillClr = fillClrCnst;
		strkClr= strkClrCnst;
		flags[drawCntlRad] = true;
	    cntlPts = new cntlPt[0];
	    interpCntlPts = new myPoint[0];
		_offset = new myNormOffset(_pa);
		flags[usesCntlPts] = true;
		flags[interpStroke] = true;
		ptsDerived = false;
		flags[cntlWInvRad] = false;			//whether slow drawing makes rad larger or smaller
	}

	//as drawing, add points to -cntlPoints-, not pts array.
	public void addPt(myPoint p){
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>(Arrays.asList(cntlPts));
		int i = tmp.size()-1;
		if(i > 0 ){tmp.get(i).w = calcCntlWeight(p,tmp.get(i),tmp.get(i-1));}//previous point's weight 
		cntlPt tmpPt = new cntlPt(pa, p);
		tmp.add(tmpPt);
		setCPts(tmp);
	}//
	
	public void finalizeDrawing(boolean procPts){
		//by here we have the drawn points from user input we want use the offset to determine the actual points of the curve we want to put this in a function so that any changes to the 
		//cntlpoints can be cascaded down to the actual loop		
		buildPointsUsingOffset(procPts, numReps);
		//calculate line points from control points
		//find loop around stroke line by cntl points' radii once loop is built, treat as poly loop
		flags[isMade] = true;
		flags[drawCntlRad] = false;
		//build array of weights from 
		buildInterpAra();
	}//finalize
	
	//build pts array using cntlpoints and offset chosen
	public void buildPointsUsingOffset(boolean procPts, int repCnt){
		if(procPts){
		    finalizeCntlW();
		    for(int i=0;i<cntlPts.length;++i){cntlPts[i].calcRadFromWeight(cntl_len/cntlPts.length, flags[cntlWInvRad]);}           //sets all radii based on weights
		    processCntlPts(flags[interpStroke] ? numIntCntlPts : numCntlPts, repCnt);
	    }
		buildCntlFrameVecAras();
		buildPtsFromCntlPts();
	}
	//build interpolants based upon weight of each cntl point.  these can be used to determine how far the edge moves (velocity) and how much it rotates per frame
	//when this is done, we have multipliers to apply to each tangent myVectortor to determine displacement for each of the frames
	public void buildInterpAra(){
		cntlPtIntrps = new float[numIntCntlPts];			//interpolant from each control point - weights
		
		float[] cmyPointInterps = new float[numIntCntlPts];
		interpCntlPts = new myPoint[numIntCntlPts];
		drawnCntlPts = new myPoint[numIntCntlPts];
		float sumWts = 0;
		for(int i=0;i<cntlPts.length;++i){sumWts += cntlPts[i].w;cmyPointInterps[i] = cntlPts[i].w;cntlPtIntrps[i]=sumWts;}
		//for(int i=0;i<cntlPts.length;++i){sumWts += cntlPts[i].w;cntlPtIntrps[i]=cntlPts[i].w;}
		//System.out.println("total weight = " + sumWts);
		for(int i=0;i<cntlPts.length;++i){cntlPtIntrps[i]/=sumWts;cmyPointInterps[i]/=sumWts;}
		//smooth interpolants now
		cntlPtIntrps = dualFloats(cntlPtIntrps);
		cmyPointInterps = dualFloats(cmyPointInterps);
		
		interpCntlPts[0] = pa.P(cntlPts[0]);			//set first point
		drawnCntlPts[0] = pa.P(cntlPts[0]);
		double distStToEnd = myPoint._dist(cntlPts[0], cntlPts[cntlPts.length-1]);			//distance from first to last point
		
		for(int i=1;i<cntlPts.length;++i){
			interpCntlPts[i] = pa.P(interpCntlPts[i-1],pa.W(distStToEnd * cmyPointInterps[i], c_tAra[i]));		
			drawnCntlPts[i] = pa.P(cntlPts[i]);
		}	
		for(int i= cntlPts.length; i<numIntCntlPts; ++i){
			interpCntlPts[i] = pa.P(interpCntlPts[i-1],pa.W(distStToEnd * cmyPointInterps[i], c_tAra[c_tAra.length-1]));		
			drawnCntlPts[i] = pa.P(cntlPts[cntlPts.length-1]);
			
		}
		d_interpPts = getPtDist(interpCntlPts, false);	
		interpLen = length(interpCntlPts, false);
		d_drwnCntlPts = getPtDist(drawnCntlPts, false);	
		drwnCntlLen = length(drawnCntlPts, false);		
		//smooth/equi-space interpolated cntl points
		processInterpPts( numIntCntlPts, numReps);
	}	
	//subdivide, tuck, respace, resample, etc. pts of this curve
	public void processInterpPts(int numPts, int numReps){
		//setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));										//makes 1 extra vert  equilspaced between each vert, to increase resolution of curve
		for(int i = 0; i < numReps; ++i){
			if(i % 2 == 0){setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));}
			setInterpPts(procInterpPts(_tuck, interpCntlPts, .5f, interpLen));
			setInterpPts(procInterpPts(_tuck, interpCntlPts, -.5f, interpLen));
		}		//smooth curve - J4
		setInterpPts(procInterpPts(_equaldist, interpCntlPts, .5f, interpLen));
		for(int i = 0; i < numReps; ++i){
			//if(i % 2 == 0){setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));}
			setInterpPts(procInterpPts(_tuck, interpCntlPts, .5f, interpLen));
			setInterpPts(procInterpPts(_tuck, interpCntlPts, -.5f, interpLen));
		}		//smooth curve - J4
		setInterpPts(procInterpPts(_resample, interpCntlPts, numPts, interpLen));	
	}	
	
	//return appropriate ara of points based on using velocities or not
	public myPoint[] getDrawnPtAra(boolean useVels){	return (useVels ? interpCntlPts : cntlPts);}
	
	//remake drawn trajectory after edit
	public void remakeDrawnTraj(boolean useVels){
		if(useVels){
			for(int i = 0; i < 10; ++i){
				if(i % 2 == 0){setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));}
				setInterpPts(procInterpPts(_tuck, interpCntlPts, .5f, interpLen));
				setInterpPts(procInterpPts(_tuck, interpCntlPts, -.49f, interpLen));
			}		//smooth curve - J4
			setInterpPts(procInterpPts(_equaldist, interpCntlPts, .5f, interpLen));
			setInterpPts(procInterpPts(_resample, interpCntlPts, numIntCntlPts, interpLen));	
		} else {
			for(int i = 0; i < 10; ++i){
				//setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));
				setCPts(procCntlPt(_tuck, cntlPts, .5f, cntlPts.length));
				setCPts(procCntlPt(_tuck, cntlPts, -.49f, cntlPts.length));
			}		//smooth curve - J4
			setCPts(procCntlPt(_equaldist, cntlPts, .5f, cntlPts.length));
			setCPts(procCntlPt(_resample, cntlPts, numIntCntlPts, cntlPts.length));
		}	
	}//remakeDrawnTraj
	
	//modify curve via editing
	public void handleMouseDrag(myVector dispVec, int drawnTrajPickedIdx){		
		if((drawnTrajPickedIdx == 0) || (drawnTrajPickedIdx == pts.length-1)) {return;}	//rough bounds checking
		myPoint[] pts = getDrawnPtAra(false);
		int minBnd = PApplet.max(drawnTrajPickedIdx - pa.drawnTrajEditWidth, 0),
			maxBnd = PApplet.min(drawnTrajPickedIdx + pa.drawnTrajEditWidth, pts.length-1);		
		//System.out.println("drag in drag zone inside disp calc -> idx bounds : " + minBnd + " | " + maxBnd);
		float modAmt, invdistLow = 1.0f/(drawnTrajPickedIdx - minBnd), invdistHigh = 1.0f/(maxBnd - drawnTrajPickedIdx);
		for(int idx = minBnd; idx < maxBnd; ++idx){
			float divMultVal = (idx > drawnTrajPickedIdx) ? invdistHigh:invdistLow;
			modAmt = pa.trajDragScaleAmt* PApplet.cos((idx-drawnTrajPickedIdx) * PConstants.HALF_PI * divMultVal);//trajDragScaleAmt/abs(1 + (idx-drawnTrajPickedIdx));
			//modAmt *= modAmt;
			pts[idx]._add(pa.W(modAmt,dispVec));
		}
	}
	//scale points to be a scaleAmt * current distance from line of myPoint a -> myPoint b
	public void scalePointsAboveAxis(myPoint a, myPoint b, myVector perpVec, double scaleAmt){
		myPoint[] pts = getDrawnPtAra(false);//, newPts = new myPoint[pts.length];\
		int numPoints = pts.length;
//		if((Double.isNaN(pts[pts.length-1].x)) || 
//				(Double.isNaN(pts[pts.length-1].y)) ||
//				(Double.isNaN(pts[pts.length-1].z))){
//			//pa.outStr2Scr("NaN pts size at start of scalePointsAboveAxis",true);
//			numPoints--;										//toss last NaN Point
//		}
		myPoint[] newPts = new myPoint[numPoints];
		double dist;
		//pa.outStr2Scr("cntlPts size at scalePointsAboveAxis : " + pts.length,true);
		for(int i =0; i<numPoints; ++i){
			dist = pa.distToLine(pts[i], a,b);
			//if(Double.isNaN(dist)){dist = 0;}
			myPoint pointOnLine = pa.projectionOnLine(pts[i], a,b);
			myVector resVec = pa.V(dist*scaleAmt,pa.U(pointOnLine,pts[i]));
			//pa.outStr2Scr("cntlPts : st : dist*scale : "+ (dist*scaleAmt)+" dist : "+ (dist)+" scale : "+ (scaleAmt)+" stPoint : " + pts[i].toStrBrf() + " | linePt : " + pointOnLine.toStrBrf() + " | resVec : " +resVec.toStrBrf() ,true);
			pts[i].set(pa.P(pointOnLine, resVec));
		}
		//pa.outStr2Scr("cntlPts size at scalePointsAboveAxis end : " + pts.length,true);
//		for(int i =0; i<newPts.length; ++i){
//			pts[i].set(newPts[i]);
//		}
	}
	
	//sets required info for points array - points and dist between pts, length, etc
	protected void setInterpPts(ArrayList<myPoint> tmp){
		interpCntlPts = tmp.toArray(new myPoint[0]);
		d_interpPts = getPtDist(interpCntlPts, false);	
		interpLen=length(interpCntlPts, false);
	}//setPts	
	//tuck untuck float values
	public float[] dualFloats(float[] src){
		float[] res = new float[src.length],res1 = new float[src.length];
		res1[0]=src[0];
		res1[src.length-1]=src[src.length-1];
		for(int i=1; i<src.length-1;++i){
			res1[i]=_Interp(src[i],.5f,_Interp(src[i-1],.5f,src[i+1],lnI_Typ),lnI_Typ);
		}
		res[0]=res1[0];
		res[src.length-1]=res1[src.length-1];
		for(int i=1; i<res1.length-1;++i){
			res[i]=_Interp(res1[i],-.5f,_Interp(res1[i-1],.5f,res1[i+1],lnI_Typ),lnI_Typ);
		}			
		return res;		
	}
	
	/**
	 * process all points using passed algorithm on passed array of points - not all args are used by all algs.
	 * @param _typ type of point processing
	 * @param _pts array to be processed
	 * @param val quantity used by variou processing : subdivision-> # of new pts +1, tuck-> amt to be tucked,  resample-> # of new verts
	 * @param len length of segment described by points, including ends if closed
	 * @return arraylist of processed points
	 */
	public ArrayList<myPoint> procInterpPts(int _typ, myPoint[] _pts, float val, float _len){
		ArrayList<myPoint> tmp = new ArrayList<myPoint>(); // temporary array
		switch(_typ){
			case _subdivide	:{
			    for(int i = 0; i < _pts.length-1; ++i){tmp.add(_pts[i]); for(int j=1;j<val;++j){tmp.add(makeNewPoint(_pts,new int[]{i,i+1}, (j/(val))));}}
			    tmp.add(_pts[_pts.length-1]);				
			    return tmp;}
			case _tuck		:{
				tmp.add(0,_pts[0]);
			    for(int i = 1; i < _pts.length-1; ++i){	tmp.add(i,makeNewPoint(_pts,new int[]{i,i-1,i+1}, val));   }
		    	tmp.add(_pts[_pts.length-1]);		
		    	return tmp;}
			case _equaldist	:{
				float ratio = _len/(1.0f * _pts.length),curDist = 0;					 //new distance between each vertex, iterative dist travelled so far			 
				for(int i =0; i<_pts.length; ++i){tmp.add(at_I(curDist/_len));curDist+=ratio;}	
				tmp.add(_pts[_pts.length-1]);				
				return tmp;}	
			case _resample	:{
				float ratio = _pts.length/(1.0f * (val-1)),f;					//distance between each vertex		 
				int idx, newIdx=0;		
				for(float i = 0; i<_pts.length-1; i+=ratio){idx = (int)i;	f = i-idx;tmp.add(newIdx++,makeNewPoint(_pts,new int[]{idx,idx+1},f));}
				tmp.add(_pts[_pts.length-1]);			//always add another point if open line/loop - want to preserve end point
				break;}	
			default :
		}
		
		return tmp;
	}

	public myPoint at_I(float t){return at(t,new float[1], interpLen, interpCntlPts, d_interpPts);}//put interpolant between adjacent points in s ara if needed
	public myPoint at_I(float t, float[] s){	return at(t,s, interpLen, interpCntlPts, d_interpPts);}//put interpolant between adjacent points in s ara if needed
	
	
	private void buildPtsFromCntlPts(){
		ArrayList<myPoint> tmp =  _offset.calcOffset(cntlPts, c_bAra, c_tAra) ;
		pts = tmp.toArray(new myPoint[0]);		
		ptsDerived = true;		
	}//buildPtsFromCntlPts

	//calculate initial weights for last few points of drawn cntlpoint line
	private void finalizeCntlW(){
		if(cntlPts.length == 0){return ;}
		cntlPts[0].w = cntlPts.length < 2 ? 1 : cntlPts[1].w;
		if(cntlPts.length == 1){return ;}
		cntlPts[cntlPts.length-2].w = cntlPts.length < 3 ? cntlPts[0].w : calcCntlWeight(cntlPts[cntlPts.length-3],cntlPts[cntlPts.length-2],cntlPts[cntlPts.length-3]);
		cntlPts[cntlPts.length-1].w = cntlPts.length < 2 ? cntlPts[0].w : cntlPts[cntlPts.length-2].w;
	}
	
	//build poly loop points using offsets from control points radii
	public void rebuildPolyPts(){
		flags[reCalcPoints]=false;
		//buildPointsUsingOffset(true,1);
		buildPtsFromCntlPts();
		flags[isMade] = true;
	}
	
	//print out all trajectory point locations for debugging
	public void dbgPrintAllPoints(boolean useDrawnVels){
    	if(useDrawnVels){
    		pa.outStr2Scr("Drawn Vels Traj :\n");
			for(int i = 0; i < interpCntlPts.length; ++i){
				pa.outStr2Scr("\tpt " + i +" : " + interpCntlPts[i].toStrBrf());
			}
    	} else {			
       		pa.outStr2Scr("Drawn Traj :\n");
			for(int i = 0; i < cntlPts.length; ++i){
				pa.outStr2Scr("\tpt " + i +" : " + cntlPts[i].toStrBrf());
			}
    	}
	}
	
	public void drawMe(boolean useDrawnVels, boolean flat){
		pa.pushMatrix();
		pa.pushStyle();
			pa.setColorValFill(fillClr);
			pa.setColorValStroke(strkClr);
			pa.strokeWeight(1);
        	if(useDrawnVels){
        		int clrInt = 0;
    			for(int i = 0; i < interpCntlPts.length; ++i){
    	        	clrInt = (int)(i/(1.0f * interpCntlPts.length) * 255.0f);
    	            //pa.fill(clrInt,255,(255 - clrInt),255);  
    	            //pa.stroke(clrInt,255,(255 - clrInt),255); 
    				pa.show(interpCntlPts[i],trajPtRad,-1,-1, flat);
    				if(flags[drawCntlRad]){pa.circle(this.interpCntlPts[i], this.cntlPts[i].r,this.c_bAra[i], this.c_tAra[i],20);}
    			}
        	} else {			
				for(int i = 0; i < cntlPts.length; ++i){
					pa.show(cntlPts[i],trajPtRad,fillClr,strkClr, flat);
				}
				if(flags[drawCntlRad]){this._offset.drawCntlPts(this.cntlPts, this.c_bAra, this.c_tAra, ptsDerived);}
        	}
		pa.popStyle();			
		pa.popMatrix();		
	}//
	
	//scale the points - for when the window is resized
	public void scaleMeY(boolean useDrawnVels, float scAmtY, float borderOffset){
		myPoint[] pts = getDrawnPtAra(useDrawnVels);
  			for(int i = 0; i < pts.length; ++i){
				pts[i].y -= borderOffset;
				pts[i].y *= scAmtY;
				pts[i].y += borderOffset;
    	}
	}//scaleMeY
	
	public myPoint[] moveVelCurveToEndPoints(myPoint startPt, myPoint endPt, boolean flip){
		int numPoints = interpCntlPts.length;

		myPoint[] destCurve = new myPoint[numPoints];

		if(numPoints == 0){return destCurve;}
		myPoint origin = interpCntlPts[0];
		myPoint end = interpCntlPts[numPoints - 1];

		//edge params		
		myVector drawnAxis = pa.V(origin, end);
		myVector edgeAxis =  pa.V(startPt, endPt);		//angle between these two is the angle to rotate everyone
		
		//transformation params
		myVector dispToStart = pa.V(origin, startPt);			//displacement myVectortor between start of drawn curve and edge 1.

		double alpha =  -pa.angle(drawnAxis,edgeAxis);			//angle to rotate everyone
		double scaleRatio = edgeAxis._mag()/drawnAxis._mag();	//ratio of distance from start to finish of drawn traj to distance between edges - multiply all elements in drawn traj by this
	
		//displace to align with start
		destCurve = movePoints(dispToStart, interpCntlPts);
		//build displacement vector - flip across axis if flipped curve
		myVector[] dispVecAra = new myVector[numPoints];
		dispVecAra[0] = new myVector();
		for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
			dispVecAra[myPointItr] = pa.V(destCurve[0],destCurve[myPointItr]);
		}			
		if((flip) || flags[isFlipped]){
			myVector udAxis = pa.U(drawnAxis);
			myVector normPt, tanPt;
			for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
				tanPt = myVector._mult(udAxis, dispVecAra[myPointItr]._dot(udAxis));			//component in udAxis dir
				normPt = myVector._sub(dispVecAra[myPointItr],tanPt);
				normPt._mult(2);
				dispVecAra[myPointItr]._sub(normPt);
			}
			flags[isFlipped] = flip;
		}

		//displace every point to be scaled distance from start of curve equivalent to scale of edge distances to drawn curve
		for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
			destCurve[myPointItr].set(pa.P(destCurve[0],scaleRatio, dispVecAra[myPointItr]));//start point displaced by scaleRatio * myVectortor from start to original location of myPoint
		}
		//rotate every point around destCurve[0] by alpha
		for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
			destCurve[myPointItr] = pa.R(destCurve[myPointItr],alpha,this.c_bAra[0], this.c_tAra[0],destCurve[0]);
		}
		//pa.outStr2Scr("interpCntlPts size : " + interpCntlPts.length,true);
		interpCntlPts = destCurve;
		//pa.outStr2Scr("interpCntlPts size : " + interpCntlPts.length,true);
		return destCurve;
	}//
		
	public cntlPt[] moveCntlCurveToEndPoints(myPoint startPt,myPoint endPt, boolean flip){
		int numPoints = cntlPts.length;
//		if((Double.isNaN(cntlPts[cntlPts.length-1].x)) || 
//				(Double.isNaN(cntlPts[cntlPts.length-1].y)) ||
//				(Double.isNaN(cntlPts[cntlPts.length-1].z))){
//			//pa.outStr2Scr("NaN cntlPts size at start ",true);
//			numPoints--;										//toss last NaN Point
//		}
		cntlPt[] destCurve = new cntlPt[numPoints];
		//pa.outStr2Scr("cntlPts size at start : " + cntlPts.length,true);
		//pa.outStr2Scr("first and last cntlPoint : " + cntlPts[0].toStrBrf() + " | " + cntlPts[numPoints-1].toStrBrf() );
		//drawn curve params
		if(numPoints == 0){return destCurve;}
		myPoint origin = cntlPts[0], end = cntlPts[numPoints - 1];
		//edge params		
		myVector drawnAxis = pa.V(origin, end);
		myVector edgeAxis =  pa.V(startPt, endPt);		//angle between these two is the angle to rotate everyone
		
		//transformation params
		myVector dispToStart = pa.V(origin, startPt);			//displacement myVectortor between start of drawn curve and edge 1.
		double alpha =  -pa.angle(drawnAxis,edgeAxis);			//angle to rotate everyone
		double scaleRatio = edgeAxis._mag()/drawnAxis._mag();	//ratio of distance from start to finish of drawn traj to distance between edges - multiply all elements in drawn traj by this

		//displace to align with start
		destCurve = movePoints(dispToStart, cntlPts);
		//build displacement vector - flip across axis if flipped curve
		myVector[] dispVecAra = new myVector[numPoints];
		dispVecAra[0] = new myVector();
		for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
			dispVecAra[myPointItr] = pa.V(destCurve[0],destCurve[myPointItr]);
		}			
		if((flip) || flags[isFlipped]){
			myVector udAxis = pa.U(drawnAxis);
			myVector normPt, tanPt;
			for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
				tanPt = myVector._mult(udAxis, dispVecAra[myPointItr]._dot(udAxis));			//component in udAxis dir
				normPt = myVector._sub(dispVecAra[myPointItr],tanPt);
				normPt._mult(2);
				dispVecAra[myPointItr]._sub(normPt);
			}
			flags[isFlipped] = flip;
		}

		//displace every point to be scaled distance from start of curve equivalent to scale of edge distances to drawn curve
		for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
			destCurve[myPointItr].set(pa.P(destCurve[0],scaleRatio, dispVecAra[myPointItr]));//start point displaced by scaleRatio * myVectortor from start to original location of myPoint
		}
		//rotate every point around destCurve[0] by alpha
		for(int myPointItr = 1; myPointItr < numPoints ; ++myPointItr){
			destCurve[myPointItr] = pa.R(destCurve[myPointItr],alpha,this.c_bAra[0], this.c_tAra[0],destCurve[0]);
		}
		//pa.outStr2Scr("cntlPts size : " + cntlPts.length,true);
		cntlPts = destCurve;
		//pa.outStr2Scr("cntlPts size : " + cntlPts.length,true);
		return destCurve;
	}//
	
	//move points by the passed myVectortor 
	public cntlPt[] movePoints(myVector move, cntlPt[] _pts){for(int i =0; i<_pts.length; ++i){	_pts[i]._add(move);	}	return _pts;}
	
	//calculate the weight of each point by determining the distance from its two neighbors - radius is inversely proportional to weight
	public float calcCntlWeight(myPoint a, myPoint p, myPoint b){	return (float)(myPoint._dist(a,p) + myPoint._dist(p,b));}
	
	//override this for cntrl-point driven constructs
	public int findClosestPt(myPoint p, double[] d){	
		return findClosestPt(p, d,cntlPts);
	}

	//drag point represented by passed idx in passed array - either point or cntl point
	public void dragPicked(myVector disp, int idx) {dragPicked(disp, idx, cntlPts);}
	//drag all points by finding displacement to mouse for COV and moving all points by same dislacement
	public void dragAll(myVector disp) {dragAll(disp, cntlPts);}	

	
	public String toString(){
		String res = "Interpolating Spline Stroke : Offset calc : "+ _offset;
		res += super.toString();
		return res;	
	}
	
}//myVariStroke
class cntlPt extends myPoint {
	public static MPM_SimMain pa;
	public int ID;
	public static int IDincr = 0;
	public static final float maxR = 75, 
			minR = 1,
			baseRad = 20;			//default radius for control points
	public float r, w;				//weight is calculated based on the distance to neighboring cntl myPoints when cntl myPoints are drawn
	public static int[][] clr = new int[][]{{0,0,255,255}, {111,111,111,255}};
	
	public cntlPt(MPM_SimMain _pa, myPoint _p, float _r, float _w){ super(_p.x,_p.y, _p.z); pa=_pa; ID=IDincr++;r=_r; w=_w; }
	public cntlPt(MPM_SimMain _pa, myPoint _p, float _w){this(_pa, _p, baseRad, _w);}
	public cntlPt(MPM_SimMain _pa, myPoint _p){this(_pa, _p, baseRad, baseRad);}
	public cntlPt(MPM_SimMain _pa){this(_pa, _pa.P(),1);}
	public cntlPt(cntlPt _p){this(cntlPt.pa, cntlPt.pa.P(_p),_p.w); r = _p.r; w = _p.w;ID = _p.ID;}		
	public static cntlPt L(cntlPt A, float s, cntlPt B){	return new cntlPt(cntlPt.pa, cntlPt.pa.L((myPoint)A, s, (myPoint)B), capInterpR(A.r, s, B.r), (1-s)*A.w + (s)*B.w);}//(1-s)*A.r + (s)*B.r,
	public static cntlPt P(cntlPt A, cntlPt B){	float s = .5f;return L(A, s, B);}
	public myPoint set(myPoint P){super.set(P); return (myPoint)this;}
	private static float capInterpR(float a, float s, float b){ float res = (1-s)*a + (s)*b; res = (res < minR ? minR : res > maxR ? maxR : res); return res;}
	public void drawMe(int cIdx, boolean flat){	pa.fill(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);  pa.stroke(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);		pa.show(this,2,-1,-1, flat);}		
	public void drawRad(int cIdx,myVector I, myVector J){
        pa.fill(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);  
        pa.stroke(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]); 
        pa.circle(this, r, I,J,20);
    }
	public void drawRad(myVector I, myVector J){
        pa.circle(this, r, I,J,20);
    }
	public void drawBall(int cIdx,myVector I, myVector J) {
	    float rhalf = this.r*0.5f;
	    myPoint center1 = pa.P(this);center1._add(pa.V(I)._mult(rhalf));
	    myPoint center2 = pa.P(this);center2._add(pa.V(I)._mult(-rhalf));
	    pa.fill(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);  
	    pa.stroke(0,0,0,255); 
        pa.circle(center1, rhalf, I,J,20);
        pa.circle(center2, rhalf, I,J,20);
        pa.show(center1,1);
        pa.show(center2,1);
    }
	public void drawNorm(int cIdx,myVector I, myVector J) {
	    myPoint p1 = pa.P(this);p1._add(pa.V(I)._mult(r));
        myPoint p2 = pa.P(this);p2._add(pa.V(I)._mult(-r));
        pa.stroke(0,0,0,255);
        pa.line(p1, p2); 
	}
	public void calcRadFromWeight(float lenRatio, boolean inv){r = Math.min(maxR, Math.max(minR, baseRad * (inv ?  (lenRatio/w) : (pa.wScale*w/(lenRatio*lenRatio)))));  }
	public void modRad(float modAmt){float oldR = r; r += modAmt; r = (r < minR ? minR : r > maxR ? maxR : r); w *= oldR/r; }
	public String toString(){String res = "Cntl Pt ID:"+ID+" p:"+super.toString()+" r:"+r+" w:"+w;return res;}
}//class cntlPoint
//class to hold functionality to calculate offset "sidewalks"
//3 types - normal, ball and radial, where normal is normal to 
//stroke line, radial is normal to resultant curve (by centering
//the ball on the center line) and ball is normal to both, by
//centering the ball at a particular radius away from the stroke
//line.

abstract class myOffset {
	public static MPM_SimMain pa;
	public int ID;
	public static int IDcnt = 0;
	public String name;
	public int capSize = 20;
	public boolean endCaps;

	public myOffset(MPM_SimMain _pa, boolean _ec){
		pa = _pa;
		ID = IDcnt++;
		endCaps = _ec;
	}	
	
	public myOffset(MPM_SimMain _pa){this(_pa, true);}
	
	/**
	 * calculate the offset points for the drawn stroke line contained in _obj
	 * @param _obj drawn stroke to build offset myPoints from
	 */
	public abstract ArrayList<myPoint> calcOffset(cntlPt[] cntlPts, myVector[] nAra, myVector[] tAra);				
	public abstract void drawCntlPts(cntlPt[] myPoints, myVector[] nAra, myVector[] tAra, boolean derived);
	
	/**
	 * build an array of points that sweeps around c clockwise in plane of norm and tan, with starting radius c.r * norm
	 * @param c control point
	 * @param norm normal of point (binormal in world frame, this is direction of offset)
	 * @param tan tangent at point
	 * @return myPoint array of sequence of points in an arc for an endcap
	 */
	public ArrayList<myPoint> buildCapPts (cntlPt c, myVector norm, myVector tan, float mult){
		ArrayList<myPoint> tmp = new ArrayList<myPoint>();
		float angle = PConstants.PI/(1.0f*capSize), sliceA = angle;			//10 slices
		tmp.add(pa.P((myPoint)c, mult * -c.r, norm));
		for(int i=1;i<capSize-1;++i){	tmp.add(pa.R(tmp.get(i-1), sliceA, norm, tan, c));}
		tmp.add(pa.P((myPoint)c, mult * c.r, norm));
		return tmp;	
	}//buildCapPts
	
	public String toString(){
		String res = name + "Offset ID : "+ID;		
		return res;
	}
}//myOffset

/**
* calculates normal offset - distance r, normal from stroke line
* @author john
*
*/
//make other classes to use different offset mechanism
class myNormOffset extends myOffset{
	myNormOffset(MPM_SimMain _pa){super(_pa); name = "Normal offset";}

	@Override
	public  ArrayList<myPoint> calcOffset(cntlPt[] cntlPts, myVector[] nAra, myVector[] tAra) {
		if(nAra.length != cntlPts.length){return  new ArrayList<myPoint>();}	
		ArrayList<myPoint> tmp = new ArrayList<myPoint>();
		int numCmyPointsM1 = cntlPts.length-1;		
		//start at first point and build endcap
		if(endCaps){tmp.addAll(buildCapPts(cntlPts[0], nAra[0], tAra[0], 1));}
		for(int i = 0; i<cntlPts.length;++i){	tmp.add(pa.P((myPoint)cntlPts[i], cntlPts[i].r, nAra[i]));}//add cntl point + rad offset from norm
		//build endcap on last cntlpoint
		if(endCaps){tmp.addAll(buildCapPts(cntlPts[numCmyPointsM1], nAra[numCmyPointsM1], tAra[numCmyPointsM1],-1));}
		for(int i = numCmyPointsM1; i>=0;--i){	tmp.add(pa.P((myPoint)cntlPts[i], -cntlPts[i].r, nAra[i]));}//add cntl point + rad offset from norm negated, in backwards order, so all points are added properly
		return tmp;
	}
	
	public String toString(){
		String res = name +super.toString();		
		return res;
	}
	
  @Override
//  public void drawCntlPts(cntlPt[] myPoints, myVector[] nAra, myVector[] tAra, boolean derived) {
//      for(int i = 0; i < myPoints.length; ++i){
//          myPoints[i].drawNorm((derived ? 0 : 1), nAra[i], tAra[i]);
//      }
//  }
  public void drawCntlPts(cntlPt[] myPoints, myVector[] nAra, myVector[] tAra, boolean derived) {
  	pa.pushStyle();
  	int clrInt = 0;
      for(int i = 0; i < myPoints.length; ++i){
      	clrInt = (int)(i/(1.0f * myPoints.length) * 255.0f);
          pa.fill(clrInt,0,(255 - clrInt),255);  
          pa.stroke(clrInt,0,(255 - clrInt),255); 
          myPoints[i].drawRad(nAra[i], tAra[i]);
      }
      pa.popStyle();
  }
      



}//myNormOffset

