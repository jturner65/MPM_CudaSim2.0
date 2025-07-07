package MPM.BaseSim.sim;

import java.util.ArrayList;

import MPM.BaseSim.material.Base_MPMMaterial;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.windowUI.base.Base_DispWindow;
import base_Utils_Objects.io.messaging.MessageObject;

/**
 * Base class holding shared MPM simulation functionality regardless of solve method
 * @author John Turner
 *
 */
public abstract class Base_MPMSim {
    /**
     * Render interface to render results
     */
    public static IRenderInterface ri;
    /**
     * Owning window
     */
    protected Base_MPMSimWindow win;
    /**
     * Message object for logging
     */
    protected MessageObject msgObj;
    /**
     * Name of instancing sim
     */
    public final String simName;
    /**
     * Material quantities of particle matter
     * TODO : support multiple materials?
     */
    public Base_MPMMaterial mat;        
    /**
     * Current ui values describing variables used in the simulation
     */
    protected MPM_SimUpdateFromUIData currUIVals;
    
    public Base_MPMSimFlags simFlags; 
    
    ////////////////////////////////////////////////////
    //Sim instance variables, populated from currUIVals structure on creation/ui update
    
    /**
     * # of snowballs - all particles evenly distributed amongst this many snowballs
     */
    protected int numSnowballs;
    /**
     * # of particles total in the sim
     */
    protected int numParts;
    /**
     * sim iterations per frame
     */
    protected int simStepsPerFrame;    
    /**
     * # of cells per side - cube so same in all 3 dims
     * TODO : handle non-cube
     */
    protected int gridSideCount;
    /**
     * total # of grid cells in grid
     */
    protected int ttlGridCount;
    /**
     * Fraction of points to draw (every x'th point will be drawn)
     */
    protected int drawPointIncr;
    
    /**
     * timestep of simulation - 
     */
    protected float deltaT;
    /**
     * Initial particle velocities
     */
    protected float initVel;
    /**
     * mass of particles
     */
    protected float particleMass;
    /**
     * size of single dimension of grid cell
     */
    protected float cellSize;
    /**
     * Volume of a grid cell
     */
    protected float cellVolume;
    /**
     * simulation boundaries - symmetric cube, only need min and max, grid length per dim
     */
    //TODO handle non-cube/non-symmetric
    protected float minSimBnds, maxSimBnds, gridDim;
    /**
     * friction coefficients of wall
     */
    private float wallFric;
    /**
     * TODO: friction coefficients of colliders
     */
    private float collFric;
    /**
     * Scaling value to scale drawn vectors
     */
    protected float vecLengthScale;    
    /**
     * Minimum scale amount for visualization to fill cube frame in 3d world
     */
    protected float minSclAmt;
    /**
     * Scaling array based on dimensions of enclosing cube
     */
    protected float[] sclAmt;

    /**
     * const gravity vector for calculations
     */
    private final myVectorf gravity;
    
    /**
     * Holds precomputed gridlines for rendering grid
     * dispGridLines[0,1,2 for i-j, i-k, j-k][first idx][second idx][0: start, 1: end]
     */
    protected ArrayList<ArrayList<ArrayList<myVectorf[]>>> dispGridLines;
    /**
     * Every gridIncr'th grid line should be drawn
     */
    private final int gridIncr = 1;
    /**
     * Constructor
     * @param _pa
     * @param _win
     * @param _simName
     * @param _gravity
     * @param _currUIVals
     */
    public Base_MPMSim(IRenderInterface _pa, Base_MPMSimWindow _win, String _simName, float[] _gravity, MPM_SimUpdateFromUIData _currUIVals) {
        ri=_pa;win=_win;simName = _simName;    
        msgObj = win.getMsgObj();
        currUIVals = new MPM_SimUpdateFromUIData(win);
        gravity = new myVectorf(_gravity[0],_gravity[1],_gravity[2]);
        //mat's quantities are managed by UI - only need to instance once
        mat = new Base_MPMMaterial(_currUIVals);
        //setup flag array
        simFlags = buildSimFlags();    
    }//ctor
    
    /**
     * Instancing class simulation flags structure
     * @return
     */
    protected abstract Base_MPMSimFlags buildSimFlags();
    /**
     * Reset simulation environment
     * @param rebuildSim Should entire sim environment be rebuilt?
     */
    public final void resetSim(SimResetProcess rebuildSim) {
        if (rebuildSim == SimResetProcess.DoNothing) {
            return;
        }
        //stop simulation and reset
        Base_DispWindow.AppMgr.setSimIsRunning(false);
        simFlags.setSimIsBuilt(false);    
        //if only simulation/kernel parameters, don't rebuild simulation environment
        if (rebuildSim != SimResetProcess.RemakeKernel) {
            buildSimEnvironment(rebuildSim);
        }       
        //instance-specific reset - rebuild simulation kernel/solver
        resetSim_Indiv(rebuildSim);

        simFlags.setSimIsBuilt(true);
        msgObj.dispDebugMessage("Base_MPMSim:"+simName, "resetSim","Finished resetting/rebuilding sim");
    }//resetSim
    
    /**
     * Instance-specific reset code
     * @param rebuildSim
     */
    protected abstract void resetSim_Indiv(SimResetProcess rebuildSim);
    
    /**
     * Called by simFlags structure, when debug is set or cleared
     * @param val
     */
    public final void handleSimFlagsDebug(boolean val) {    
        handleSimFlagsDebug_Indiv(val);
    }
    
    /**
     * Instance class-specific debug handling
     * Called by simFlags structure, when debug is set or cleared
     * @param val
     */
    public abstract void handleSimFlagsDebug_Indiv(boolean val);
    
    /**
     * Functionality to rebuild the simulation environment
     * @param rebuildSim
     */
    protected final void buildSimEnvironment(SimResetProcess rebuildSim) {
        //initialize all particle values
        //Build the particle layout for this simulation
        initPartArrays();     
        
        // build initial grid vis structure
        buildGridVis();
        //determine sim-specific particle layouts
        if (rebuildSim == SimResetProcess.RebuildSim) {
            //call this if we wish to rebuild simulation layout
            buildPartLayouts();
        } else {
            //call this if we want to reinitialize existing simulation configuration
            reinitSimObjects();    
        } 

    }//buildSimEnvironment
    /**
     * Initialize environmental layout particle holders/arrays
     */
    protected abstract void initPartArrays();

    /**
     * build initial layout for particles for this simulation
     * @param partVals [OUT] map of particle locs, initial velocities and min/max vals being constructed
     */
    protected abstract void buildPartLayouts();    
    
    /**
     * Reinitialize existing sim - will resynthesize sample points but will not rederive locations of sim objs.
     * @param partVals
     */
    protected abstract void reinitSimObjects();    
    
    /**
     * Update the simulator values whenever UI values change
     * @param upd
     */
    public final void updateSimVals_FromUI(MPM_SimUpdateFromUIData upd) {
        
        //Determine level of reset required
        //Some UI changes should not require any changes to sim, others might require only remaking the kernels
        SimResetProcess procToDo = checkValuesForChanges(upd);
        ///////////////////////
        // Require complete rebuild
        // SimResetProcess.RebuildSim
        
        //# of grid cells per side of cube
        int oldGridSideCount = gridSideCount;
        gridSideCount = upd.getGridCellsPerSide();
        if (oldGridSideCount != gridSideCount) {
            buildGridVis();
        }
        //# of snowballs to make
        numSnowballs = upd.getNumSnowballs();
        //cell size
        cellSize = upd.getGridCellSize();    
        //cell volume
        cellVolume = cellSize*cellSize*cellSize;
        
        ////////////////////////////
        // Require remaking objects, but can retain # of objects, centers and base velocities
        // SimResetProcess.ResetSim
        //# of particles requested to have in simulator
        numParts = upd.getNumParticles();
        //Initial velocity of spheres
        initVel = upd.getInitVel();
                
        ////////////////////////////
        // Requires kernels to be remade but no changes to sim environment
        // SimResetProcess.RemakeKernel
        //time step
        deltaT = upd.getTimeStep();        
        //particle mass to use
        particleMass = upd.getPartMass();    
        // wall friction
        wallFric = upd.getWallFricCoeff();
        // non-wall collider friction : TODO not yet supported
        collFric = upd.getCollFricCoeff();
        // update materials to match ui values
        mat.updateMatVals_FromUI(upd);                
                
        ////////////////////////////
        // Requires no special handling or sim modification
        // SimResetProcess.DoNothing
        //# of simulation steps to perform between renders
        simStepsPerFrame = upd.getSimStepsPerFrame();
        //every x'th point to draw - use higher values to draw fewer points
        drawPointIncr = upd.getDrawPtIncr();
        //drawn vector/scalar quantities scaling
        vecLengthScale = upd.getDrawnVecScale() * gridSideCount/100.0f;        
        
        ////////////////////////////
        // calculated derived values dependent on UI values
        maxSimBnds = (gridSideCount*cellSize)/2.0f;
        minSimBnds = -maxSimBnds;
        gridDim = maxSimBnds - minSimBnds;    
        //total grid size       
        ttlGridCount=gridSideCount*gridSideCount*gridSideCount;
        //scale amount to fill gridX x gridY x gridZ visualization cube        
        sclAmt = Base_DispWindow.AppMgr.getScaled3dGridDims(1.0f/(gridSideCount * cellSize));
        
        //The minimum of all scaled grid dimensions
        minSclAmt = MyMathUtils.min(sclAmt);
        //Sim specific values
        updateSimVals_FromUI_Indiv(upd);
        //copy UI data to local var - copy to local last so that values that have changed can be observed
        currUIVals.setAllVals(upd);
        
        //UI changes forced reset of the simulator
        resetSim(procToDo);    
        
    }//updateSimVals_FromUI
    
    /**
     * Set UI values that are specific to an instancing sim type (either CUDA or CPU-based)
     * @param upd
     */
    protected abstract void updateSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd); 
    
    /**
     * 
     * @param upd
     * @return
     */    
    protected SimResetProcess checkValuesForChanges(MPM_SimUpdateFromUIData upd) {
        boolean rebuildSim = upd.checkSimRebuild(currUIVals);
        if(rebuildSim) {
            msgObj.dispDebugMessage("Base_MPMSim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.RebuildSim");
            return SimResetProcess.RebuildSim;
        }
        boolean resetSim = upd.checkSimReset(currUIVals);
        if(resetSim) {
            msgObj.dispDebugMessage("Base_MPMSim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.ResetSim");
            return SimResetProcess.ResetSim;
        }
        boolean matsHaveChanged = upd.haveMaterialValsChanged(currUIVals);
        if(matsHaveChanged) {
            msgObj.dispDebugMessage("Base_MPMSim:"+simName, "checkValuesForChanges","Materials have changed; Specifying SimResetProcess.RemakeKernel");
            return SimResetProcess.RemakeKernel;            
        }    
        boolean remakeKernel = upd.checkSimKernelRebuilt(currUIVals);
        if(remakeKernel) {
            msgObj.dispDebugMessage("Base_MPMSim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.RemakeKernel");
            return SimResetProcess.RemakeKernel;
        }
        msgObj.dispDebugMessage("Base_MPMSim:"+simName, "checkValuesForChanges","Specifying SimResetProcess.DoNothing - nothing pertinent has changed.");
        return SimResetProcess.DoNothing;
    }
    
    /**
     * Intialize all particle-based values upon sim reset
     */
    protected abstract void initValues_Parts();
    
    /**
     * Intialize all grid-based values upon sim reset
     */
    protected abstract void initValues_Grid();

    //cube wall normals
    protected final myVectorf[] wallNorms = new myVectorf[] {
        new myVectorf( 1.0, 0.0, 0.0),new myVectorf(-1.0, 0.0, 0.0),
        new myVectorf( 0.0, 1.0, 0.0),new myVectorf( 0.0,-1.0, 0.0),
        new myVectorf( 0.0, 0.0, 1.0),new myVectorf( 0.0, 0.0,-1.0)
    };    
    
    //collision detection for wall collisions, returns idx in wallNorms of collision
    public int checkWallCollision(myVectorf pos) {
        if(pos.x<=minSimBnds) {        return 0;    } else if(pos.x>=maxSimBnds) {    return 1;    }            
        if(pos.y<=minSimBnds) {        return 2;    } else if(pos.y>=maxSimBnds) {    return 3;    }
        if(pos.z<=minSimBnds) {        return 4;    } else if(pos.z>=maxSimBnds) {    return 5;    }
        return -1;
    }//checkWallCollision
    
    /**
     * Get a random location within +/- bound to place sphere. bound should 
     * take into account desired sphere's radius, so as to not breach collider box.
     * @param bound furthest +/- value per axis to locate center of sphere.
     * @return a 3D center coord
     */
    protected myPointf getRandSphereCenter(float bound) {        
        myPointf ctr = new myPointf(
                MyMathUtils.randomDouble(-1,1), 
                MyMathUtils.randomDouble(-1,1),
                MyMathUtils.randomDouble(-1,1));
        ctr._mult(bound);
        return ctr;
    }//getRandSphereCenter
    
        
    /**
     * Color value to correspond to val location in span between min and max, with a minimum allowed value
     * @param val
     * @param min
     * @param max
     * @return
     */
    protected int[] getClrValInt(float[] val, float[] min, float[] max) {
        int[] res = new int[3];
        for (int i=0;i<3;++i) {
            float denom = (max[i]-min[i]);
            if(denom <=0) {res[i] = 255;}
            res[i] = (int) (55.0f + 200.0f * (val[i] - min[i])/denom);    
        }
        return res;
    }//getClrValInt
    
    /**
     * Grayscale value to correspond to val location in span between min and max, with a minimum allowed value
     * @param val
     * @param min
     * @param max
     * @return
     */
    protected int[] getGreyValInt(float[] val, float[] min, float[] max) {
        int sum = 0;
        for (int i=0;i<3;++i) {
            float denom = (max[i]-min[i]);
            if(denom <=0) {sum += 255;}
            sum += (int) (55.0f + 200.0f * (val[i] - min[i])/denom);    
        }
        int avg = 50 + sum/3;
        //Give a slight blue tint, for snow
        return new int[] {avg,avg+5,avg+11} ;
    }//getClrValInt
    
    /**
     * Execute simStepsPerFrame step of simulation
     * @param modAmtMillis is in milliseconds, counting # of ms since last sim call 
     * @return true when simulation run is complete - true turns run sim flag off
     */        
    public final boolean simMe(float modAmtMillis) {
        boolean finished = false;
         for (int i=0;i<simStepsPerFrame;++i) {             
             finished |= simMe_Indiv(modAmtMillis);
         }    
        
        finished |= simMePost_Indiv(modAmtMillis);
        return finished;
    }
    
    /**
     * Instance-specific per-sim cycle code
     * @param modAmtMillis
     * @return
     */
    protected abstract boolean simMe_Indiv(float modAmtMillis);
    
    /**
     * Instance-specific post-sim cyle code
     * @param modAmtMillis
     * @return
     */
    protected abstract boolean simMePost_Indiv(float modAmtMillis);
    
    
    /**
     * sim method to show execution time and debug information for each sim step
     * @param modAmtMillis
     * @return
     */
    public abstract boolean simMeDebug(float modAmtMillis);    //simMeDebug
    
    //////////////////////////////////////////
    // Draw routines
    
    /**
     * draw 1 frame of results    //animTimeMod is in seconds, counting # of seconds since last draw
     * @param animTimeMod
     */
    public final void drawMe(float animTimeMod) {
        if(!simFlags.getSimIsBuilt()) {return;}//if not built yet, don't try to draw anything
        //render all particles - TODO determine better rendering method
        //minimum magnitude to use for rendering cutoff
        float minMag = MyMathUtils.EPS_F/vecLengthScale;
        ri.pushMatState();
        //set stroke values and visual scale
            ri.setStrokeWt(2.0f/minSclAmt);
            ri.scale(sclAmt[0], sclAmt[1], sclAmt[2]);    
            
            //point-based rendering
            if(simFlags.getShowParticles()){    _drawParts(animTimeMod, simFlags.getShowLocColors());}            
            if(simFlags.getShowPartVels()){        _drawPartVel(animTimeMod, minMag, drawPointIncr);}
            
            //draw colliders, if exist
            if(simFlags.getShowCollider()){        _drawColliders(animTimeMod);}
            
            //Grid-based rendering
            if(simFlags.getShowGrid()) {        _drawGrid();}
            if(simFlags.getShowGridVel()) {        _drawGridVel(animTimeMod, minMag);}
            if(simFlags.getShowGridAccel()){    _drawGridAccel(animTimeMod, minMag);}
            if(simFlags.getShowGridMass()) {    _drawGridMass(animTimeMod, minMag);}
        ri.popMatState();
    }//drawMe
    
    /**
     * Draw instance class particles
     * @param animTimeMod
     * @param showLocColors
     */
    protected abstract void _drawParts(float animTimeMod, boolean showLocColors);
    
    /**
     * Draw instance class particle velocities
     * @param animTimeMod
     * @param minMag minimum magnitude per axis to draw vector
     * @param pincr
     */
    protected abstract void _drawPartVel(float animTimeMod, float minMag, int pincr);
    
    /**
     * draw internal-to-sim colliders, if they exist
     * @param animTimeMod
     */
    protected abstract void _drawColliders(float animTimeMod);
    
    /**
     * Draw instance class grid velocities - use _drawGridVec method
     * @param animTimeMod
     * @param minMag minimum magnitude per axis to draw vector
     */
    protected abstract void _drawGridVel(float animTimeMod, float minMag);

    /**
     * Draw instance class grid accelerations - use _drawGridVec method
     * @param animTimeMod
     * @param minMag minimum magnitude per axis to draw vector
     */
    protected abstract void _drawGridAccel(float animTimeMod, float minMag);
    
    /**
     * Draw instance class grid masses - use _drawGridScalar method
     * @param animTimeMod
     * @param minMag minimum magnitude to draw scalar mass
     */
    protected abstract void _drawGridMass(float animTimeMod, float minMag);    
    
    /**
     * Build a graphical representation of the computational grid, for rendering.
     */
    protected final void buildGridVis() {
        dispGridLines = new ArrayList<ArrayList<ArrayList<myVectorf[]>>>();        
        ArrayList<ArrayList<myVectorf[]>> ijGrid = new ArrayList<ArrayList<myVectorf[]>>(),
                ikGrid = new ArrayList<ArrayList<myVectorf[]>>(),
                jkGrid = new ArrayList<ArrayList<myVectorf[]>>();         

        for (int i=0; i<=gridSideCount;i+=gridIncr) {
            ArrayList<myVectorf[]> jLines = new ArrayList<myVectorf[]>(); 
            float iLoc = i*cellSize;
            for(int j=0;j<=gridSideCount;j+=gridIncr) {
                myVectorf startPos=new myVectorf(iLoc,j*cellSize,0.0f);
                jLines.add(new myVectorf[] {startPos, new myVectorf(iLoc, startPos.y,gridDim)});
            }
            ijGrid.add(jLines);
            ArrayList<myVectorf[]> kLines = new ArrayList<myVectorf[]>(); 
            for(int k=0;k<=gridSideCount;k+=gridIncr) {
                myVectorf startPos=new myVectorf(iLoc,0.0f, k*cellSize);
                kLines.add(new myVectorf[] {startPos, new myVectorf(iLoc,gridDim,startPos.z)});
            }
            ikGrid.add(kLines);
        }
        for(int j=0;j<=gridSideCount;j+=gridIncr) {
            float jLoc = j*cellSize;
            ArrayList<myVectorf[]> jkLines = new ArrayList<myVectorf[]>(); 
            for(int k=0;k<=gridSideCount;k+=gridIncr) {
                myVectorf startPos=new myVectorf(0.0f,jLoc,k*cellSize);
                jkLines.add(new myVectorf[] {startPos, new myVectorf(gridDim,jLoc,startPos.z)});
            }
            jkGrid.add(jkLines);
        }
        dispGridLines.add(ijGrid);
        dispGridLines.add(ikGrid);
        dispGridLines.add(jkGrid);
    }//buildGrid
    
    /**
     * Draw a representation of the eulerian grid
     */
    protected final void _drawGrid() {
        ri.pushMatState();        
            ri.setStroke(255,255,255,20);
            ri.translate(minSimBnds,minSimBnds,minSimBnds);
            for(var dispGridAra : dispGridLines) {
                for (var dispGridLinesAra : dispGridAra) {
                    for(myVectorf[] lineAra : dispGridLinesAra) {
                        ri.drawLine(lineAra[0],lineAra[1]);
                    }                    
                }
            }
        ri.popMatState();        
    }//_drawGrid()

    /**
     * Draw grid vector quantities
     * @param clr
     * @param val 2d array of grid values w/ first idx is x/y/z and 2nd is quantity
     * @param yVal y dim value to draw
     * @param zVal z dim value to draw
     * @param grid_pos 2d array of grid positions w/ first idx is x/y/z and 2nd is quantity
     * @param minMag minimum magnitude in any direction to display vector
     */
    protected final void _drawGridVec(int[] clr, float[][] val, float[][] grid_pos, float minMag) {
        ri.pushMatState();    
        ri.setStroke(clr,180);
        ri.translate(minSimBnds,minSimBnds,minSimBnds);
        for (int i=0; i<ttlGridCount;++i) {            
            if(        (Math.abs(val[0][i]) > minMag) || 
                    (Math.abs(val[1][i]) > minMag) || 
                    (Math.abs(val[2][i]) > minMag)) {
                ri.pushMatState();    
                ri.translate(grid_pos[0][i], grid_pos[1][i], grid_pos[2][i]);
                ri.drawLine(0,0,0, vecLengthScale*val[0][i],vecLengthScale*val[1][i],vecLengthScale*val[2][i]);
                ri.popMatState();
            }
        }
        ri.popMatState();
    }//_drawGridVec
    /**
     * Draw grid scalar quantities
     * @param clr
     * @param xVal value to draw
     * @param grid_pos 2d array of grid positions w/ first idx is x/y/z and 2nd is quantity
     * @param minMag minimum magnitude in any direction to display scalar quanatity
     */
    protected final void _drawGridScalar(int[] clr, float[] xVal, float[][] grid_pos, float minMag) {
        ri.pushMatState();    
        ri.setSphereDetail(3);
        ri.setStroke(clr,111);
        ri.translate(minSimBnds,minSimBnds,minSimBnds);
        for (int i=0; i<ttlGridCount;++i) {            
            if((Math.abs(xVal[i]) > minMag)) {
                ri.pushMatState();    
                ri.translate(grid_pos[0][i], grid_pos[1][i], grid_pos[2][i]);
                ri.drawSphere(xVal[i]*vecLengthScale);
                ri.popMatState();
            }
        }
        ri.popMatState();        
    }
    
    public final float getCellSize() {return cellSize;}
    /**
     * Number of grid cells per side of euler projection grid
     * @return
     */
    public final int getGridSideCount() {return gridSideCount;}
    /**
     * Minimum dimension of cube centered at origin serving as wall boundary
     * @return
     */
    public final myVectorf getMinSimBnds() {return new myVectorf(minSimBnds,minSimBnds,minSimBnds);}
    /**
     * Maximum dimension of cube centered at origin serving as wall boundary
     * @return
     */
    public final myVectorf getMaxSimBnds() {return new myVectorf(maxSimBnds,maxSimBnds,maxSimBnds);}
    /**
     * Sim Boundary wall friction
     * @return
     */
    public float getWallFric() {        return wallFric;    }

    /**
     * Current timestep
     * @return
     */
    public final float getDeltaT() {return deltaT;}
    
    /**
     * Constant gravity vector
     * @return
     */
    public myVectorf getGravity() {        return gravity;    }
    /**
     * Collider friction. TODO support multiple colliders
     * @return
     */
    public float getCollFric() {        return collFric;    }

    /**
     * Retrieve message object
     * @return
     */
    public MessageObject getMsgObj() {return msgObj;}
    
}//class Base_MPMSim
