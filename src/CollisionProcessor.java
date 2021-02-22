package comp559.lcp;

import java.util.ArrayList;
import java.util.List;
import java.lang.Math;
import java.util.Collections;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;

/**
 * Class for detecting and resolving collisions.  Currently this class uses penalty forces between rigid bodies.
 * @author kry
 */
public class CollisionProcessor {

    private List<RigidBody> bodies;
    
    /**
     * The current contacts that resulted in the last call to process collisions
     */
    public ArrayList<Contact> contacts = new ArrayList<Contact>();
    
    /**
     * Creates this collision processor with the provided set of bodies
     * @param bodies
     */
    public CollisionProcessor( List<RigidBody> bodies ) {
        this.bodies = bodies;
    }
    
    /** keeps track of the time used for collision detection on the last call */
    double collisionDetectTime = 0;
    
    /** keeps track of the time used to solve the LCP based velocity update on the last call */
    double collisionSolveTime = 0;
    
    /**
     * Processes all collisions 
     * @param dt time step
     */
    public void processCollisions( double dt ) {
        contacts.clear();
        Contact.nextContactIndex = 0;
        
        long now = System.nanoTime();
        broadPhase();
        collisionDetectTime = ( System.nanoTime() - now ) * 1e-9;
                
        if ( contacts.size() > 0  && doLCP.getValue() ) {
            now = System.nanoTime();

            double bounce = restitution.getValue();
            double mu = friction.getValue();
            // obj 3
            //TODO: Compute velocity update with iterative solve of contact constraint matrix.
            double[] lambda = new double[contacts.size()*2];
            double[] b = new double[lambda.length];
            double[] Dii = new double [lambda.length];
            double[] deltaV = new double [3*bodies.size()];
            for (Contact c : contacts) {
            	RigidBody b1 = c.body1;
            	RigidBody b2 = c.body2;
            	int i = c.index;
            	//b = J*(u + dt*minv*force)
            	double [] u = new double[] {b1.v.x, b1.v.y, b1.omega, b2.v.x,b2.v.y,b2.omega};
            	double [] f = new double[] {b1.force.x, b1.force.y, b1.torque, b2.force.x, b2.force.y, b2.torque};
            	double [] m = new double[] {b1.minv, b1.minv, b1.jinv, b2.minv, b2.minv, b2.jinv};
            	for (int j=0;j<u.length;j++) {
            		b[2*i]+=c.J1[j]*(u[j]+dt*f[j]*m[j]);
            		//bounce
            		b[2*i]+=c.J1[j]*bounce*u[j];
            		b[2*i+1]+=c.J2[j]*(u[j]+dt*f[j]*m[j]);
//            		Dii = Ji0^2*mAinv + ji1^2*mAinv + ji2^2*jAinv +...
            		Dii[2*i]+=c.J1[j]*c.J1[j]*m[j];
            		Dii[2*i+1]+=c.J2[j]*c.J2[j]*m[j];
            	}
            	if (warmStart.getValue()&&b1.BlockHash.containsValue(b2.index)) {
                	//if using warm start, initialize lambda with previously stored lambdas
        			lambda[2*i] = b1.lN;
        			lambda[2*i+1] = b1.lF;
        			c.iter = true;
        			for (int k=0;k<3;k++) {
	        			deltaV[3*b1.index+k] +=m[k]*c.J1[k]*b1.lambdaNorm;
	        			deltaV[3*b1.index+k] +=m[k]*c.J2[k]*b1.lambdaTan;
	        			deltaV[3*b2.index+k] +=m[k+3]*c.J1[k+3]*b2.lambdaNorm;
	        			deltaV[3*b2.index+k] +=m[k+3]*c.J2[k+3]*b2.lambdaTan;
	        		}
                }
            }
            for (int i = 0;i<iterations.getValue();i++) {
            	if (randomization.getValue())
            		Collections.shuffle(contacts);
	        	for (Contact c : contacts) {
	        		int j = c.index;
//	        	lambda_i_new = lambda_i + (-b-J*v)/Dii
	        		double lKplus1 = -b[2*j]/Dii[2*j] + lambda[2*j];
					for (int k = 0; k < 3; k ++) {
						lKplus1 -= c.J1[k] * deltaV[3* c.body1.index + k] / Dii[2* j];
						lKplus1 -= c.J1[k + 3]  * deltaV[3* c.body2.index + k]/ Dii[2* j];
					}
					lKplus1 = Math.max(0,lKplus1);
	        		double delLamb = lKplus1 - lambda[2*j];
	        		lambda[2*j] = lKplus1;
	        		double high = mu* lKplus1;
	        		
//	        		update deltaV = minv*J*delta_Lambda
	        		double [] m = new double[] {c.body1.minv, c.body1.minv, c.body1.jinv, c.body2.minv, c.body2.minv, c.body2.jinv};
//	        		if (!c.iter) {
	        		for (int k=0;k<3;k++) {
	        			deltaV[3*c.body1.index+k] +=m[k]*c.J1[k]*delLamb;
	        			deltaV[3*c.body2.index+k] +=m[k+3]*c.J1[k+3]*delLamb;
	        		}
//	        		}
	        		if(i==iterations.getValue()-1 && warmStart.getValue()) {
	        			//last iteration I update the lambda attributes
	        			c.body1.lN = lambda[2*j];
	        			c.body2.lN = lambda[2*j];
	        			c.body1.lambdaNorm = delLamb;
	        			c.body2.lambdaNorm = delLamb;
	        		}
	        		double lambFric = -b[2*j+1]/Dii[2*j+1] + lambda[2*j+1];
	        		for (int k = 0; k < 3; k ++) {
	        			lambFric -= c.J2[k] * deltaV[3*c.body1.index+k] / Dii[2*j+1];
	        			lambFric -= c.J2[k+3]  * deltaV[3*c.body2.index+k]/ Dii[2*j+1];
					}
//	        		projection
	        		lambFric = Math.max(lambFric, -high);
	        		lambFric = Math.min(lambFric, high);
	        		delLamb = lambFric - lambda[2*j+1];
	        		lambda[2*j+1] = lambFric;
//	        		if (!c.iter) {
	        		for (int k=0;k<3;k++) {
	        			deltaV[3*c.body1.index+k] +=m[k]*c.J2[k]*delLamb;
	        			deltaV[3*c.body2.index+k] +=m[k+3]*c.J2[k+3]*delLamb;
	        		}
//	        		}
	        		if(i==iterations.getValue()-1 && warmStart.getValue()) {
	        			//last iteration I update the lambda hash and position hash
	        			c.body1.lF = lambda[2*j+1];
	        			c.body1.lambdaTan = delLamb;
	        			c.body2.lF = lambda[2*j+1];
	        			c.body2.lambdaTan = delLamb;
	        		}
	        	}
            }
	        for (RigidBody bd : bodies) {
	        	int i = bd.index;
	        	bd.v.x += deltaV[3*i];
	        	bd.v.y += deltaV[3*i+1];
//	        	bd.omega += deltaV[3*i+2];
	        }
            collisionSolveTime = (System.nanoTime() - now) * 1e-9;
        }
    }
    
    /**
     * Checks for collisions between bodies.  Note that you can optionaly implement some broad
     * phase test such as spatial hashing to reduce the n squared body-body tests.
     * Currently this does the naive n squared collision check.
     */
    private void broadPhase() {
        // Naive n squared body test.. might not be that bad for small number of bodies 
        visitID++;
        for ( RigidBody b1 : bodies ) {
            for ( RigidBody b2 : bodies ) { // not so inefficient given the continue on the next line
                if ( b1.index >= b2.index ) continue;
                if ( b1.pinned && b2.pinned ) continue;                
                narrowPhase( b1, b2 );                
            }
        }        
    }
    
    /**
     * Checks for collision between boundary blocks on two rigid bodies.
     * TODO: This needs to be improved as the n-squared block test is too slow!
     * @param body1
     * @param body2
     */
    private void narrowPhase( RigidBody body1, RigidBody body2 ) {
        if ( ! useBVTree.getValue() ) {
            for ( Block b1 : body1.blocks ) {
                for ( Block b2 : body2.blocks ) {
                    processCollision( body1, b1, body2, b2 );
                }
            }
        } else {
            // object 2
        	//TODO: implement code to use hierarchical collision detection on body pairs
        	detection (body1,body2,body1.root,body2.root);
        }
    }
    //traverse recursively
    public void detection(RigidBody body1, RigidBody body2, BVNode bn1, BVNode bn2) {
    	Disc first = bn1.boundingDisc;    	
    	Disc second = bn2.boundingDisc;
    	if (bn1.visitID!=visitID) {
    		bn1.visitID = visitID;
    		bn1.boundingDisc.updatecW();
    	}
    	if (bn2.visitID != visitID){
    		bn2.visitID = visitID;
    		bn2.boundingDisc.updatecW();
    	}
    	if (first.intersects(second)) {
    		//If nodes are leaves we can process the collisions directly.
    		BVNode [] lst = {bn1.child1,bn1.child2,bn2.child1,bn2.child2};
    		if (bn1.isLeaf()&&bn2.isLeaf()) {		
    			processCollision(body1,bn1.leafBlock,body2,bn2.leafBlock);
    			if (warmStart.getValue() && (bn1.leafBlock !=null)) 
    				body1.BlockHash.put(bn1.leafBlock, body2.index); 
    			else if (warmStart.getValue() && (bn2.leafBlock !=null))
    				body2.BlockHash.put(bn2.leafBlock, body1.index);
    		} 
    		else if (bn1.isLeaf()) {
    			detection(body1,body2,bn1,lst[2]);
    			detection(body1,body2,bn1, lst[3]);
    		} 
    		else if (bn2.isLeaf()) {
    			detection(body1,body2,lst[0],bn2);
    			detection(body1,body2,lst[1],bn2);
    		} 
    		else {
    			detection(body1,body2,lst[0],lst[2]);
    			detection(body1,body2,lst[0],lst[3]);
    			detection(body1,body2,lst[1],lst[2]);
    			detection(body1,body2,lst[1],lst[3]);    			
    		}
    	}
    }	    
    /** 
     * The visitID is used to tag boundary volumes that are visited in 
     * a given time step.  Marking boundary volume nodes as visited during
     * a time step allows for a visualization of those used, but it can also
     * be used to more efficiently update the centeres of bounding volumes
     * (i.e., call a BVNode's updatecW method at most once on any given timestep)
     */
    int visitID = 0;
    
    /**
     * Resets the state of the collision processor by clearing all
     * currently identified contacts, and reseting the visitID for
     * tracking the bounding volumes used
     */
    public void reset() {
        contacts.clear();
        Contact.nextContactIndex = 0;
        visitID = 0;            
    }
    
    // some working variables for processing collisions
    private Point2d tmp1 = new Point2d();
    private Point2d tmp2 = new Point2d();
    private Point2d contactW = new Point2d();
    private Vector2d force = new Vector2d();
    private Vector2d contactV1 = new Vector2d();
    private Vector2d contactV2 = new Vector2d();
    private Vector2d relativeVelocity = new Vector2d();
    private Vector2d normal = new Vector2d();
        
    /**
     * Processes a collision between two bodies for two given blocks that are colliding.
     * Currently this implements a penalty force
     * @param body1
     * @param b1
     * @param body2
     * @param b2
     */
    private void processCollision( RigidBody body1, Block b1, RigidBody body2, Block b2 ) {        
        double k = contactSpringStiffness.getValue();
        double c1 = contactSpringDamping.getValue();
        double threshold = separationVelocityThreshold.getValue();
        boolean useSpring = enableContactSpring.getValue();
        boolean useDamping = enableContactDamping.getValue();
        
        body1.transformB2W.transform( b1.pB, tmp1 );
        body2.transformB2W.transform( b2.pB, tmp2 );
        double distance = tmp1.distance(tmp2);
        if ( distance < Block.radius * 2 ) {
            // contact point at halfway between points 
            // NOTE: this assumes that the two blocks have the same radius!
            contactW.interpolate( tmp1, tmp2, .5 );
            // contact normal
            normal.sub( tmp2, tmp1 );
            normal.normalize();
            // create the contact
            Contact contact = new Contact( body1, body2, contactW, normal);
            // simple option... add to contact list...
            contacts.add( contact );
            if ( ! doLCP.getValue() ) {
                // compute relative body velocity at contact point
                body1.getSpatialVelocity( contactW, contactV1 );
                body2.getSpatialVelocity( contactW, contactV2 );
                relativeVelocity.sub( contactV1, contactV2 );
                if ( -relativeVelocity.dot( normal ) < threshold ) {
                    if ( useSpring ) {
                        // spring force
                        double interpenetration = distance - Block.radius * 2; // a negative quantity
                        force.scale( -interpenetration * k, normal );
                        body2.applyContactForceW(contactW, force);
                        force.scale(-1);
                        body1.applyContactForceW(contactW, force);
                    }
                    if ( useDamping ) {
                        // spring damping forces!
                        // vertical
                        force.scale( relativeVelocity.dot(normal) * c1, normal );                    
                        body2.applyContactForceW( contactW, force );
                        force.scale(-1);
                        body1.applyContactForceW( contactW, force );
                    }
                }
            }
        }
    }
   
    /** Stiffness of the contact penalty spring */
    private DoubleParameter contactSpringStiffness = new DoubleParameter("penalty contact stiffness", 1e3, 1, 1e5 );
    
    /** Viscous damping coefficient for the contact penalty spring */
    private DoubleParameter contactSpringDamping = new DoubleParameter("penalty contact damping", 10, 1, 1e4 );
    
    /** Threshold for the relative velocity in the normal direction, for determining if spring force will be applied. */
    private DoubleParameter separationVelocityThreshold = new DoubleParameter( "penalty separation velocity threshold (controls bounce)", 1e-9, 1e-9, 1e3 );
    
    /** Enables the contact penalty spring */
    private BooleanParameter enableContactSpring = new BooleanParameter("enable penalty contact spring", true );
    
    /** Enables damping of the contact penalty spring */
    private BooleanParameter enableContactDamping = new BooleanParameter("enable penalty contact damping", true );
    
    /** Restitution parameter for contact constraints */
    public DoubleParameter restitution = new DoubleParameter( "restitution (bounce)", 0, 0, 1 );
    
    /** Coulomb friction coefficient for contact constraint */
    public DoubleParameter friction = new DoubleParameter("Coulomb friction", 0.33, 0, 2 );
    
    /** Number of iterations to use in projected Gauss Seidel solve */
    public IntParameter iterations = new IntParameter("iterations for GS solve", 10, 1, 500);
    
    /** Flag for switching between penalty based contact and contact constraints */
    private BooleanParameter doLCP = new BooleanParameter( "do LCP solve", false );
    
    /** Flag for enabling the use of hierarchical collision detection for body pairs */
    private BooleanParameter useBVTree = new BooleanParameter( "use BVTree", false );
    
    private BooleanParameter randomization = new BooleanParameter ("randomization",false);
    
    private BooleanParameter warmStart = new BooleanParameter ("warm start", false);
    /**
     * @return controls for the collision processor
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Collision Processing Controls") );
        vfp.add( useBVTree.getControls() );
        vfp.add( doLCP.getControls() );
        vfp.add(randomization.getControls());
        vfp.add(warmStart.getControls());
        vfp.add( iterations.getSliderControls() );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( friction.getSliderControls(false) );
        
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("penalty method controls") );
        vfp2.add( contactSpringStiffness.getSliderControls(true) );
        vfp2.add( contactSpringDamping.getSliderControls(true) );
        vfp2.add( separationVelocityThreshold.getSliderControls( true ) );
        vfp2.add( enableContactDamping.getControls() );
        vfp2.add( enableContactSpring.getControls() );
        
        CollapsiblePanel cp = new CollapsiblePanel(vfp2.getPanel());
        cp.collapse();
        vfp.add( cp );        
        return vfp.getPanel();
    }
    
}
