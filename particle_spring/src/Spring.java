package comp559.particle;

import javax.vecmath.Vector2d;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.DenseMatrix;
import javax.vecmath.Vector2d;

/**
 * Spring class for 599 assignment 1
 * @author kry
 */
public class Spring {

    Particle p1 = null;
    Particle p2 = null;
    
    /** Spring stiffness, sometimes written k_s in equations */
    public static double k = 1;
    /** Spring damping (along spring direction), sometimes written k_d in equations */
    public static double c = 20;
    /** Rest length of this spring */
    double l0 = 0;
    
    /**
     * Creates a spring between two particles
     * @param p1
     * @param p2
     */
    public Spring( Particle p1, Particle p2 ) {
        this.p1 = p1;
        this.p2 = p2;
        recomputeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
    }
    
    /**
     * Computes and sets the rest length based on the original position of the two particles 
     */
    public void recomputeRestLength() {
        l0 = p1.p0.distance( p2.p0 );
    }
    
    /**
     * Applies the spring force by adding a force to each particle
     */
    public void apply() {
        // TODO: Objective 1, FINISH THIS CODE!
    	double distNorm = p1.p.distance(p2.p);
    	Vector2d fa = new Vector2d();
    	fa.sub(p1.p,p2.p);
    	fa.normalize();
    	fa.scale(-k*(distNorm-l0));
    	p1.addForce(fa);
    	fa.negate();
    	p2.addForce(fa);
    	//spring damp
	    Vector2d sDamp = new Vector2d();
    	sDamp.sub(p1.p,p2.p);
    	sDamp.normalize();
    	Vector2d vdiff = new Vector2d();
    	vdiff.sub(p1.v,p2.v);
    	Vector2d xdiff = new Vector2d();
    	xdiff.sub(p1.p,p2.p);
    	double scale = vdiff.dot(xdiff);
    	scale*=distNorm;
    	scale*=-c;
    	sDamp.scale(scale);
    	p1.addForce(sDamp);
    	sDamp.negate();
    	p2.addForce(sDamp);
    }
   
    /** TODO: the functions below are for the backwards Euler solver */
    
    /**
     * Computes the force and adds it to the appropriate components of the force vector.
     * (This function is something you might use for a backward Euler integrator)
     * @param f
     */
    public void addForce( Vector f) {
        // TODO: Objective 8, FINISH THIS CODE for backward Euler method (probably very simlar to what you did above)
    	double distNorm = p1.p.distance(p2.p);
    	Vector2d fa = new Vector2d();
    	fa.sub(p1.p,p2.p);
    	fa.normalize();
    	fa.scale(-k*(distNorm-l0));
    	f.add(p1.index*2,fa.x);
    	f.add(p1.index*2+1,fa.y);
    	fa.negate();
    	f.add(p2.index*2,fa.x);
    	f.add(p2.index*2+1,fa.y);
	    Vector2d sDamp = new Vector2d();
    	sDamp.sub(p1.p,p2.p);
    	sDamp.normalize();
    	Vector2d vdiff = new Vector2d();
    	vdiff.sub(p1.v,p2.v);
    	Vector2d xdiff = new Vector2d();
    	xdiff.sub(p1.p,p2.p);
    	double scale = vdiff.dot(xdiff);
    	scale*=distNorm;
    	scale*=-c;
    	sDamp.scale(scale);
    	f.add(p1.index*2,sDamp.x);
    	f.add(p1.index*2+1,sDamp.y);
    	sDamp.negate();
    	f.add(p2.index*2,sDamp.x);
    	f.add(p2.index*2+1,sDamp.y);
    	
    }
    
    /**
     * Adds this springs contribution to the stiffness matrix
     * @param dfdx
     */
    public void addDfdx( Matrix dfdx ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward euler integration
        DenseVector tmpV = new DenseVector(2);
        DenseMatrix tmpM = new DenseMatrix (2,2);
        DenseMatrix tmp = new DenseMatrix(2,2);
        DenseMatrix I = new DenseMatrix (2,2);
        Vector2d l  = new Vector2d(p1.p.x-p2.p.x, p1.p.y-p2.p.y);//A-B
        tmpV.set(0,l.x/l.length());
        tmpV.set(1,l.y/l.length());//tmpV = A-B/||A-B||
        tmp.set(0,0,l.x);
        tmp.set(1,0,l.y);
        tmpM.rank1(tmpV);//tmpVtmpV^T
        I.add(0,0,1.0);
        I.add(1,1,1.0);
        I.scale(1-l0/l.length());
        
        tmpM.mult(tmp,tmpM);
        tmpM.mult(I, tmpM);
        tmpM.scale(-k);
        for (int i = 0;i<2;i++) {
        	dfdx.add(p1.index*2, p1.index*2+i, tmpM.get(0,i));
        	dfdx.add(p1.index*2+1, p1.index*2, tmpM.get(1,i));
        	dfdx.add(p2.index*2, p2.index*2+i, tmpM.get(0,i));
        	dfdx.add(p2.index*2+1, p2.index*2, tmpM.get(1,i));
        }
        tmpM.scale(-1);
        for (int i = 0;i<2;i++) {
        	dfdx.add(p1.index*2, p2.index*2+i, tmpM.get(0,i));
        	dfdx.add(p1.index*2+1, p2.index*2, tmpM.get(1,i));
        	dfdx.add(p2.index*2, p1.index*2+i, tmpM.get(0,i));
        	dfdx.add(p2.index*2+1, p1.index*2, tmpM.get(1,i));
        }
    }   
 
    /**
     * Adds this springs damping contribution to the implicit damping matrix
     * @param dfdv
     */
    public void addDfdv( Matrix dfdv ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward Euler integration
        DenseVector tmpV = new DenseVector(3);
        DenseMatrix tmpM = new DenseMatrix (3,3);
        Vector2d I = new Vector2d (p1.p.x-p2.p.x,p1.p.y-p2.p.y);
        tmpV.set(0,I.x/I.length());
        tmpV.set(1,I.y/I.length());
        tmpM.rank1(tmpV);
        tmpM.scale(-c);
        for (int i = 0;i<2;i++) {
        	dfdv.add(p1.index*2, p1.index*2+i, tmpM.get(0,i));
        	dfdv.add(p1.index*2+1, p1.index*2, tmpM.get(1,i));
        	dfdv.add(p2.index*2, p2.index*2+i, tmpM.get(0,i));
        	dfdv.add(p2.index*2+1, p2.index*2, tmpM.get(1,i));
        }
        tmpM.scale(-1);
        for (int i = 0;i<2;i++) {
        	dfdv.add(p1.index*2, p2.index*2+i, tmpM.get(0,i));
        	dfdv.add(p1.index*2+1, p2.index*2, tmpM.get(1,i));
        	dfdv.add(p2.index*2, p1.index*2+i, tmpM.get(0,i));
        	dfdv.add(p2.index*2+1, p1.index*2, tmpM.get(1,i));
        }
        
    } 
    
}
