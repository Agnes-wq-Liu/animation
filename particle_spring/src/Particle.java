package comp559.particle;

import java.util.ArrayList;

import javax.vecmath.Color3f;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

/**
 * Particle class that contains particle properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this particle.
 * @author kry
 */
public class Particle {
    
    /** Identifies this particles position in the particle list */
    int index;
    
    public boolean pinned = false;
    
    Color3f color = new Color3f(0,.95f,0);
    
    float size = 10;
    
    double mass = 1;
        
    Point2d p = new Point2d();
    
    Vector2d v = new Vector2d();
    
    Point2d p0 = new Point2d();
    
    Vector2d v0 = new Vector2d();
    
    Vector2d f = new Vector2d();

    /**
     * A list of springs that use this particle.  This list is only needed
     * to adjust rest lengths when dragging particles around.
     * This is only used for UI... it is probably not needed for simulation 
     */
    ArrayList<Spring> springs = new ArrayList<Spring>();
    
    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    public Particle( double x, double y, double vx, double vy ) {
        p0.set(x,y);
        v0.set(vx,vy);
        reset();
    }
    
    /**
     * Resets the position of this particle
     */
    public void reset() {
        p.set(p0);
        v.set(v0);
        f.set(0,0);
    }
    
    /**
     * Clears all forces acting on this particle
     */
    public void clearForce() {
        f.set(0,0);
    }
    
    /**
     * Adds the given force to this particle
     * @param force
     */
    public void addForce( Vector2d force ) {
        f.add(force);
    }
    
    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    public double distance( double x, double y ) {
        Point2d tmp = new Point2d( x, y );
        return tmp.distance(p);
    }
    
    @Override
    public String toString() {
        return f.toString();
    }
    
}
