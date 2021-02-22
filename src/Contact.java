package comp559.lcp;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

/**
 * Implementation of a contact constraint.
 * @author kry
 */
public class Contact {

    /** Next available contact index, used for determining which rows of the jacobian a contact uses */
    static public int nextContactIndex = 0;
    
    /** Index of this contact, determines its rows in the jacobian */
    int index;
    
    /** Store whether to iterate or not for warm start*/
    boolean iter = false;
    /** First RigidBody in contact */
    RigidBody body1;
    
    /** Second RigidBody in contact */
    RigidBody body2;
    
    /** Contact normal in world coordinates */
    Vector2d normal = new Vector2d();
    
    /** Position of contact point in world coordinates */
    Point2d contactW = new Point2d();
    
    Vector2d tangent = new Vector2d();
    double[] J1 = new double[6];
    double[] J2 = new double[6];
    
    /**
     * Creates a new contact, and assigns it an index
     * @param body1
     * @param body2
     * @param contactW
     * @param normal
     */
    public Contact( RigidBody body1, RigidBody body2, Point2d contactW, Vector2d normal ) {
        this.body1 = body1;
        this.body2 = body2;
        this.contactW.set( contactW );
        this.normal.set( normal );        
        index = nextContactIndex++;        
        // objective 3 TODO: you may want to add code here to compute and store the contact Jacobian
//        tangent perp norm
        this.tangent = new Vector2d(-normal.y, normal.x);
        Point2d ra = new Point2d(body1.x);
        ra.negate();
        ra.add(this.contactW);
        Point2d rb = new Point2d(body2.x);
        rb.negate();
        rb.add(this.contactW);
        Vector2d raPerp = new Vector2d (-ra.y, ra.x);
        Vector2d rbPerp = new Vector2d (-rb.y, rb.x);
        double tmp1 = raPerp.dot(this.normal);
        double tmp2 = rbPerp.dot(this.normal);
      //Jrow1 = [-nx,-ny, -n.'*ra_perp, nx,ny,n.'*rb_perp]
        this.J1 = new double[] { -this.normal.x, -this.normal.y, -tmp1, this.normal.x, this.normal.y,tmp2};
        tmp1 = raPerp.dot(tangent);
        tmp2 = rbPerp.dot(tangent);
        //Jrow2 = [-tx,-ty, -t.'*ra_perp, tx,ty,t.'*rb_perp]
        this.J2 = new double[] {-this.tangent.x, -this.tangent.y, -tmp1, this.tangent.x, this.tangent.y, tmp2};
        
    }
    
    /**
     * Draws the contact points
     * @param drawable
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        gl.glPointSize(3);
        gl.glColor3f(.7f,0,0);
        gl.glBegin( GL.GL_POINTS );
        gl.glVertex2d(contactW.x, contactW.y);
        gl.glEnd();
    }
    
    /**
     * Draws the connections between bodies to visualize the 
     * the adjacency structure of the matrix as a graph.
     * @param drawable
     */
    public void displayConnection( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        // draw a line between the two bodies but only if they're both not pinned
        if ( !body1.pinned && ! body2.pinned ) {
            gl.glLineWidth(2);
            gl.glColor4f(0,.3f,0, 0.5f);
            gl.glBegin( GL.GL_LINES );
            gl.glVertex2d(body1.x.x, body1.x.y);
            gl.glVertex2d(body2.x.x, body2.x.y);
            gl.glEnd();
        }
    }
    
}
