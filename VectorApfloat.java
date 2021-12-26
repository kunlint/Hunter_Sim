/* *****************************************************************************
 *  Name:    Alan Turing
 *  NetID:   aturing
 *  Precept: P00
 *
 *  Description:  Prints 'Hello, World' to the terminal window.
 *                By tradition, this is everyone's first program.
 *                Prof. Brian Kernighan initiated this tradition in 1974.
 *
 **************************************************************************** */

import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;

public class VectorApfloat {
    private final int n;         // length of the VectorApfloat
    private Apfloat[] data;       // array of VectorApfloat's components

    // create the zero VectorApfloat of length n
    public VectorApfloat(int n) {
        this.n = n;
        this.data = new Apfloat[n];
    }

    // create a VectorApfloat from an array
    public VectorApfloat(Apfloat[] data) {
        n = data.length;

        // defensive copy so that client can't alter our copy of data[]
        this.data = new Apfloat[n];
        for (int i = 0; i < n; i++)
            this.data[i] = data[i];
    }

    // return the length of the VectorApfloat
    public int length() {
        return n;
    }

    // return the inner product of this VectorApfloat a and b
    public Apfloat dot(VectorApfloat that) {
        if (this.length() != that.length())
            throw new IllegalArgumentException("dimensions disagree");
        Apfloat sum = new Apfloat(0);
        for (int i = 0; i < n; i++)
            sum = sum.add(this.data[i].multiply(that.data[i]));
        return sum;
    }

    // return the Euclidean norm of this VectorApfloat
    public Apfloat magnitude() {
        return ApfloatMath.sqrt(this.dot(this));
    }

    // return the Euclidean distance between this and that
    public Apfloat distanceTo(VectorApfloat that) {
        if (this.length() != that.length())
            throw new IllegalArgumentException("dimensions disagree");
        return this.minus(that).magnitude();
    }

    // return this + that
    public VectorApfloat plus(VectorApfloat that) {
        if (this.length() != that.length())
            throw new IllegalArgumentException("dimensions disagree");
        VectorApfloat c = new VectorApfloat(n);
        for (int i = 0; i < n; i++)
            c.data[i] = this.data[i].add(that.data[i]);
        return c;
    }

    // return this - that
    public VectorApfloat minus(VectorApfloat that) {
        if (this.length() != that.length())
            throw new IllegalArgumentException("dimensions disagree");
        VectorApfloat c = new VectorApfloat(n);
        for (int i = 0; i < n; i++)
            c.data[i] = this.data[i].subtract(that.data[i]);
        return c;
    }

    // return the corresponding coordinate
    public Apfloat cartesian(int i) {
        return data[i];
    }

    // create and return a new object whose value is (this * factor)
    @Deprecated
    public VectorApfloat times(Apfloat factor) {
        VectorApfloat c = new VectorApfloat(n);
        for (int i = 0; i < n; i++)
            c.data[i] = factor.multiply(data[i]);
        return c;
    }

    // create and return a new object whose value is (this * factor)
    public VectorApfloat scale(Apfloat factor) {
        VectorApfloat c = new VectorApfloat(n);
        for (int i = 0; i < n; i++)
            c.data[i] = factor.multiply(data[i]);
        return c;
    }


    // return a string representation of the VectorApfloat
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append('(');
        for (int i = 0; i < n; i++) {
            s.append(data[i]);
            if (i < n - 1) s.append(", ");
        }
        s.append(')');
        return s.toString();
    }

    public static void main(String[] args) {
        StdOut.print('a');

    }
}
