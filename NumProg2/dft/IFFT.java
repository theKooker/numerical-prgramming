package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     *
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static Complex[] ifft(Complex[] c) {
        Complex[] v = new Complex[c.length];
        if(c.length == 1)
            v[0] = c[0];
        else {
            int m = c.length/2;
            Complex[] z1 = new Complex[m];
            Complex[] z2 = new Complex[m];
            for (int  i = 0; i < m; i++) {
                z1[i] = c[i*2];
                z2[i] = c[i*2+1];
            }
            z1 = ifft(z1);
            z2 = ifft(z2);
            Complex omega = Complex.fromPolar(1, 2*Math.PI/ c.length);
            for (int i = 0; i < m; i++) {
                v[i] = z1[i].add(omega.power(i).mul(z2[i]));
                v[m+i] = z1[i].sub(omega.power(i).mul(z2[i]));
            }
        }
        return v;
    }
}
