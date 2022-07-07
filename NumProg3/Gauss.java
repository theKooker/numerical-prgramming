import java.util.Arrays;

public class Gauss {
    private static class Tuple<U, V> {
        U matrix;
        V result;

        Tuple(U matrix, V result) {
            this.matrix = matrix;
            this.result = result;
        }

        public U getMatrix() {
            return matrix;
        }

        public V getResult() {
            return result;
        }
    }

    /**
     * Diese Methode soll die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution ermitteln.
     * PARAMETER:
     * R: Eine obere Dreiecksmatrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
        int n = b.length;
        double[] x = new double[n];
        x[n - 1] = b[n - 1] / R[n - 1][n - 1];

        for (int i = n - 2; i >= 0; i--) {
            x[i] = b[i];
            for (int j = n - 1; j > i; j--) {
                x[i] += -R[i][j] * x[j];
            }
            x[i] = x[i] / R[i][i];
        }

        return x;
    }

    private static Tuple<double[][], double[]> pivotOperationOnColumn(double[][] A, int column, double[] b) {
        // wir führen Selection-Sort duch
        for (int j = column; j < A.length - 1; j++) {
            if (Math.abs(A[j][column]) < Math.abs(A[j + 1][column])) {
                double[] temp = A[j];
                A[j] = A[j + 1];
                A[j + 1] = temp;
                double bTemp = b[j];
                b[j] = b[j + 1];
                b[j + 1] = bTemp;
            }
        }

        return new Tuple<double[][], double[]>(A, b);
    }

    /**
     * Diese Methode soll die Loesung x des LGS A*x=b durch Gauss-Elimination mit
     * Spaltenpivotisierung ermitteln. A und b sollen dabei nicht veraendert werden.
     * PARAMETER: A:
     * Eine regulaere Matrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] solve(double[][] A, double[] b) {
        for (int i = 0; i < A.length; i++) {
            var t = pivotOperationOnColumn(A, i, b);
            A = t.getMatrix();
            b = t.getResult();

            for (int j = i + 1; j < A.length; j++) {
                double scale = A[j][i] / A[i][i];
                if (A[j][i] != 0) {
                    for (int k = i; k < A.length; k++) {
                        A[j][k] -= A[i][k] * scale;
                    }
                    b[j] -= b[i] * (scale);
                }

            }

        }
        return backSubst(A, b);
    }

    /**
     * Diese Methode soll eine Loesung p!=0 des LGS A*p=0 ermitteln. A ist dabei
     * eine nicht invertierbare Matrix. A soll dabei nicht veraendert werden.
     *
     * Gehen Sie dazu folgendermassen vor (vgl.Aufgabenblatt):
     * -Fuehren Sie zunaechst den Gauss-Algorithmus mit Spaltenpivotisierung
     * solange durch, bis in einem Schritt alle moeglichen Pivotelemente
     * numerisch gleich 0 sind (d.h. <1E-10)
     * -Betrachten Sie die bis jetzt entstandene obere Dreiecksmatrix T und
     * loesen Sie Tx = -v durch Rueckwaertssubstitution
     * -Geben Sie den Vektor (x,1,0,...,0) zurueck
     *
     * Sollte A doch intvertierbar sein, kann immer ein Pivot-Element gefunden
     * werden(>=1E-10).
     * In diesem Fall soll der 0-Vektor zurueckgegeben werden.
     * PARAMETER:
     * A: Eine singulaere Matrix der Groesse n x n
     */
    public static double[] solveSing(double[][] A) {

        // Initialisierung
        int n = A.length;
        double[][] Ar = new double[n][n];// A soll nicht veraendert werden
        double[] b = new double[n];// b wird automatisch mit 0 gefüllt
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Ar[i][j] = A[i][j];
            }
        }
        // Gauss-Algorithmus mit Spaltenpivotisierung
        int newSize = 0;
        A: for (int i = 0; i < n; i++) {
            var t = pivotOperationOnColumn(Ar, i, b);
            Ar = t.getMatrix();
            b = t.getResult();

            for (int j = i + 1; j < n; j++) {

                double scale = Ar[j][i] / Ar[i][i];
                if (Ar[j][i] != 0) {
                    for (int k = i; k < n; k++) {
                        Ar[j][k] -= Ar[i][k] * scale;
                    }
                    b[j] -= b[i] * (scale);
                }
                if (Math.abs(Ar[j][j]) < 1E-10) {
                    newSize = i + 1;
                    break A;
                }

            }
        }
        if (newSize == 0) {
            b = new double[n];
            return b;
        }
        // Konstruktion der neuen Matrix
        // T ist bereits eine Dreiecksmatrix
        double[][] T = new double[newSize][newSize];
        double[] v = new double[newSize];
        for (int i = 0; i < newSize; i++) {
            for (int j = i; j < newSize; j++) {
                T[i][j] = Ar[i][j];
            }
            v[i] = -Ar[i][newSize];
        }

        double[] x = backSubst(T, v);
        for (int i = 0; i < newSize; i++) {
            b[i] = x[i];
        }
        b[newSize] = 1;
        for (int i = newSize + 1; i < n; i++) {
            b[i] = 0;
        }

        return b;
    }

    /**
     * Diese Methode berechnet das Matrix-Vektor-Produkt A*x mit A einer nxm
     * Matrix und x einem Vektor der Laenge m. Sie eignet sich zum Testen der
     * Gauss-Loesung
     */
    public static double[] matrixVectorMult(double[][] A, double[] x) {
        int n = A.length;
        int m = x.length;

        double[] y = new double[n];

        for (int i = 0; i < n; i++) {
            y[i] = 0;
            for (int j = 0; j < m; j++) {
                y[i] += A[i][j] * x[j];
            }
        }

        return y;
    }

    // public static void main(String[] args) {
    //     // TEST 1
    //     System.out.println("#############TEST 1#############");
    //     double[][] A = { { 4, 2, 3 }, { 2, 2, 1 }, { 2, 2, 2 } };
    //     double[][] B = { { 2, -2, 1 }, { 1, 3, -2 }, { 3, -1, -1 } };
    //     double[][] C = { { 1, -2, 1 }, { 0, 1, -2 }, { 0, 1, -2 } };
    //     double[] r = { 5, -3, 0 };
    //     double[] r2 = { -3, 1, 2 };
    //     double[] secondResult = solve(A, r);
    //     System.out.println(Arrays.toString(secondResult));
    //     double[] thirdResult = solveSing(C);
    //     System.out.println(Arrays.toString(thirdResult));
    //     // double[] test = { 3, 2, 1 };
    //     // System.out.println(Arrays.toString(matrixVectorMult(C, test)));
    //     System.out.println("#############TEST 2#############");
    //     double[][] D = { { 1, 1, 1 }, { 2, 3, 3 }, { 3, 5, 8 } };
    //     double[] r3 = { -1, 1, 0 };
    //     System.out.println(Arrays.toString(solve(D, r3)));
    //     System.out.println("#############TEST 3#############");
    //     double[][] E = { {-0.001,1},  {2,1}};
    //     double[] r4 = { 1, 0 };
    //     System.out.println(Arrays.toString(solve(E, r4)));

    // }
}
