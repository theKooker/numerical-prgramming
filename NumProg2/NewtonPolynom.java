import java.util.Arrays;
import java.util.OptionalDouble;

/**
 * Die Klasse Newton-Polynom beschreibt die Newton-Interpolation. Die Klasse
 * bietet Methoden zur Erstellung und Auswertung eines Newton-Polynoms, welches
 * uebergebene Stuetzpunkte interpoliert.
 *
 * @author braeckle
 *
 */
public class NewtonPolynom implements InterpolationMethod {

    /** Stuetzstellen xi */
    double[] x;

    /**
     * Koeffizienten/Gewichte des Newton Polynoms p(x) = a0 + a1*(x-x0) +
     * a2*(x-x0)*(x-x1)+...
     */
    double[] a;

    /**
     * die Diagonalen des Dreiecksschemas. Diese dividierten Differenzen werden
     * fuer die Erweiterung der Stuetzstellen benoetigt.
     */
    double[] f;

    /**
     * leerer Konstruktore
     */
    public NewtonPolynom() {
    };

    /**
     * Konstruktor
     *
     * @param x
     *          Stuetzstellen
     * @param y
     *          Stuetzwerte
     */
    public NewtonPolynom(double[] x, double[] y) {
        this.init(x, y);
    }

    /**
     * {@inheritDoc} Zusaetzlich werden die Koeffizienten fuer das
     * Newton-Polynom berechnet.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        x = new double[n + 1];
        double h = (b - a) / n;

        for (int i = 0; i < n + 1; i++) {
            x[i] = a + i * h;
        }
        computeCoefficients(y);
    }

    /**
     * Initialisierung der Newtoninterpolation mit beliebigen Stuetzstellen. Die
     * Faelle "x und y sind unterschiedlich lang" oder "eines der beiden Arrays
     * ist leer" werden nicht beachtet.
     *
     * @param x
     *          Stuetzstellen
     * @param y
     *          Stuetzwerte
     */
    public void init(double[] x, double[] y) {
        this.x = Arrays.copyOf(x, x.length);
        computeCoefficients(y);
    }

    /**
     * computeCoefficients belegt die Membervariablen a und f. Sie berechnet zu
     * uebergebenen Stuetzwerten y, mit Hilfe des Dreiecksschemas der
     * Newtoninterpolation, die Koeffizienten a_i des Newton-Polynoms. Die
     * Berechnung des Dreiecksschemas soll dabei lokal in nur einem Array der
     * Laenge n erfolgen (z.B. spaltenweise Berechnung). Am Ende steht die
     * Diagonale des Dreiecksschemas in der Membervariable f, also f[0],f[1],
     * ...,f[n] = [x0...x_n]f,[x1...x_n]f,...,[x_n]f. Diese koennen spaeter bei
     * der Erweiterung der Stuetzstellen verwendet werden.
     *
     * Es gilt immer: x und y sind gleich lang.
     */
    private void computeCoefficients(double[] y) {
        this.a = new double[y.length];
        this.f = new double[y.length];
        // erste Spalte
        for (int i = 0; i < y.length; i++) {
            f[i] = y[i];
        }
        a[0] = f[0];
        // die restlichen Spalten
        for (int i = 1; i < y.length; i++) {
            // Spaltenweise
            for (int j = 0; j < y.length - i; j++) {
                f[j] = (f[j + 1] - f[j]) / (x[i + j] - x[j]);
                if (j == 0) {
                    a[i] = f[j];
                }
            }
        }
    }

    /**
     * Gibt die Koeffizienten des Newton-Polynoms a zurueck
     */
    public double[] getCoefficients() {
        return a;
    }

    /**
     * Gibt die Dividierten Differenzen der Diagonalen des Dreiecksschemas f
     * zurueck
     */
    public double[] getDividedDifferences() {
        return f;
    }

    /**
     * addSamplintPoint fuegt einen weiteren Stuetzpunkt (x_new, y_new) zu x
     * hinzu. Daher werden die Membervariablen x, a und f vergoessert und
     * aktualisiert . Das gesamte Dreiecksschema muss dazu nicht neu aufgebaut
     * werden, da man den neuen Punkt unten anhaengen und das alte
     * Dreiecksschema erweitern kann. Fuer diese Erweiterungen ist nur die
     * Kenntnis der Stuetzstellen und der Diagonalen des Schemas, bzw. der
     * Koeffizienten noetig. Ist x_new schon als Stuetzstelle vorhanden, werden
     * die Stuetzstellen nicht erweitert.
     *
     * @param x_new
     *              neue Stuetzstelle
     * @param y_new
     *              neuer Stuetzwert
     */
    public void addSamplingPoint(double x_new, double y_new) {
       OptionalDouble optional =  Arrays.stream(x).filter(e -> e == x_new).findAny();
       if (optional.isPresent()) {
           //das neue x existiert schon und kann nicht wider hinzugefügt werden
           return;
       }
       //neu deklarieren der arrays.
       double f_new[] = new double[f.length + 1];
       double a_new[] = new double[a.length + 1];
       double xs_new[] = new double[x.length + 1];
       for(int i = 0; i < x.length; i++) {
           xs_new[i] = x[i];
       }
       xs_new[x.length] = x_new;
       f_new[f.length] = y_new;
       for(int i = f.length - 1; i >= 0; i--) {
        int k = f_new.length - 1 - i;
        f_new[i] = (f_new[i + 1] - f[i])/(xs_new[i + k] - xs_new[i]);
       } 

       for(int i = 0; i < this.a.length; i++) {
           a_new[i] = a[i];
       }
       a_new[a_new.length - 1] = f_new[f_new.length - 1];

       this.a = a_new;
       this.f = f_new;
       this.x = xs_new;
    }

    /**
     * {@inheritDoc} Das Newton-Polynom soll effizient mit einer Vorgehensweise
     * aehnlich dem Horner-Schema ausgewertet werden. Es wird davon ausgegangen,
     * dass die Stuetzstellen nicht leer sind.
     */
    @Override
    public double evaluate(double z) {
        double result = a[a.length - 1];
        for(int i = a.length - 2; i >= 0; i--) {
            result = this.a[i] + (z - this.x[i]) * result;
        }
        return result;
    }
}
