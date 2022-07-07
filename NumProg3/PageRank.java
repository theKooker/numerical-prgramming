import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

public class PageRank {

    /**
     * Diese Methode erstellt die Matrix A~ fuer das PageRank-Verfahren
     * PARAMETER:
     * L: die Linkmatrix: Wenn Verbindung von Website j zu i 1 sonst 0
     * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
     * zufaellig irgendeine Seite zu besuchen
     */
    public static double[][] buildProbabilityMatrix(int[][] L, double rho) {
        double[][] probabilityMatrix = new double[L.length][L.length];
        double links = 0;
        for (int i = 0; i < L.length; i++) {
            for (int j = 0; j < L.length; j++) {
                if (L[i][j] == 1) {
                    links = 0;
                    for (int k = 0; k < L.length; k++) {
                        links += L[k][j];
                    }
                    probabilityMatrix[i][j] = (1 - rho) * (1 / links) + (rho / L.length);
                } else {
                    probabilityMatrix[i][j] = (rho / L.length);
                }
            }
        }
        return probabilityMatrix;
    }

    /**
     * Diese Methode berechnet die PageRanks der einzelnen Seiten,
     * also das Gleichgewicht der Aufenthaltswahrscheinlichkeiten.
     * (Entspricht dem p-Strich aus der Angabe)
     * Die Ausgabe muss dazu noch normiert sein.
     * PARAMETER:
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, zufaellig irgendeine Seite zu besuchen
     * ,anstatt einem Link zu folgen.
     *
     */
    public static double[] rank(int[][] L, double rho) {
        double[][] A_tilde = buildProbabilityMatrix(L, rho);
        IntStream.range(0, A_tilde.length).forEach(i->A_tilde[i][i]-=1);
        double[] p = Gauss.solveSing(A_tilde);
        double lambda = 1/Arrays.stream(p).sum();
        double[] p_strich = Arrays.stream(p).map(x->x*lambda).toArray();
        return p_strich;
    }

    /**
     * Diese Methode erstellt eine Rangliste der uebergebenen URLs nach
     * absteigendem PageRank.
     * PARAMETER:
     * urls: Die URLs der betrachteten Seiten
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
     * zufaellig irgendeine Seite zu besuchen
     */
    public static String[] getSortedURLs(String[] urls, int[][] L, double rho) {
        int n = L.length;

        double[] p = rank(L, rho);

        RankPair[] sortedPairs = new RankPair[n];
        for (int i = 0; i < n; i++) {
            sortedPairs[i] = new RankPair(urls[i], p[i]);
        }

        Arrays.sort(sortedPairs, new Comparator<RankPair>() {

            @Override
            public int compare(RankPair o1, RankPair o2) {
                return -Double.compare(o1.pr, o2.pr);
            }
        });

        String[] sortedUrls = new String[n];
        for (int i = 0; i < n; i++) {
            sortedUrls[i] = sortedPairs[i].url;
        }

        return sortedUrls;
    }

    /**
     * Ein RankPair besteht aus einer URL und dem zugehoerigen Rang, und dient
     * als Hilfsklasse zum Sortieren der Urls
     */
    private static class RankPair {
        public String url;
        public double pr;

        public RankPair(String u, double p) {
            url = u;
            pr = p;
        }
    }
}
