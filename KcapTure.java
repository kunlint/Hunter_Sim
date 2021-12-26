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
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

import java.awt.Color;
import java.awt.Font;
import java.awt.Paint;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

public class KcapTure {
    private Apfloat convertToeV = new Apfloat("931494320.0", 500);
    private Apfloat csGround = new Apfloat("130.9054649", 500);
    // in amu
    private Apfloat xeGround = new Apfloat("130.90508406", 500);
    // in amu
    private Apfloat electron = new Apfloat("511000.0", 500);  // electron mass in eV
    private Apfloat m1; // after first decay
    private Apfloat m2; // after second decay
    private Apfloat m3; // after third decay
    private Apfloat m4; // final stable Xenon ion
    private Apfloat m11; // some intermediate state

    public KcapTure() {
        m4 = xeGround.multiply(convertToeV).subtract(electron).
                add(new Apfloat(12));
        m3 = m4.add(new Apfloat(23)).subtract(new Apfloat(4));
        m2 = xeGround.multiply(convertToeV).add(new Apfloat(146)).subtract(new Apfloat(4));
        m1 = m2.add(new Apfloat(34430));
        m11 = m2.add(new Apfloat(29600));
    }

    public Apfloat randomTheta1() {
        double a = StdRandom.uniform();
        a = a * 2 - 1;
        Apfloat theta = new Apfloat(a);
        return ApfloatMath.acos(theta);
    }

    public Apfloat randomTheta2() {
        double a = StdRandom.uniform() * 2 * Math.PI;
        Apfloat theta = new Apfloat(a);
        return theta;
    }

    public VectorApfloat formVector(Apfloat theta1, Apfloat theta2, Apfloat p) {
        Apfloat[] array = new Apfloat[3];
        array[0] = p.multiply(ApfloatMath.sin(theta1)).multiply(ApfloatMath.cos(theta2));
        array[1] = p.multiply(ApfloatMath.sin(theta1)).multiply(ApfloatMath.sin(theta2));
        array[2] = p.multiply(ApfloatMath.cos(theta1));
        VectorApfloat u = new VectorApfloat(array);
        return u;
    }

    // return the momentum of the ion in the ion rest frame
    private Apfloat twoBody(Apfloat M1, Apfloat M2, Apfloat m) {
        Apfloat p = ApfloatMath.pow(M1, 2)
                               .subtract(ApfloatMath.pow(M2, 2)).
                                       subtract(ApfloatMath.pow(m, 2));
        Apfloat x = ApfloatMath.pow(M2, 2).multiply(ApfloatMath.pow(m, 2))
                               .multiply(new Apfloat(4));
        p = ApfloatMath.sqrt(ApfloatMath.pow(p, 2).subtract(x));
        p = p.divide(new Apfloat(2)).divide(M1);
        return p;
    }

    // return the corresponding neutrino momentum for a given neutrino mass
    // nuetrino mass in keV
    public Apfloat twoBodynu(Apfloat m) {
        m = m.multiply(new Apfloat(1000));
        Apfloat M1 = csGround.multiply(convertToeV);
        Apfloat M2 = m1;
        Apfloat p = ApfloatMath.pow(M1, 2)
                               .subtract(ApfloatMath.pow(M2, 2)).
                                       subtract(ApfloatMath.pow(m, 2));
        Apfloat x = ApfloatMath.pow(M2, 2).multiply(ApfloatMath.pow(m, 2)).
                multiply(new Apfloat(4));
        p = ApfloatMath.sqrt(ApfloatMath.pow(p, 2).subtract(x));
        p = p.divide(new Apfloat(2)).divide(M1);
        return p;
    }

    // lorentz transform back to lab frame
    // return a 4-vector
    private Apfloat lorentz(VectorApfloat p, VectorApfloat b) {
        Apfloat[] array = new Apfloat[4];
        Apfloat beta = b.magnitude();
        Apfloat gamma = new Apfloat(1)
                .divide(ApfloatMath.sqrt(new Apfloat(1).subtract(beta.multiply(beta))));
        array[0] = p.cartesian(0).multiply(gamma)
                    .add(p.cartesian(1).multiply(gamma).multiply(b.cartesian(0)).add
                            (p.cartesian(2).multiply(gamma).multiply(b.cartesian(1)))
                          .add(p.cartesian(3).multiply(gamma).multiply(b.cartesian(2))));
        return array[0];
    }

    // return the reconstructed mass
    public VectorApfloat[] kcapturep(Apfloat nu) {
        nu = nu.multiply(new Apfloat(1000));
        Apfloat p = twoBody(csGround.multiply(convertToeV), m1, nu);
        // StdOut.println(p.multiply(p).add(nu.multiply(nu)));
        Apfloat theta1 = randomTheta1();
        Apfloat theta2 = randomTheta2();
        VectorApfloat pvector1 = formVector(theta1, theta2, p);
        // StdOut.println(pvector1.magnitude().multiply(pvector1.magnitude()));
        Apfloat pint = csGround.multiply(convertToeV)
                               .subtract(ApfloatMath.sqrt(m1.multiply(m1).add(p.multiply(p))));
        // end of first decay
        Apfloat theta3 = randomTheta1(); // propogate this one
        Apfloat theta4 = randomTheta2();
        Apfloat p1 = twoBody(m1, m2, new Apfloat(0, 500));
        VectorApfloat pvector2 = pvector1.scale(new Apfloat(1).divide(m1));
        VectorApfloat pvectorx1 = formVector(theta3, theta4, p1);
        VectorApfloat pvectorx = pvectorx1.scale(new Apfloat(1).divide(m2));
        pvector2 = pvector2.minus(pvectorx);
        pvector2 = pvector2.scale(m2);
        pvectorx1 = pvector1.minus(pvector2);
        Apfloat p2 = pvector2.magnitude();
        Apfloat pint2 = csGround.multiply(convertToeV)
                                .subtract(ApfloatMath.sqrt(m2.multiply(m2).add(p2.multiply(p2))))
                                .subtract(pvectorx1.magnitude());
        Apfloat p12 = pvectorx1.plus(pvector2).magnitude();
        //  end of second decay
        Apfloat pboost = twoBody(m2, m3, electron);
        // check energy conservation in the rest frame of the atom
        Apfloat theta7 = randomTheta1();
        Apfloat theta8 = randomTheta2();
        VectorApfloat pvector3copy2 = pvector2;
        VectorApfloat pvector3 = pvector2.scale(new Apfloat(1).divide(m2));
        VectorApfloat intere = pvector3.scale(electron);
        // electron momentum in Xe restframe
        VectorApfloat pvectore1 = formVector(theta7, theta8, pboost);
        // electron momentum in lab frame
        VectorApfloat pvectore1ab = pvectore1.plus(intere);
        VectorApfloat pvectore = pvectore1.scale(new Apfloat(1).divide(m3));
        pvector3 = pvector3.minus(pvectore);
        pvector3 = pvector3.scale(m3);
        pvectore1ab = pvector3copy2.minus(pvector3);
        // end of third decay
        Apfloat theta9 = randomTheta1();
        Apfloat theta10 = randomTheta2();
        VectorApfloat pvector3copy = pvector3;
        pvector3 = pvector3.scale(new Apfloat(1).divide(m3));
        Apfloat p3 = twoBody(m3, m4, new Apfloat(0, 500));
        VectorApfloat pvectoruv = formVector(theta9, theta10, p3);
        VectorApfloat pvectoruv2 = pvectoruv.scale(new Apfloat(1).divide(m4));
        pvector3 = pvector3.minus(pvectoruv2);
        pvector3 = pvector3.scale(m4);
        pvectoruv = pvector3copy.minus(pvector3);
        // end of additional UV photon
        VectorApfloat[] answer = new VectorApfloat[4];
        answer[0] = pvector3;
        answer[1] = pvectorx1;
        answer[2] = pvectore1ab;
        answer[3] = pvectoruv;
        // System.out.println(ApfloatMath.pow(pvectore1ab.magnitude(), 2).divide(new Apfloat(2))
        //         .divide(electron));
        return answer;
    }

    // generate random polar angle for x-ray detector
    public Apfloat randomThetaX() {
        double a = StdRandom.uniform();
        if (a > 0.5) return new Apfloat(StdRandom.uniform() * 0.293176);
        else
            return new Apfloat(StdRandom.uniform() * 0.293176 + Math.PI - 0.293176);
    }

    // reconstruction from momentum and energy
    public Apfloat reconstruct(VectorApfloat[] input) {
        VectorApfloat pvector3 = input[0];
        VectorApfloat pvectorx1 = input[1];
        VectorApfloat pvectore1ab = input[2];
        VectorApfloat puv = input[3];
        Apfloat a = pvectore1ab.magnitude();
        Apfloat b = pvector3.magnitude();
        Apfloat ione = ApfloatMath.sqrt(b.multiply(b).add(m4.multiply(m4)));
        Apfloat electrone = ApfloatMath.sqrt(a.multiply(a).add(electron.multiply(electron)));
        Apfloat m = ApfloatMath
                .pow((csGround.multiply(convertToeV)).subtract(electrone)
                                                     .subtract(pvectorx1.magnitude())
                                                     .subtract(ione), 2);
        m = m.subtract(
                ApfloatMath
                        .pow(pvector3.plus(pvectorx1).plus(pvectore1ab).magnitude(), 2));
        return m.divide(new Apfloat(1000000));
    }

    public Apfloat realmomentum(Apfloat nu) {
        return twoBody(csGround.multiply(convertToeV), m2.add(new Apfloat(34430)), nu);
    }

    // 2 xray M-fill N-fill
    public Apfloat[] kcapture2(Apfloat nu) {
        m11 = m2.add(new Apfloat(550));
        m1 = m11.add(new Apfloat(33880));
        nu = nu.multiply(new Apfloat(1000));
        Apfloat p = twoBody(csGround.multiply(convertToeV), m1, nu);
        Apfloat theta1 = randomTheta1();
        Apfloat theta2 = randomTheta2();
        VectorApfloat pvector1 = formVector(theta1, theta2, p);
        // end of first decay
        Apfloat theta3 = randomThetaX();
        Apfloat theta4 = randomTheta2();
        //   Apfloat theta41 = randomTheta2();
        Apfloat p1 = twoBody(m1, m11, new Apfloat(0, 500));
        VectorApfloat pvector2 = pvector1.scale(new Apfloat(1).divide(m1));
        VectorApfloat pvectorx1 = formVector(theta3, theta4, p1);
        // VectorApfloat pvectorfalse = formVector(theta31, theta41,
        //               p1); // false x-ray from another event
        VectorApfloat pvectorx = pvectorx1.scale(new Apfloat(1).divide(m11));
        pvector2 = pvector2.minus(pvectorx);
        pvector2 = pvector2.scale(m11);
        pvectorx1 = pvector1.minus(pvector2); // conserve momentum in lab frame
        // end of second decay (first x-ray)
        Apfloat theta5 = randomTheta1();
        Apfloat theta6 = randomTheta2();
        VectorApfloat pvector2copy = pvector2;
        pvector2 = pvector2.scale(new Apfloat(1).divide(m11));
        Apfloat p2 = twoBody(m11, m2, new Apfloat(0, 500));
        VectorApfloat pvectorx12 = formVector(theta5, theta6, p2);
        VectorApfloat pvectorx2 = pvectorx12.scale(new Apfloat(1).divide(m2));
        pvector2 = pvector2.minus(pvectorx2);
        pvector2 = pvector2.scale(m2);
        //  pvectorx12 = pvector2copy.minus(pvector2);
        // end of additional decay (second x-ray)
        Apfloat pboost = twoBody(m2, m3, electron);
        Apfloat theta7 = randomTheta1();
        Apfloat theta8 = randomTheta2();
        VectorApfloat pvector3 = pvector2.scale(new Apfloat(1).divide(m2));
        VectorApfloat intere = pvector3.scale(electron);
        // electron momentum in Xe restframe
        VectorApfloat pvectore1 = formVector(theta7, theta8, pboost);
        // electron momentum in lab frame
        VectorApfloat pvectore1ab = pvectore1.plus(intere);
        VectorApfloat pvectore = pvectore1.scale(new Apfloat(1).divide(m3));
        pvector3 = pvector3.minus(pvectore);
        pvector3 = pvector3.scale(m3);
        pvectore1ab = pvector2.minus(pvector3);
        // end of third decay
        Apfloat theta9 = randomTheta1();
        Apfloat theta10 = randomTheta2();
        VectorApfloat pvector3copy = pvector3;
        pvector3 = pvector3.scale(new Apfloat(1).divide(m3));
        Apfloat p3 = twoBody(m3, m4, new Apfloat(0, 500));
        VectorApfloat pvectoruv = formVector(theta9, theta10, p3);
        VectorApfloat pvectoruv2 = pvectoruv.scale(new Apfloat(1).divide(m4));
        pvector3 = pvector3.minus(pvectoruv2);
        pvector3 = pvector3.scale(m4);
        pvectoruv = pvector3copy.minus(pvector3);
        // end of additional UV photon
        Apfloat ione = ApfloatMath.pow(pvector3.magnitude(), 2).divide(new Apfloat(2)).divide(m4);
        VectorApfloat pvectorxs = pvectorx1.scale(new Apfloat(34430).divide(pvectorx1.magnitude()));
        Apfloat electronK = ApfloatMath.sqrt(electron.multiply(electron).add(pvectore1ab.magnitude()
                                                                                        .multiply(
                                                                                                pvectore1ab
                                                                                                        .magnitude())));
        Apfloat m = ApfloatMath
                .pow(csGround.multiply(convertToeV).subtract(m4).subtract(electronK)
                             .subtract(new Apfloat(34430))
                             .subtract(ione), 2);
        m = m.subtract(
                ApfloatMath.pow((pvector3
                        .plus(pvectorxs)
                        .plus(pvectore1ab)).magnitude(), 2));
        Apfloat[] answer = new Apfloat[3];
        // answer[0] store the reconstructed mass square
        answer[0] = m.divide(new Apfloat(1000000));
        // answer[1] stores the (reconstructed) magnitude of neutrino momentum
        answer[1] = pvector3.plus(pvectorx1).plus(pvectore1ab).magnitude();
        // answer[2] stores the difference btw the real and the reconstructed neutrino momentum
        if (answer[0].compareTo(new Apfloat(0)) > 0) {
            answer[2] = answer[1]
                    .subtract(realmomentum(
                            ApfloatMath.sqrt(answer[0].multiply(new Apfloat(1000000)))));
        }
        return answer;
    }

    public static void main(String[] args) {
        KcapTure test = new KcapTure();
        int n = 160000;
        double[] massStandard = new double[n];
        for (int i = 0; i < n; i++) {
            double random = StdRandom.uniform();
            if (i < 7) {
                Apfloat[] result = test.kcapture2(new Apfloat(0));
                massStandard[i] = Double.parseDouble(result[0].toString());
            }
            else
                massStandard[i] = Double
                        .parseDouble(test.reconstruct(test.kcapturep(new Apfloat(0))).toString());
        }
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.FREQUENCY);
        int number = 100;
        // dataset.addSeries("Neutrino mass square", mSquared, number);
        dataset.addSeries("mass-squared (keV\u00b2)", massStandard, number);
        // dataset.addSeries("33.6mistaken", pfake, number);
        String plotTitle
                = "";
        String xaxis = "mass-squared (keV\u00b2)"; // ²
        String yaxis = "occurence";
        PlotOrientation orientation = PlotOrientation.VERTICAL;
        JFreeChart chart = ChartFactory.createHistogram(plotTitle, xaxis, yaxis,
                                                        dataset, orientation, false, false, false);
        chart.getTitle().setFont(new Font("TimesRoman", Font.ITALIC, 20));
        // chart.getLegend().setItemFont(new Font("TimesRoman", Font.ITALIC, 10));
        String stdvstr = String.format("mean  : %.2f stdv : %.2f", StdStats.mean(massStandard),
                                       StdStats.stddev(massStandard));
        TextTitle a = new TextTitle(stdvstr);
        a.setFont(new Font("TimesRoman", Font.ITALIC, 15));
        chart.addSubtitle(a);
        XYPlot plot = chart.getXYPlot();
        Paint[] paintArray = {
                Color.BLUE,
                new Color(0x80ff0000, true),
                new Color(0x800000ff, true)
        };
        plot.setDrawingSupplier(new DefaultDrawingSupplier(
                paintArray,
                DefaultDrawingSupplier.DEFAULT_FILL_PAINT_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE));
        plot.getDomainAxis().setLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        plot.getDomainAxis().setTickLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        plot.getRangeAxis().setLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        plot.getRangeAxis().setTickLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        domain.setRange(0, 700);
        NumberAxis range = (NumberAxis) plot.getRangeAxis();
        range.setRange(0, 10);
        int scale = 10;
        try {
            ChartUtilities
                    .writeScaledChartAsPNG(
                            new FileOutputStream(new File("mi.png")),
                            chart,
                            500, 300,
                            scale, scale);
        }
        catch (IOException e) {
            System.out.println("error");
        }
        // XYSeriesCollection dataset1 = new XYSeriesCollection();
        // XYSeries fake = new XYSeries("momentum magnitude difference(mis-id)");
        // for (int i = 0; i < n; i++) {
        //     Apfloat[] input2 = test.kcapture2(new Apfloat(0));
        //     if (input2[2] != null) {
        //         fake.add(Double.parseDouble(input2[0].toString()),
        //                  Double.parseDouble(input2[2].toString()));
        //     }
        // }
        // dataset1.addSeries(fake);
        // JFreeChart chart2 = ChartFactory.createScatterPlot(
        //         "Momentum Deviation",
        //         "mass square(keV\u00b2)", "momentum difference(eV/c)", dataset1);
        // Shape shape = new Ellipse2D.Double(0, 0, 2, 2);
        // XYPlot xyPlot = (XYPlot) chart2.getPlot();
        // XYItemRenderer renderer = xyPlot.getRenderer();
        // renderer.setBaseShape(shape);
        // renderer.setBasePaint(Color.red);
        // renderer.setSeriesShape(0, shape);
        // try {
        //     ChartUtilities
        //             .writeScaledChartAsPNG(
        //                     new FileOutputStream(new File("testing.png")),
        //                     chart2,
        //                     500, 300,
        //                     scale, scale);
        // }
        // catch (IOException e) {
        //     System.out.println("error");
        // }
        // Apfloat[] mu = new Apfloat[n];
        // for (int i = 0; i < n; i++) {
        //     mu[i] = new Apfloat(i * 300.0 / n);
        // }
        // XYSeriesCollection dataset2 = new XYSeriesCollection();
        // XYSeries massmom = new XYSeries("mass-square-momentum");
        // for (int i = 0; i < n; i++) {
        //     Apfloat[] input2 = test.kcapture2(ApfloatMath.sqrt(mu[i]));
        //     Apfloat p = input2[1];
        //     p = p.divide(new Apfloat(1000));
        //     massmom.add(input2[0],
        //                 Double.parseDouble(p.multiply(p).toString()));
        // }
        // dataset2.addSeries(massmom);
        // JFreeChart chart3 = ChartFactory.createScatterPlot(
        //         "p\u00b2-m\u00b2",
        //         "mass square(keV\u00b2)", "momentum square(keV\u00b2)", dataset2);
        // XYPlot xyPlot3 = (XYPlot) chart3.getPlot();
        // XYItemRenderer renderer3 = xyPlot3.getRenderer();
        // renderer3.setBaseShape(shape);
        // renderer3.setBasePaint(Color.red);
        // renderer3.setSeriesShape(0, shape);
        // try {
        //     ChartUtilities
        //             .writeScaledChartAsPNG(
        //                     new FileOutputStream(new File("testingmass2.png")),
        //                     chart3,
        //                     500, 300,
        //                     scale, scale);
        // }
        // catch (IOException e) {
        //     System.out.println("error");
        // }


    }
}
