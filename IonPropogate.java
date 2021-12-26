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
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.XYSeries;
import org.jfree.ui.RectangleInsets;

import java.awt.Color;
import java.awt.Font;
import java.awt.Paint;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

public class IonPropogate {
    private static double s = 0.3; // acceleration length in m 0.1
    private static double d = 0.626; // drift length in m 0.826
    private static double v = 120; // extraction field in volts 20
    private static double q = 1.6 * Math.pow(10, -19); // charge of ion
    private static double m = 1.21936 * Math.pow(10, 11) * 1.79 * Math.pow(10, -36);
    // final stable xe ion mass in kg
    private static double convertTokg = 5.36 * Math.pow(10, -28); // convert ev/c to kg m/s
    private static double convertkgToeV = 5.6095883571872 * Math.pow(10, 35); // convert kg to eV
    private static final double electron = 9.10938 * Math.pow(10, -31); // electron mass in kg
    private static final double magF = 8 * Math.pow(10, -4); // magnetic field in Tesla
    private static final double TOF_ION_LOW = 182 * Math.pow(10, -6);
    // time window ion tof (low bound) in micro s
    private static final double TOF_ION_HI = 194 * Math.pow(10, -6);
    // time window ion tof( high bound) in mirco s
    private static final double TOF_E_LOW = 100 * Math.pow(10, -9);
    private static final double TOF_E_HI = 164 * Math.pow(10, -9);

    // get time of flight from longitudinal momentum
    // p in SI unit
    public static double toF(double p) {
        return s * (-p / (q * v) + Math.sqrt(p * p / (q * q * v * v) + 2 * m / (q * v)))
                + d * Math.pow(p * p / (m * m) + 2 * q * v / m, -0.5);

    }

    // get time of flight of the electron from longitudinal momentum
    // p in SI unit
    public static double toFe(double p) {
        return (Math.sqrt(p * p + 2 * 424.0 * electron * q) - p) / (q * 424.0);
    }

    // reconstruct electron longitudinal momentum
    public static double toFde(double t) {
        return (1 - q * 424.0 * t * t / (2 * electron)) * electron / t;
    }


    // derivative of toF equation
    public static double toFd(double p) {
        return -s / (q * v) + Math.pow(p * p / (q * q + v * v) + 2 * m / (q * v), -0.5) * p * s / (q
                * q * v * v) - p * d *
                Math.pow(p * p / (m * m) + 2 * q * v / m, -1.5) / (m * m);
    }

    // get the momentum from TOF
    public static double reverseToF(double t, double pguess) {
        for (int i = 0; i < 100; i++) {
            pguess = pguess - (toF(pguess) - t) / toFd(pguess);
        }
        return pguess;
    }

    // get the hitting position of the electron
    public static double[] tranE(double px, double py, double t) {
        double[] ans = new double[2];
        ans[0] =
                2 * Math.sqrt(px * px + py * py) * Math.abs(Math.sin(q * magF * t / (electron * 2)))
                        / (q * magF);
        ans[1] = Math.atan(py / px) + ((q * magF * t / electron) - 2 * Math.PI * Math
                .floor(t / (2 * Math.PI * electron / (q * magF)))) / 2;
        if (px < 0) {
            ans[1] = ans[1] + Math.PI;
        }

        return ans;

    }

    // reconstruct the transverse momentum of the electron
    public static double[] reverseT(double r, double theta, double t) {
        double pnorm = r * q * magF / (2 * Math.abs(Math.sin(q * magF * t / (electron * 2))));
        theta = theta - ((q * magF * t / electron) - 2 * Math.PI * Math
                .floor(t / (2 * Math.PI * electron / (q * magF)))) / 2;
        double p1 = pnorm * Math.cos(theta);
        double p2 = pnorm * Math.sin(theta);
        double[] ptr = new double[2];
        ptr[0] = p1;
        ptr[1] = p2;
        return ptr;
    }

    public static void main(String[] args) {
        // try {
        //     System.setOut(
        //             new PrintStream(
        //                     new BufferedOutputStream(new FileOutputStream("m2.txt"))));
        // }
        // catch (FileNotFoundException ex) {
        //     System.out.println("error");
        // }
        KcapTure test = new KcapTure();
        int n = 10000;
        Queue<Double> mSquared = new Queue<Double>();
        Queue<Double> mSquared2 = new Queue<Double>();
        XYSeries tofp = new XYSeries('k');
        double[] tofs = new double[n];
        int size = 0;
        int size2 = 0;
        //  double[] time2 = new double[n];
        for (int i = 0; i < n; i++) {
            VectorApfloat[] input = test.kcapturep(new Apfloat(0));
            double p1 = Double
                    .parseDouble(input[0].cartesian(2).toString()); // longitudinal momentum
            double time1 = toF(p1 * convertTokg);
            VectorApfloat pElectron = input[2];
            double toFe = toFe(
                    Double.parseDouble(
                            pElectron.cartesian(2)
                                     .toString())
                            * convertTokg);
            tofs[i] = toFe * 1000000000;
            double randomLow = TOF_E_LOW - toFe;
            double randomHi = TOF_E_HI - toFe;
            // double tau = (randomHi - randomLow) * StdRandom.uniform() + randomLow;
            double tau = 0;
            if (toFe + tau > 131 * Math.pow(10, -9) && toFe + tau < 137 * Math.pow(10, -9)) {
                continue;
            }
            double pfake = reverseToF(time1 + tau,
                                      (p1 + 10000) * convertTokg); // in kg
            pfake = pfake / convertTokg; // back to eV/c
            Apfloat[] form = new Apfloat[3];
            form[0] = input[0].cartesian(0).multiply(new Apfloat(time1))
                              .divide(new Apfloat(time1 + tau));
            form[1] = input[0].cartesian(1).multiply(new Apfloat(time1))
                              .divide(new Apfloat(time1 + tau));
            form[2] = new Apfloat(pfake);
            input[0] = new VectorApfloat(form);
            Apfloat theta31 = test.randomThetaX();
            Apfloat theta4 = test.randomTheta2();
            VectorApfloat pvectorx1 = test.formVector(theta31, theta4, input[1].magnitude());
            input[1] = pvectorx1;
            double peL = toFde(toFe + tau) / convertTokg;
            double px = Double.parseDouble(pElectron.cartesian(0).toString()) * convertTokg;
            double py = Double.parseDouble(pElectron.cartesian(1).toString()) * convertTokg;
            double[] ptrans = tranE(px, py, toFe);
            double r = ptrans[0];
            if (r < 0.003 || r > 0.06) {
                continue;
            }
            double theta = ptrans[1];
            double[] pxyreconstruct = reverseT(r, theta, toFe + tau);
            double pxr = pxyreconstruct[0] / convertTokg;
            double pyr = pxyreconstruct[1] / convertTokg;
            Apfloat[] pdata = new Apfloat[3];
            pdata[0] = new Apfloat(pxr);
            pdata[1] = new Apfloat(pyr);
            pdata[2] = new Apfloat(peL);
            input[2] = new VectorApfloat(pdata);
            double mass = Double.parseDouble(test.reconstruct(input).toString());
            mSquared2.enqueue(mass);
            size2++;
            double Ps = Double.parseDouble(
                    input[0].plus(input[1]).plus(input[2]).magnitude().toString());
            if (Ps * Ps / 1000000 + mass < 320.2 * 320.2 - 9.6
                    || Ps * Ps / 1000000 + mass > 320.2 * 320.2 + 9.6) {
                continue;
            }
            // StdOut.println(mass);
            mSquared.enqueue(mass);
            size++;
        }
        StdOut.println(size);
        double[] massq = new double[size];
        for (int i = 0; i < size; i++) {
            massq[i] = mSquared.dequeue();
        }
        double[] massq2 = new double[size2];
        for (int i = 0; i < size2; i++) {
            massq2[i] = mSquared2.dequeue();
        }
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.FREQUENCY);
        int number = 100;
        // dataset.addSeries("Neutrino mass square", mSquared, number);
        dataset.addSeries("time", massq, number);
        // dataset.addSeries("kinetic energy of electron + kinetic energy of ion", kSquared, number
        //);
        // dataset.addSeries("33.6keV x-ray", time2, number);
        String plotTitle
                = "";
        String xaxis = "mass-squared (keV\u00b2)";
        String yaxis = "occurence";
        PlotOrientation orientation = PlotOrientation.VERTICAL;
        JFreeChart chart = ChartFactory.createHistogram(plotTitle, xaxis, yaxis,
                                                        dataset, orientation, false, false, false);
        chart.setPadding(new RectangleInsets(10, 10, 10, 10));
        chart.getTitle().setFont(new Font("TimesRoman", Font.ITALIC, 20));
        String stdvstr = String
                .format("mean : %.2f, stdv : %.2f", StdStats.mean(massq), StdStats.stddev(massq));
        TextTitle a = new TextTitle(stdvstr);
        a.setFont(new Font("TimesRoman", Font.ITALIC, 15));
        chart.addSubtitle(a);
        XYPlot plot = chart.getXYPlot();
        // NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        // NumberAxis range = (NumberAxis) plot.getRangeAxis();
        // domain.setRange(0, 10000);
        // range.setRange(0, 30);
        Paint[] paintArray = {
                Color.BLUE,
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
        int scale = 3;
        try {
            ChartUtilities
                    .writeScaledChartAsPNG(
                            new FileOutputStream(new File("accid1.png")),
                            chart,
                            500, 300,
                            scale, scale);
        }
        catch (IOException e) {
            System.out.println("error");
        }
        HistogramDataset dataset2 = new HistogramDataset();
        dataset2.setType(HistogramType.FREQUENCY);
        // dataset.addSeries("Neutrino mass square", mSquared, number);
        dataset2.addSeries("time", massq2, number);
        // dataset.addSeries("kinetic energy of electron + kinetic energy of ion", kSquared, number
        //);
        // dataset.addSeries("33.6keV x-ray", time2, number);
        JFreeChart chart2 = ChartFactory.createHistogram(plotTitle, xaxis, yaxis,
                                                         dataset2, orientation, false, false,
                                                         false);
        chart2.setPadding(new RectangleInsets(10, 10, 10, 10));
        chart2.getTitle().setFont(new Font("TimesRoman", Font.ITALIC, 20));
        String stdvstr2 = String
                .format("mean : %.2f, stdv : %.2f", StdStats.mean(massq2), StdStats.stddev(massq2));
        TextTitle a2 = new TextTitle(stdvstr2);
        a2.setFont(new Font("TimesRoman", Font.ITALIC, 15));
        chart2.addSubtitle(a2);
        XYPlot plot2 = chart2.getXYPlot();
        // NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        // NumberAxis range = (NumberAxis) plot.getRangeAxis();
        // domain.setRange(0, 10000);
        // range.setRange(0, 30);
        plot2.setDrawingSupplier(new DefaultDrawingSupplier(
                paintArray,
                DefaultDrawingSupplier.DEFAULT_FILL_PAINT_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
                DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE));
        plot2.getDomainAxis().setLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        plot2.getDomainAxis().setTickLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        plot2.getRangeAxis().setLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        plot2.getRangeAxis().setTickLabelFont(new Font("TimesRoman", Font.ITALIC, 10));
        try {
            ChartUtilities
                    .writeScaledChartAsPNG(
                            new FileOutputStream(new File("accid2.png")),
                            chart2,
                            500, 300,
                            scale, scale);
        }
        catch (IOException e) {
            System.out.println("error");
        }
        // double[] testmass = new double[n];
        // for (int i = 0; i < n; i++) {
        //     testmass[i] = i * 20.0 / n; // neutrino mass in keV
        // }
        // XYSeriesCollection testset = new XYSeriesCollection();
        // XYSeries pm = new XYSeries("Momentum square and Mass square");
        // XYSeries f1 = new XYSeries("-9.6");
        // XYSeries f2 = new XYSeries("+9.6");
        // for (int i = 0; i < n; i++) {
        //     VectorApfloat[] momenta = test.kcapturep(new Apfloat(testmass[i]));
        //     double plong = Double.parseDouble(momenta[0].cartesian(2).toString());
        //     double TOF = toF(plong * convertTokg);
        //     double preconstruct = reverseToF(TOF, (plong + 10000) * convertTokg);
        //     preconstruct = preconstruct / convertTokg; // back to eV/c
        //     Apfloat[] p = new Apfloat[3];
        //     p[0] = momenta[0].cartesian(0);
        //     p[1] = momenta[0].cartesian(1);
        //     p[2] = new Apfloat(preconstruct);
        //     momenta[0] = new VectorApfloat(p);
        //     double pnu = Double.parseDouble(
        //             momenta[0].plus(momenta[1]).plus(momenta[2]).magnitude().toString());
        //     pnu = pnu * pnu / 1000000;
        //     // System.out.print(pnu);
        //     // System.out.print("   ");
        //     double msquare = testmass[i] * testmass[i];
        //     double msquare1 = Double.parseDouble(test.reconstruct(momenta).toString());
        //     // System.out.println(msquare1);
        //     // pm.add(msquare1, pnu);
        //     //   if (pnu > 320.178 * 320.178 - 9.6 - msquare && 320.178 * 320.178 + 9.6 - msquare > pnu)
        //     //       size++;
        //     // f1.add(msquare1, 320.1718 * 320.1718 - 9.6 - msquare1);
        //     // f2.add(msquare1, 320.1718 * 320.1718 + 9.6 - msquare1);
        // }
        // // testset.addSeries(pm);
        // // testset.addSeries(f1);
        // // testset.addSeries(f2);
        // testset.addSeries(tofp);
        // JFreeChart chart2 = ChartFactory.createXYLineChart(
        //         "Q value",
        //         "mass square(keV\u00b2)", "p square(keV\u00b2)", testset);
        // // chart2.addSubtitle(new TextTitle(String.valueOf(size)));
        // Shape shape = new Ellipse2D.Double(0.2, 0.2, 0.5, 0.5);
        // XYPlot xyPlot = (XYPlot) chart2.getPlot();
        // NumberAxis domain2 = (NumberAxis) xyPlot.getDomainAxis();
        // //  domain2.setRange(-1000, 1000);
        // // NumberAxis range = (NumberAxis) xyPlot.getRangeAxis();
        // // range.setRange(102100, 102800);
        // XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
        // renderer2.setSeriesLinesVisible(0, false);
        // renderer2.setSeriesShapesVisible(0, true);
        // renderer2.setSeriesShapesVisible(1, false);
        // renderer2.setSeriesShapesVisible(2, false);
        // renderer2.setSeriesShape(0, shape);
        // xyPlot.setRenderer(renderer2);
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


    }
}
