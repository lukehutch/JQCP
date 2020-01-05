import org.joml.Matrix3d;
import org.joml.Vector3d;

/*******************************************************************************
 * Code translated from C to Java by Luke Hutchison,
 * and ported to JOML linear algebra library ( https://github.com/JOML-CI/JOML )
 * 
 * Original code: https://theobald.brandeis.edu/qcp/
 * 
 * -----------------------------------------------------------------------------
 * 
 * <pre>
 * 
 *  -/_|:|_|_\-
 *
 *  File:           qcprot.c
 *  Version:        1.5
 *
 *  Function:       Rapid calculation of the least-squares rotation using a
 *                  quaternion-based characteristic polynomial and
 *                  a cofactor matrix
 *
 *  Author(s):      Douglas L. Theobald
 *                  Department of Biochemistry
 *                  MS 009
 *                  Brandeis University
 *                  415 South St
 *                  Waltham, MA  02453
 *                  USA
 *
 *                  dtheobald@brandeis.edu
 *
 *                  Pu Liu
 *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                  665 Stockton Drive
 *                  Exton, PA  19341
 *                  USA
 *
 *                  pliu24@its.jnj.com
 *
 *
 *    If you use this QCP rotation calculation method in a publication, please
 *    reference:
 *
 *      Douglas L. Theobald (2005)
 *      "Rapid calculation of RMSD using a quaternion-based characteristic
 *      polynomial."
 *      Acta Crystallographica A 61(4):478-480.
 *
 *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *      "Fast determination of the optimal rotational matrix for macromolecular
 *      superpositions."
 *      Journal of Computational Chemistry 31(7):1561-1563.
 *
 *  Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *    of conditions and the following disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *    endorse or promote products derived from this software without specific prior written
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    2009/04/13      Started source
 *    2010/03/28      Modified FastCalcRMSDAndRotation() to handle tiny qsqr
 *                    If trying all rows of the adjoint still gives too small
 *                    qsqr, then just return identity matrix. (DLT)
 *    2010/06/30      Fixed prob in assigning A[9] = 0 in InnerProduct()
 *                    invalid mem access
 *    2011/02/21      Made CenterCoords use weights
 *    2011/05/02      Finally changed CenterCoords declaration in qcprot.h
 *                    Also changed some functions to static
 *    2011/07/08      Put in fabs() to fix taking sqrt of small neg numbers, fp error
 *    2012/07/26      Minor changes to comments and main.c, more info (v.1.4)
 *    2016/07/13      Fixed normalization of RMSD in FastCalcRMSDAndRotation(), should divide by
 *                    sum of weights (thanks to Geoff Skillman)
 *    2020/01/04      Ported to Java by Luke Hutchison
 * </pre>
 ******************************************************************************/

public class JQCP_JOML {
    /**
     * Calculate the inner product of two structures. If weight array is not null, calculate the weighted inner
     * product.
     * 
     * Warning:
     * 
     * 1. You MUST center the structures, coords1 and coords2, before calling this function.
     * 
     * 2. Please note how the structure coordinates are stored in the double **coords arrays. They are 3xN arrays,
     * not Nx3 arrays as is also commonly used (where the x, y, z axes are interleaved). The difference is something
     * like this for storage of a structure with 8 atoms:
     *
     * <pre>
     *   Nx3: xyzxyzxyzxyzxyzxyzxyzxyz
     *   3xN: xxxxxxxxyyyyyyyyzzzzzzzz
     * </pre>
     * 
     * The functions can be easily modified, however, to accomodate any data format preference. I chose this format
     * because it is readily used in vectorized functions (SIMD, Altivec, MMX, SSE2, etc.).
     * 
     * @param coords1 reference structure
     * @param coords2 candidate structure
     * @param len     the size of the system
     * @param weight  the weight array of size len, or null to use unweighted coordinates.
     * 
     * @param AOut    out parameter (double[9]): the inner product matrix
     * 
     * @return (G1 + G2) * 0.5; used as E0 in function
     *         {@link #FastCalcRMSDAndRotation(double[], double[], double[], double, double, double)}.
     */
    private static double innerProduct(Vector3d[] coords1, Vector3d[] coords2, int len, double[] weight,
            double[] AOut) {
        double G1 = 0.0, G2 = 0.0;

        AOut[0] = AOut[1] = AOut[2] = AOut[3] = AOut[4] = AOut[5] = AOut[6] = AOut[7] = AOut[8] = 0.0;

        if (weight != null) {
            for (int i = 0; i < len; ++i) {
                double x1 = weight[i] * coords1[i].x;
                double y1 = weight[i] * coords1[i].y;
                double z1 = weight[i] * coords1[i].z;

                G1 += x1 * coords1[i].x + y1 * coords1[i].y + z1 * coords1[i].z;

                double x2 = coords2[i].x;
                double y2 = coords2[i].y;
                double z2 = coords2[i].z;

                G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

                AOut[0] += (x1 * x2);
                AOut[1] += (x1 * y2);
                AOut[2] += (x1 * z2);

                AOut[3] += (y1 * x2);
                AOut[4] += (y1 * y2);
                AOut[5] += (y1 * z2);

                AOut[6] += (z1 * x2);
                AOut[7] += (z1 * y2);
                AOut[8] += (z1 * z2);
            }
        } else {
            for (int i = 0; i < len; ++i) {
                double x1 = coords1[i].x;
                double y1 = coords1[i].y;
                double z1 = coords1[i].z;

                G1 += x1 * x1 + y1 * y1 + z1 * z1;

                double x2 = coords2[i].x;
                double y2 = coords2[i].y;
                double z2 = coords2[i].z;

                G2 += (x2 * x2 + y2 * y2 + z2 * z2);

                AOut[0] += (x1 * x2);
                AOut[1] += (x1 * y2);
                AOut[2] += (x1 * z2);

                AOut[3] += (y1 * x2);
                AOut[4] += (y1 * y2);
                AOut[5] += (y1 * z2);

                AOut[6] += (z1 * x2);
                AOut[7] += (z1 * y2);
                AOut[8] += (z1 * z2);
            }
        }

        return (G1 + G2) * 0.5;
    }

    /**
     * Calculate the RMSD, and/or the optimal rotation matrix.
     * 
     * @param A        the inner product of two structures (double[9])
     * @param E0       (G1 + G2) * 0.5
     * @param len      the size of the system
     * @param minScore if( minScore > 0 && rmsd < minScore) then calculate only the rmsd; otherwise, calculate both
     *                 the RMSD & the rotation matrix
     * 
     * @param rotOut   out parameter: the rotation matrix
     * @param rmsdOut  out parameter (double[1]): the RMSD value
     * 
     * @return a value less than or equal to 0 if only the rmsd was calculated (including if the rotation was
     *         smaller than the threshold of precision); or a value greater than 0 if both the RMSD and rotational
     *         matrix were calculated.
     */
    public static int fastCalcRMSDAndRotation(double[] A, double E0, double len, double minScore, Matrix3d rotOut,
            double[] rmsdOut) {

        double evecprec = 1e-6;
        double evalprec = 1e-11;
        int maxIter = 50;

        double Sxx = A[0];
        double Sxy = A[1];
        double Sxz = A[2];
        double Syx = A[3];
        double Syy = A[4];
        double Syz = A[5];
        double Szx = A[6];
        double Szy = A[7];
        double Szz = A[8];

        double Sxx2 = Sxx * Sxx;
        double Syy2 = Syy * Syy;
        double Szz2 = Szz * Szz;

        double Sxy2 = Sxy * Sxy;
        double Syz2 = Syz * Syz;
        double Sxz2 = Sxz * Sxz;

        double Syx2 = Syx * Syx;
        double Szy2 = Szy * Szy;
        double Szx2 = Szx * Szx;

        double SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
        double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

        double C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
        double C1 = 8.0 * (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx - Sxx * Syy * Szz - Syz * Szx * Sxy
                - Szy * Syx * Sxz);

        double SxzpSzx = Sxz + Szx;
        double SyzpSzy = Syz + Szy;
        double SxypSyx = Sxy + Syx;
        double SyzmSzy = Syz - Szy;
        double SxzmSzx = Sxz - Szx;
        double SxymSyx = Sxy - Syx;
        double SxxpSyy = Sxx + Syy;
        double SxxmSyy = Sxx - Syy;
        double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

        double C0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
                + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
                + (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz))
                        * (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz))
                + (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz))
                        * (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz))
                + (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz))
                        * (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz))
                + (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz))
                        * (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz));

        /* Newton-Raphson */
        double mxEigenV = E0;
        double oldg = 0.0;
        for (int i = 0; i < maxIter; ++i) {
            oldg = mxEigenV;
            double x2 = mxEigenV * mxEigenV;
            double b = (x2 + C2) * mxEigenV;
            double a = b + C1;
            double delta = ((a * mxEigenV + C0) / (2.0 * x2 * mxEigenV + b + a));
            mxEigenV -= delta;
            /*
             * System.out.printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV);
             */
            if (Math.abs(mxEigenV - oldg) < Math.abs(evalprec * mxEigenV)) {
                break;
            }
            if (i == maxIter - 1) {
                System.err.printf("\nMore than %d iterations needed!\n", maxIter);
            }
        }

        /* the abs() is to guard against extremely small, but *negative* numbers due to floating point error */
        double rms = Math.sqrt(Math.abs(2.0 * (E0 - mxEigenV) / len));
        rmsdOut[0] = rms;
        /* System.out.printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len); */

        if (minScore > 0) {
            if (rms < minScore) {
                return -1; // Don't bother with rotation.
            }
        }

        double a11 = SxxpSyy + Szz - mxEigenV;
        double a12 = SyzmSzy;
        double a13 = -SxzmSzx;
        double a14 = SxymSyx;
        double a21 = SyzmSzy;
        double a22 = SxxmSyy - Szz - mxEigenV;
        double a23 = SxypSyx;
        double a24 = SxzpSzx;
        double a31 = a13;
        double a32 = a23;
        double a33 = Syy - Sxx - Szz - mxEigenV;
        double a34 = SyzpSzy;
        double a41 = a14;
        double a42 = a24;
        double a43 = a34;
        double a44 = Szz - SxxpSyy - mxEigenV;
        double a3344_4334 = a33 * a44 - a43 * a34;
        double a3244_4234 = a32 * a44 - a42 * a34;
        double a3243_4233 = a32 * a43 - a42 * a33;
        double a3143_4133 = a31 * a43 - a41 * a33;
        double a3144_4134 = a31 * a44 - a41 * a34;
        double a3142_4132 = a31 * a42 - a41 * a32;
        double q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
        double q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
        double q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
        double q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;

        double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

        /*
         * The following code tries to calculate another column in the adjoint matrix when the norm of the current
         * column is too small. Usually this block will never be activated. To be absolutely safe this should be
         * uncommented, but it is most likely unnecessary.
         */
        if (qsqr < evecprec) {
            q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
            q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
            q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
            q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;
            qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

            if (qsqr < evecprec) {
                double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
                double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
                double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

                q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
                q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
                q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
                q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
                qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

                if (qsqr < evecprec) {
                    q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                    q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                    q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                    q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

                    if (qsqr < evecprec) {
                        /* if qsqr is still too small, return the identity matrix. */
                        rotOut.identity();
                        return 0;
                    }
                }
            }
        }

        double normq = Math.sqrt(qsqr);
        q1 /= normq;
        q2 /= normq;
        q3 /= normq;
        q4 /= normq;

        double a2 = q1 * q1;
        double x2 = q2 * q2;
        double y2 = q3 * q3;
        double z2 = q4 * q4;

        double xy = q2 * q3;
        double az = q1 * q4;
        double zx = q4 * q2;
        double ay = q1 * q3;
        double yz = q3 * q4;
        double ax = q1 * q2;

        // Rotation matrix elements in row-major order
        double rotOut0 = a2 + x2 - y2 - z2; // xx
        double rotOut1 = 2 * (xy + az); // xy
        double rotOut2 = 2 * (zx - ay); // xz
        double rotOut3 = 2 * (xy - az); // yx
        double rotOut4 = a2 - x2 + y2 - z2; // yy
        double rotOut5 = 2 * (yz + ax); // yz
        double rotOut6 = 2 * (zx + ay); // zx
        double rotOut7 = 2 * (yz - ax); // zy
        double rotOut8 = a2 - x2 - y2 + z2; // zz

        // JOML elements are in column-major order
        rotOut.set(rotOut0, rotOut3, rotOut6, rotOut1, rotOut4, rotOut7, rotOut2, rotOut5, rotOut8);

        return 1;
    }

    /**
     * Center the coordinates.
     * 
     * Warning: If you are doing a full superposition (the usual least squares way), you MUST center each structure
     * first. That is, you must translate each structure so that its centroid is at the origin. You can use this
     * method for this.
     * 
     * @param coords the coordinates to center.
     * @param len    the number of coordinates to center.
     * @param weight the weight array of size len, or null to use unweighted coordinates.
     */
    public static void centerCoords(Vector3d[] coords, int len, double[] weight) {
        double xsum = 0.0, ysum = 0.0, zsum = 0.0;
        if (weight != null) {
            double wsum = 0.0;
            for (int i = 0; i < len; ++i) {
                xsum += weight[i] * coords[i].x;
                ysum += weight[i] * coords[i].y;
                zsum += weight[i] * coords[i].z;

                wsum += weight[i];
            }

            if (Math.abs(wsum) < 1.0e-30) {
                throw new IllegalArgumentException("Weights sum to zero");
            }

            xsum /= wsum;
            ysum /= wsum;
            zsum /= wsum;
        } else {
            for (int i = 0; i < len; ++i) {
                xsum += coords[i].x;
                ysum += coords[i].y;
                zsum += coords[i].z;
            }

            if (len == 0) {
                throw new IllegalArgumentException("len == 0");
            }

            xsum /= len;
            ysum /= len;
            zsum /= len;
        }

        for (int i = 0; i < len; ++i) {
            coords[i].x -= xsum;
            coords[i].y -= ysum;
            coords[i].z -= zsum;
        }
    }

    /**
     * Calculate the RMSD & rotational matrix. Superposition coords2 onto coords1 -- in other words, coords2 is
     * rotated, coords1 is held fixed.
     * 
     * @param coords1 reference structure -- will be centered in-place
     * @param coords2 candidate structure -- will be centered in-place
     * @param len     the size of the system
     * @param weight  the weight array of size len, or null to use unweighted coordinates.
     * 
     * @param rotOut  out parameter: the rotation matrix
     * 
     * @return RMSD value
     */
    public static double calcRMSDRotationalMatrix(Vector3d[] coords1, Vector3d[] coords2, int len, double[] weight,
            Matrix3d rotOut) {
        /* center the structures -- if precentered you can omit this step */
        centerCoords(coords1, len, weight);
        centerCoords(coords2, len, weight);

        double wsum;
        if (weight == null) {
            wsum = len;
        } else {
            wsum = 0.0;
            for (int i = 0; i < len; ++i) {
                wsum += weight[i];
            }
        }

        /* calculate the (weighted) inner product of two structures */
        double[] A = new double[9];
        double E0 = innerProduct(coords1, coords2, len, weight, A);

        /* calculate the RMSD & rotational matrix */
        double[] rmsd = new double[1];
        fastCalcRMSDAndRotation(A, E0, wsum, -1, rotOut, rmsd);

        return rmsd[0];
    }

    private static void printCoords(Vector3d[] coords, int len) {
        for (int i = 0; i < len; ++i) {
            System.out.printf("\n % 8.3f % 8.3f % 8.3f", coords[i].x, coords[i].y, coords[i].z);
        }
        System.out.println();
    }

    private static void Mat3Print(Matrix3d matrix) {
        System.out.printf("\n");
        System.out.printf(" [ % 14.8f % 14.8f % 14.8f ]\n", matrix.m00, matrix.m10, matrix.m20);
        System.out.printf(" [ % 14.8f % 14.8f % 14.8f ]\n", matrix.m01, matrix.m11, matrix.m21);
        System.out.printf(" [ % 14.8f % 14.8f % 14.8f ]\n", matrix.m02, matrix.m12, matrix.m22);
    }

    /**
     * Sample code to use the routine for fast RMSD & rotational matrix calculation. Note that we superposition
     * frag_b onto frag_a. For the example provided below, the minimum least-squares RMSD for the two 7-atom
     * fragments should be [according to the documentation] 0.719106 A (although both the C code and this java port
     * give 0.745016).
     * 
     * [According to the documentation] the rotation quaternion should be:
     * 
     * <pre>
        -8.620063e-01   3.435505e-01   1.242953e-01  -3.513814e-01
     * </pre>
     * 
     * And [according to the documentation] the corresponding 3x3 rotation matrix should be:
     * 
     * <pre>
        [     0.72216358     0.69118937    -0.02714790 ]
        [    -0.52038257     0.51700833    -0.67963547 ]
        [    -0.45572112     0.50493528     0.73304748 ]
     * </pre>
     * 
     * but actually both the C code and this java port give:
     * 
     * <pre>
       [     0.77227551     0.63510272    -0.01533190 ]
       [    -0.44544846     0.52413614    -0.72584914 ]
       [    -0.45295276     0.56738509     0.68768304 ]
     * </pre>
     */
    public static void main(String[] args) {
        int len = 7;
        Vector3d[] frag_a = new Vector3d[len];
        Vector3d[] frag_b = new Vector3d[len];

        frag_a[0] = new Vector3d(-2.803, -15.373, 24.556);
        frag_a[1] = new Vector3d(0.893, -16.062, 25.147);
        frag_a[2] = new Vector3d(1.368, -12.371, 25.885);
        frag_a[3] = new Vector3d(-1.651, -12.153, 28.177);
        frag_a[4] = new Vector3d(-0.440, -15.218, 30.068);
        frag_a[5] = new Vector3d(2.551, -13.273, 31.372);
        frag_a[6] = new Vector3d(0.105, -11.330, 33.567);

        frag_b[0] = new Vector3d(-14.739, -18.673, 15.040);
        frag_b[1] = new Vector3d(-12.473, -15.810, 16.074);
        frag_b[2] = new Vector3d(-14.802, -13.307, 14.408);
        frag_b[3] = new Vector3d(-17.782, -14.852, 16.171);
        frag_b[4] = new Vector3d(-16.124, -14.617, 19.584);
        frag_b[5] = new Vector3d(-15.029, -11.037, 18.902);
        frag_b[6] = new Vector3d(-18.577, -10.001, 17.996);

        double[] weight = new double[len];
        for (int i = 0; i < len; ++i) {
            weight[i] = i + 1.0;
        }

        System.out.printf("\nCoords before centering:\n");

        printCoords(frag_a, len);
        printCoords(frag_b, len);

        Matrix3d rotmat = new Matrix3d();
        double rmsd = calcRMSDRotationalMatrix(frag_a, frag_b, len, weight, rotmat);

        System.out.printf("\nCoords after centering:\n");

        printCoords(frag_a, len);
        printCoords(frag_b, len);

        System.out.printf("\nQCP rmsd: %f\n", rmsd);

        System.out.printf("\nQCP Rotation matrix:\n");
        Mat3Print(rotmat);

        /* apply rotation matrix */
        for (int i = 0; i < len; ++i) {
            rotmat.transform(frag_b[i]);
        }

        /* calculate euclidean distance */
        double euc_dist = 0.0;
        double wtsum = 0.0;
        for (int i = 0; i < len; ++i) {
            wtsum += weight[i];
            double tmp = frag_a[i].distanceSquared(frag_b[i]);
            euc_dist += weight[i] * tmp;
        }

        System.out.printf("\nCoords 2 after rotation:\n");
        printCoords(frag_b, len);

        System.out.printf("\nExplicit RMSD calculated from transformed coords: %f\n\n",
                Math.sqrt(euc_dist / wtsum));
    }
}
