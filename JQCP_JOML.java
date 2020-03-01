package slam.reg.rotation;

import java.util.Arrays;

import org.joml.AxisAngle4f;
import org.joml.Matrix3d;
import org.joml.Quaterniond;
import org.joml.Quaternionf;
import org.joml.Vector3d;
import org.joml.Vector3f;

/*******************************************************************************
 * Code translated from C to Java by Luke Hutchison, and ported to JOML linear algebra library.
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
     * product. You must calculate the centroid of the structures coords1 and coords2, before calling this function.
     * 
     * @param coords1      reference structure
     * @param centroid1    the centroid of coords1
     * @param coords2      candidate structure
     * @param centroid2    the centroid of coords2
     * @param idxsToUse    the indices in the coords1, coords2 and weight arrays to use.
     * @param numIdxsToUse the number of entries in idxsToUse.
     * @param weight       the weight array of size len, or null to use unweighted coordinates.
     * 
     * @param AOut         out parameter: the inner product matrix
     * 
     * @return (G1 + G2) * 0.5; used as E0 in function
     *         {@link #FastCalcRMSDAndRotation(double[], double[], double[], double, double, double)}.
     */
    private static double innerProduct(Vector3d[] coords1, Vector3d centroid1, Vector3d[] coords2,
            Vector3d centroid2, int[] idxsToUse, int numIdxsToUse, double[] weight, Matrix3d AOut) {
        double G1 = 0.0, G2 = 0.0;

        AOut.zero();

        if (weight != null) {
            for (int i = 0; i < numIdxsToUse; ++i) {
                int idx = idxsToUse[i];
                double x1 = coords1[idx].x - centroid1.x;
                double y1 = coords1[idx].y - centroid1.y;
                double z1 = coords1[idx].z - centroid1.z;

                double wx1 = weight[idx] * x1;
                double wy1 = weight[idx] * y1;
                double wz1 = weight[idx] * z1;

                G1 += wx1 * x1 + wy1 * y1 + wz1 * z1;

                double x2 = coords2[idx].x - centroid2.x;
                double y2 = coords2[idx].y - centroid2.y;
                double z2 = coords2[idx].z - centroid2.z;

                G2 += weight[idx] * (x2 * x2 + y2 * y2 + z2 * z2);

                AOut.m00 += (wx1 * x2);
                AOut.m10 += (wx1 * y2);
                AOut.m20 += (wx1 * z2);

                AOut.m01 += (wy1 * x2);
                AOut.m11 += (wy1 * y2);
                AOut.m21 += (wy1 * z2);

                AOut.m02 += (wz1 * x2);
                AOut.m12 += (wz1 * y2);
                AOut.m22 += (wz1 * z2);
            }
        } else {
            for (int i = 0; i < numIdxsToUse; ++i) {
                int idx = idxsToUse[i];

                double x1 = coords1[idx].x - centroid1.x;
                double y1 = coords1[idx].y - centroid1.y;
                double z1 = coords1[idx].z - centroid1.z;

                G1 += x1 * x1 + y1 * y1 + z1 * z1;

                double x2 = coords2[idx].x - centroid2.x;
                double y2 = coords2[idx].y - centroid2.y;
                double z2 = coords2[idx].z - centroid2.z;

                G2 += (x2 * x2 + y2 * y2 + z2 * z2);

                AOut.m00 += (x1 * x2);
                AOut.m10 += (x1 * y2);
                AOut.m20 += (x1 * z2);

                AOut.m01 += (y1 * x2);
                AOut.m11 += (y1 * y2);
                AOut.m21 += (y1 * z2);

                AOut.m02 += (z1 * x2);
                AOut.m12 += (z1 * y2);
                AOut.m22 += (z1 * z2);
            }
        }

        return (G1 + G2) * 0.5;
    }

    /**
     * Calculate the RMSD, and/or the optimal rotation matrix.
     * 
     * @param A        the inner product of two structures
     * @param E0       (G1 + G2) * 0.5
     * @param len      the size of the system
     * @param minScore if( minScore > 0 && rmsd < minScore) then calculate only the rmsd; otherwise, calculate both
     *                 the RMSD & the rotation matrix
     * 
     * @param rotOut   out parameter: the recovered rotation
     * @param rmsdOut  out parameter (double[1]): the RMSD value
     * 
     * @return a value less than or equal to 0 if only the rmsd was calculated (including if the rotation was
     *         smaller than the threshold of precision); or a value greater than 0 if both the RMSD and rotational
     *         matrix were calculated.
     */
    public static int fastCalcRMSDAndRotation(Matrix3d A, double E0, double len, double minScore,
            Quaterniond rotOut, double[] rmsdOut) {

        double evecprec = 1e-6;
        double evalprec = 1e-11;
        int maxIter = 50;

        double Sxx = A.m00;
        double Sxy = A.m10;
        double Sxz = A.m20;
        double Syx = A.m01;
        double Syy = A.m11;
        double Syz = A.m21;
        double Szx = A.m02;
        double Szy = A.m12;
        double Szz = A.m22;

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

        // Newton-Raphson
        double mxEigenV = E0;
        double oldg = 0.0;
        for (int i = 0; i < maxIter; ++i) {
            oldg = mxEigenV;
            double x2 = mxEigenV * mxEigenV;
            double b = (x2 + C2) * mxEigenV;
            double a = b + C1;
            double delta = ((a * mxEigenV + C0) / (2.0 * x2 * mxEigenV + b + a));
            mxEigenV -= delta;
            System.out.printf("diff[%3d]: %16g %16g %16g\n", i, mxEigenV - oldg, evalprec * mxEigenV, mxEigenV);
            if (Math.abs(mxEigenV - oldg) < Math.abs(evalprec * mxEigenV)) {
                break;
            }
            if (i == maxIter - 1) {
                System.err.printf("More than %d iterations needed!\n", maxIter);
            }
        }

        // The abs() is to guard against extremely small but *negative* numbers due to floating point error
        double rms = Math.sqrt(Math.abs(2.0 * (E0 - mxEigenV) / len));
        rmsdOut[0] = rms;
        /* System.out.printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len); */

        if (minScore > 0) {
            if (rms < minScore) {
                rotOut.identity();
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

        // The recovered quaternion maps frag_1 to frag_2, so need to invert quaternion to map frag_2 onto frag_1
        // (i.e. use -q1 rather than q1 for w)
        rotOut.set(q2, q3, q4, -q1);

        return 1;
    }

    /**
     * Calculate centroid of points. If you are doing a full superposition (the usual least squares way), you MUST
     * calculate the centroid of each structure first.
     * 
     * @param coords       the coordinates to find the centroid for.
     * @param idxsToUse    the indices in the coords and weight arrays to use.
     * @param numIdxsToUse the number of entries in idxsToUse.
     * @param weight       the weight array of size len, or null to use unweighted coordinates.
     * @param centroidOut  out parameter: the centroid
     */
    private static void calcCentroid(Vector3d[] coords, int[] idxsToUse, int numIdxsToUse, double[] weight,
            Vector3d centroidOut) {
        centroidOut.zero();
        if (weight != null) {
            double wsum = 0.0;
            for (int i = 0; i < numIdxsToUse; ++i) {
                int idx = idxsToUse[i];

                centroidOut.x += weight[idx] * coords[idx].x;
                centroidOut.y += weight[idx] * coords[idx].y;
                centroidOut.z += weight[idx] * coords[idx].z;

                wsum += weight[idx];
            }

            if (Math.abs(wsum) < 1.0e-30) {
                throw new IllegalArgumentException("Weights sum to zero");
            }

            centroidOut.x /= wsum;
            centroidOut.y /= wsum;
            centroidOut.z /= wsum;
        } else {
            for (int i = 0; i < numIdxsToUse; ++i) {
                int idx = idxsToUse[i];

                centroidOut.x += coords[idx].x;
                centroidOut.y += coords[idx].y;
                centroidOut.z += coords[idx].z;
            }

            if (numIdxsToUse == 0) {
                throw new IllegalArgumentException("len == 0");
            }

            centroidOut.x /= numIdxsToUse;
            centroidOut.y /= numIdxsToUse;
            centroidOut.z /= numIdxsToUse;
        }
    }

    /**
     * Calculate the RMSD & rotational matrix. Superposition coords2 onto coords1 -- in other words, coords2 is
     * rotated, coords1 is held fixed.
     * 
     * <b>Important:</b> this rotation is calculated about the <i>centroid</> of coords2, so before applying the
     * rotation matrix, you must subtract the centroid from points in coords2, then add the centroid back again
     * after rotation. There will also be a translational offset of (centroid2 - centroid1) required to map coords2
     * to coords1, since the centroid is also subtracted from coords1 before registration.
     * 
     * @param coords1      reference structure.
     * @param coords2      candidate structure.
     * @param idxsToUse    the indices in the coords1, coords2 and weight arrays to use.
     * @param numIdxsToUse the number of entries in idxsToUse.
     * @param weight       the weight array of size len, or null to use unweighted coordinates.
     * 
     * @param rotOut       out parameter: the rotation
     * @param centroid1Out out parameter: the centroid of coords1.
     * @param centroid2Out out parameter: the centroid of coords2.
     * 
     * @return RMSD value
     */
    public static double calcRMSDRotation(Vector3d[] coords1, Vector3d[] coords2, int[] idxsToUse, int numIdxsToUse,
            double[] weight, Quaterniond rotOut, Vector3d centroid1Out, Vector3d centroid2Out) {
        // Find centroid of the structures
        calcCentroid(coords1, idxsToUse, numIdxsToUse, weight, centroid1Out);
        calcCentroid(coords2, idxsToUse, numIdxsToUse, weight, centroid2Out);

        double wsum;
        if (weight == null) {
            wsum = numIdxsToUse;
        } else {
            wsum = 0.0;
            for (int i = 0; i < numIdxsToUse; ++i) {
                int idx = idxsToUse[i];
                wsum += weight[idx];
            }
        }

        // Calculate the (weighted) inner product of two structures
        Matrix3d A = new Matrix3d();
        double E0 = innerProduct(coords1, centroid1Out, coords2, centroid2Out, idxsToUse, numIdxsToUse, weight, A);

        // Calculate the RMSD & rotation
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
        Vector3d[] frag_1 = new Vector3d[len];
        Vector3d[] frag_2 = new Vector3d[len];

        frag_1[0] = new Vector3d(-2.803, -15.373, 24.556);
        frag_1[1] = new Vector3d(0.893, -16.062, 25.147);
        frag_1[2] = new Vector3d(1.368, -12.371, 25.885);
        frag_1[3] = new Vector3d(-1.651, -12.153, 28.177);
        frag_1[4] = new Vector3d(-0.440, -15.218, 30.068);
        frag_1[5] = new Vector3d(2.551, -13.273, 31.372);
        frag_1[6] = new Vector3d(0.105, -11.330, 33.567);

        frag_2[0] = new Vector3d(-14.739, -18.673, 15.040);
        frag_2[1] = new Vector3d(-12.473, -15.810, 16.074);
        frag_2[2] = new Vector3d(-14.802, -13.307, 14.408);
        frag_2[3] = new Vector3d(-17.782, -14.852, 16.171);
        frag_2[4] = new Vector3d(-16.124, -14.617, 19.584);
        frag_2[5] = new Vector3d(-15.029, -11.037, 18.902);
        frag_2[6] = new Vector3d(-18.577, -10.001, 17.996);

        double[] weight = new double[len];
        int[] idxsToUse = new int[len];
        for (int i = 0; i < len; ++i) {
            weight[i] = i + 1.0;
            idxsToUse[i] = i;
        }

        System.out.printf("\nCoords before centering:\n");

        printCoords(frag_1, len);
        printCoords(frag_2, len);

        Quaterniond rot = new Quaterniond();
        Vector3d centroid1 = new Vector3d();
        Vector3d centroid2 = new Vector3d();
        double rmsd = calcRMSDRotation(frag_1, frag_2, idxsToUse, len, weight, rot, centroid1, centroid2);

        System.out.printf("\nCoords after centering:\n");

        printCoords(frag_1, len);
        printCoords(frag_2, len);

        System.out.printf("\nQCP rmsd: %f\n", rmsd);

        System.out.printf("\nQCP Rotation matrix:\n");
        Mat3Print(rot.get(new Matrix3d()));

        System.out.printf("\nQCP quaternion:\n" + rot.w + "\t" + rot.x + "\t" + rot.y + "\t" + rot.z);

        // Apply rotation matrix
        for (int i = 0; i < len; ++i) {
            // Modify frag_2 in place to register it to frag1 
            rot.transform(frag_2[i].sub(centroid2)).add(centroid1);
        }

        /* calculate euclidean distance */
        double euc_dist = 0.0;
        double wtsum = 0.0;
        for (int i = 0; i < len; ++i) {
            wtsum += weight[i];
            double tmp = frag_1[i].distanceSquared(frag_2[i]);
            euc_dist += weight[i] * tmp;
        }

        System.out.printf("\nCoords 2 after rotation:\n");
        printCoords(frag_2, len);

        System.out.printf("\nExplicit RMSD calculated from transformed coords: %f\n\n",
                Math.sqrt(euc_dist / wtsum));
    }
}
