#include "pybind11/pybind11.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <iostream>
#include <numeric>
#include <cmath>
#include <utility>
#include <iostream>
#include <tuple>

#define PI 3.14159265358979323846264338327950288

namespace py = pybind11;

using VectorD = xt::pytensor<double, 1>;
using VectorI = xt::pytensor<int, 1>;
using MatrixD = xt::pytensor<double, 2>;

float absolute(float x)
{
    if (x > 0.0)
    { return x; }
    else
    { return -x; }
}

int get_detector_coordinates(VectorD& base_point, VectorD& horizontal, VectorD& vertical, VectorD& normal, 
        VectorD& h_out, VectorD& v_out, VectorD& wl_out, VectorD& distance_out, xt::pytensor<bool, 1>& hits_it, 
        const MatrixD& beam, int array_size, float n_dot_base, int height, int width, float wl_min, float wl_max) {

        float az, elev;
        float bx,by,bz;
        float x,y,z, temp;
        float h,v,d;
        float diff_x, diff_y, diff_z;

        //some vars
        float base_point_x = base_point(0);
        float base_point_y = base_point(1);
        float base_point_z = base_point(2);
        float horizontal_x = horizontal(0);
        float horizontal_y = horizontal(1);
        float horizontal_z = horizontal(2);
        float vertical_x = vertical(0);
        float vertical_y = vertical(1);
        float vertical_z = vertical(2);
        float nx = normal(0);
        float ny = normal(1);
        float nz = normal(2);
        float n_dot_base_f = n_dot_base;

        int i;
        int error_count = 0;
        int bad_beam = 0;
        float projection, beam_length,  wavelength;

        for (i=0; i<array_size; i++)
        {
            //Good beam, nice beam.
            bad_beam = 0;

            // Non-normalized beam direction
            bx=beam(0, i);
            by=beam(1, i);
            bz=beam(2, i);

            // So we normalize it
            beam_length = sqrt(bx*bx + by*by + bz*bz);
            bx = bx/beam_length;
            by = by/beam_length;
            bz = bz/beam_length;

            //Check if the wavelength is within range
            wavelength = 6.2831853071795862/beam_length;

            //If there are any nan's in the beam direction, this next check will return false.
            if ((wavelength <= wl_max) && (wavelength >= wl_min))
            {
                //Wavelength in range! Keep going.

                //Make sure the beam points in the same direction as the detector, not opposite to it
                // project beam onto detector's base_point
                projection = (base_point_x*bx)+(base_point_y*by)+(base_point_z*bz);
                if (projection > 0)
                {
                    //beam points towards the detector

                    //This beam coincides with the origin (0,0,0)
                    //Therefore the line equation is x/bx = y/by = z/bz

                    //Now we look for the intersection between the plane of normal nx,ny,nz and the given angle.

                    //Threshold to avoid divide-by-zero
                    float min = 1e-6;
                    if ((absolute(bz) > min)) // && (absolute(nz) > min))
                    {
                        z = n_dot_base_f / ((nx*bx)/bz + (ny*by)/bz + nz);
                        temp = (z / bz);
                        y = by * temp;
                        x = bx * temp;
                    }
                    else if ((absolute(by) > min)) //  && (absolute(ny) > min))
                    {
                        y = n_dot_base_f / (nx*bx/by + ny + nz*bz/by);
                        temp = (y / by);
                        x = bx * temp;
                        z = bz * temp;
                    }
                    else if ((absolute(bx) > min)) //  && (absolute(nx) > min))
                    {
                        x = n_dot_base_f / (nx + ny*by/bx + nz*bz/bx);
                        temp = (x / bx);
                        y = by * temp;
                        z = bz * temp;
                    }
                    else
                    {
                        // The scattered beam is 0,0,0
                        error_count += 1;
                        bad_beam = 1;
                    }
                }
                else
                {
                    //The projection is <0
                    // means the beam is opposite the detector. BAD BEAM! No cookie!
                    bad_beam = 1;
                }
            }
            else
            {
                //Wavelength is out of range. Can't measure it!
                bad_beam = 1;
            }


            if (bad_beam)
            {
                //A bad beam means it does not hit, for sure.
                h_out(i) = NAN;
                v_out(i) = NAN;
                wl_out(i) = wavelength; //This may be NAN too, for NAN inputs.
                hits_it(i) = 0;
            }
            else
            {
                //Valid beam calculation
                //Difference between this point and the base point (the center)
                diff_x = x - base_point_x;
                diff_y = y - base_point_y;
                diff_z = z - base_point_z;

                //Project onto horizontal and vertical axes by doing a dot product
                h = diff_x*horizontal_x + diff_y*horizontal_y + diff_z*horizontal_z;
                v = diff_x*vertical_x + diff_y*vertical_y + diff_z*vertical_z;

                // Save to matrix
                h_out(i) = h;
                v_out(i) = v;

                // the scattered beam is 1/wl long.
                wl_out(i) = wavelength;

                //What was the distance to the detector spot?
                distance_out(i) = sqrt(x*x + y*y + z*z);

                // And do we hit that detector?
                // Detector is square and our h,v coordinates are relative to the center of it.
                hits_it(i) = (v > -height/2) && (v < height/2) && (h > -width/2) && (h < width/2);
            }
        }
        return error_count;
}

MatrixD make_volume_symmetry_map(const MatrixD& B, const MatrixD& invB, MatrixD& symm, double qres, double qlim, int n, int order, xt::pytensor<double, 3>& table) {
    //-- Calculate  the hkl array ---
    int ix, iy, iz;
    int eix, eiy, eiz, eindex;
    int index, ord;
    double qx, qy, qz;
    double eqx, eqy, eqz;
    double h, k, l;
    double eh, ek, el;
    for (ix=0; ix<n; ix++)
    {
        qx = ix*qres - qlim;
        for (iy=0; iy<n; iy++)
        {
            qy = iy*qres - qlim;
            for (iz=0; iz<n; iz++)
            {
                qz = iz*qres - qlim;
                index = iz + iy*n + ix*n*n;
                //Ok, now we matrix multiply invB.hkl to get all the HKLs as a column array
                h = qx * invB(0,0) + qy * invB(0,1) + qz * invB(0,2);
                k = qx * invB(1,0) + qy * invB(1,1) + qz * invB(1,2);
                l = qx * invB(2,0) + qy * invB(2,1) + qz * invB(2,2);

                //Now go through each equivalency table.
                for (ord=0; ord<order; ord++)
                {
                    //Do TABLE.hkl to find a new equivalent hkl
                    eh = h * table(ord, 0,0) + k * table(ord, 0,1) + l * table(ord, 0,2);
                    ek = h * table(ord, 1,0) + k * table(ord, 1,1) + l * table(ord, 1,2);
                    el = h * table(ord, 2,0) + k * table(ord, 2,1) + l * table(ord, 2,2);
                    //Now, matrix mult B . equiv_hkl to get the other q vector
                    eqx = eh * B(0,0) + ek * B(0,1) + el * B(0,2);
                    eqy = eh * B(1,0) + ek * B(1,1) + el * B(1,2);
                    eqz = eh * B(2,0) + ek * B(2,1) + el * B(2,2);

                    //Ok, now you have to find the index into QSPACE
                    eix = round( (eqx+qlim)/qres ); if ((eix >= n) || (eix < 0)) eix = -1;
                    eiy = round( (eqy+qlim)/qres ); if ((eiy >= n) || (eiy < 0)) eiy = -1;
                    eiz = round( (eqz+qlim)/qres ); if ((eiz >= n) || (eiz < 0)) eiz = -1;

                    if ((eix < 0) || (eiy < 0) || (eiz < 0))
                    {
                        //One of the indices was out of bounds.
                        //Put this marker to mean NO EQUIVALENT
                        symm(index, ord) = -1;
                    }
                    else
                    {
                        //No problem!, Now I put it in there
                        eindex = eiz + eiy*n + eix*n*n;
                        //This pixel (index) has this equivalent pixel index (eindex) for this order transform ord.
                        symm(index, ord) = eindex;
                    }

                }

            }
        }
    }
    return symm;
}

VectorD getq(double wl_output, double az, double elevation, double pi, const MatrixD& rot_matrix)
{
    double wl_input = wl_output;
    double r2, x, y, z;

    // The scattered beam emanates from the centre of this spher.
    // Find the intersection of the scattered beam and the sphere, in XYZ
    // We start with an Ewald sphere of radius 1/wavelength

    //Assuming azimuth of zero points to z positive = same direction as incident radiation.
    r2 = cos(elevation) / wl_output;
    z=cos(az) * r2;
    x=sin(az) * r2;

    // Assuming elevation angle is 0 when horizontal, positive to y positive:
    y=sin(elevation) / wl_output;

    //And here is the incident beam direction: Along the z-axis
    double incident_z = 1.0 / wl_input;

    //The vector difference between the two is the q vector
    float qx, qy, qz;
    qx = 2 * pi * x;
    qy = 2 * pi * y;
    qz = 2 * pi * (z - incident_z);

    //#Now we switch to the coordinate system of the crystal.
    //#The scattered beam direction (the detector location) is rotated relative to the crystal because the sample is rotated.
    //#So is the incident beam direction.
    //#Therefore, the q-vector measured is simply rotated

    //Here we perform the rotation by doing a matrix multiplication.
    VectorD q = {
        qx * rot_matrix(0,0) + qy * rot_matrix(0,1) + qz * rot_matrix(0,2),
        qx * rot_matrix(1,0) + qy * rot_matrix(1,1) + qz * rot_matrix(1,2),
        qx * rot_matrix(2,0) + qy * rot_matrix(2,1) + qz * rot_matrix(2,2)
    };
    return q; 
}


std::pair<VectorD, VectorD> getq_inelastic(double wl_input, double wl_output, double az, double elevation, double pi, const MatrixD& rot_matrix) {
    double r2, x, y, z;
    //Assuming azimuth of zero points to z positive = same direction as incident radiation.
    r2 = cos(elevation) / wl_output;
    z=cos(az) * r2;
    x=sin(az) * r2;

    // Assuming elevation angle is 0 when horizontal, positive to y positive:
    y=sin(elevation) / wl_output;

    //And here is the incident beam direction: Along the z-axis
    double incident_z = 1.0 / wl_input;

    //The vector difference between the two is the q vector
    float qx, qy, qz;
    qx = 2 * pi * x;
    qy = 2 * pi * y;
    qz = 2 * pi * (z - incident_z);

    //Here we perform the rotation by doing a matrix multiplication.
    VectorD q = {
        qx * rot_matrix(0,0) + qy * rot_matrix(0,1) + qz * rot_matrix(0,2),
        qx * rot_matrix(1,0) + qy * rot_matrix(1,1) + qz * rot_matrix(1,2),
        qx * rot_matrix(2,0) + qy * rot_matrix(2,1) + qz * rot_matrix(2,2)
    };

    VectorD qlab = { qx, qy, qz };
    return std::make_pair(q, qlab);
}

VectorD apply_volume_symmetry(const VectorD& old_q, VectorD& qspace_flat, int numpix, int order, MatrixD& symm) {
    int pix, ord, index;
    for (pix=0; pix<numpix; pix++)
    {
        //Go through each pixel
        for (ord=0; ord<order; ord++)
        {
            //Now go through each equivalent q.
            index = symm(pix, ord);
            if (index >= 0)
            {
                //Valid index.
                qspace_flat(pix) += old_q(index);
            }
        }
    }
    return qspace_flat;
}

std::tuple<VectorI, VectorI, VectorI, VectorI, VectorI, VectorI> calculate_coverage_stats(const VectorD& qspace, const VectorD& qspace_radius, double q_step, double qlim, VectorI& total_points, int qspace_size, int num, VectorI& covered_points0, VectorI& covered_points1, VectorI& covered_points2, VectorI& covered_points3) {
    int i, j;
    int slice;
    int val;
    int overall_points = 0;
    int overall_covered_points = 0;
    int overall_redundant_points = 0;

    for (i=0; i<qspace_size; i++)
    {
        //Coverage value at this points
        val = qspace(i);

        //Do the overall stats
        if (qspace_radius(i) < qlim)
        {
            //But only within the sphere
            overall_points++;
            if (val > 0)
            {
                overall_covered_points++;
                if (val > 1)
                {
                    overall_redundant_points++;
                }
            }
        }

        //Which slice are we looking at?
        slice = static_cast<int>(qspace_radius(i) / q_step);
        if ((slice < num) && (slice >= 0))
        {
            total_points(slice)++;
            if (val>0)
            {
                covered_points0(slice)++;
                if (val>1)
                {
                    covered_points1(slice)++;
                    if (val>2)
                    {
                        covered_points2(slice)++;
                        if (val>3)
                        {
                            covered_points3(slice)++;
                        }
                    }
                }
            }
        }
    }

    VectorI stats = { overall_points, overall_covered_points, overall_redundant_points };
    return std::make_tuple(stats, total_points, covered_points0, covered_points1, covered_points2, covered_points3);
}

float fitness_function_generic(float phi, float chi, float omega, const VectorD& params)
{
    float phi_min = params(0);
    float phi_max = params(1);
    float chi_min = params(2);
    float chi_max = params(3);
    float omega_min = params(4);
    float omega_max = params(5);

    float phi_mid = (phi_min + phi_max) / 2;
    float chi_mid = (chi_min + chi_max) / 2;
    float omega_mid = (omega_min + omega_max) / 2;

    float fitness =  absolute(chi - chi_mid) + absolute(omega - omega_mid) + absolute(phi - phi_mid);

    // Big penalties for being out of the range
    if (phi < phi_min) fitness += (phi_min - phi) * 1.0;
    if (phi > phi_max) fitness += (phi - phi_max) * 1.0;
    if (chi < chi_min) fitness += (chi_min - chi) * 1.0;
    if (chi > chi_max) fitness += (chi - chi_max) * 1.0;
    if (omega < omega_min) fitness += (omega_min - omega) * 1.0;
    if (omega > omega_max) fitness += (omega - omega_max) * 1.0;

    // if (phi < phi_min || phi > phi_max) fitness += 10;
    // if (chi < chi_min || chi > chi_max) fitness += 10;
    // if (omega < omega_min || omega > omega_max) fitness += 10;

    return fitness;
}

float fitness_function_MANDI(float phi, float chi, float omega, const VectorD& params)
{
    float fitness = absolute(phi);
    return fitness;
}

float fitness_function_SXD(float phi, float chi, float omega, const VectorD& params)
{
    float fitness = absolute(phi);
    return fitness;
}

float fitness_function_SNAP(float phi, float chi, float omega, const VectorD& params)
{
    float phi_min = params(0);
    float phi_max = params(1);
    float chi_mid = params(2);
    float phi_mid = (phi_min + phi_max) / 2;

    float fitness = absolute(chi - chi_mid)*10.0 + absolute(phi - phi_mid)/10.0;

    // Big penalties for being out of the range
    if (phi < phi_min) fitness += (phi_min - phi) * 1.0;
    if (phi > phi_max) fitness += (phi - phi_max) * 1.0;

    return fitness;
}

float fitness_function_MANDIVaryOmega(float phi, float chi, float omega, const VectorD& params)
{
    float phi_min = params(0);
    float phi_max = params(1);
    float omega_min = params(2);
    float omega_max = params(3);

    float phi_mid = (phi_min + phi_max) / 2;
    float chi_mid = params(4);
    float omega_mid = (omega_min + omega_max) / 2;

    float fitness = absolute(chi - chi_mid)*10.0 + absolute(omega - omega_mid)/10.0 + absolute(phi - phi_mid)/10.0;

    // Big penalties for being out of the range
    if (phi < phi_min) fitness += (phi_min - phi) * 1.0;
    if (phi > phi_max) fitness += (phi - phi_max) * 1.0;
    if (omega < omega_min) fitness += (omega_min - omega) * 1.0;
    if (omega > omega_max) fitness += (omega - omega_max) * 1.0;

    return fitness;
}


float fitness_function_HB3A(float phi, float chi, float omega, const VectorD& params)
{
    double center = 3.14159*25.0/180.0;
    double omegadiff = omega - center;
    if (omegadiff < 0) omegadiff = -omegadiff;

    //if (omegadiff > center)
    // omegadiff = omegadiff + (omega-center) * 10.0;

    return absolute(chi) + omegadiff + absolute(phi)/1000.0;
    //return absolute(chi) + omegadiff + absolute(phi);
}


float fitness_function(float phi, float chi, float omega, const std::string& func_name, const VectorD& params) {
    if (func_name == "general") {
        return fitness_function_generic(phi, chi, omega, params);
    } else if (func_name == "mandi") {
        return fitness_function_MANDI(phi, chi, omega, params);
    } else if (func_name == "sxd") {
        return fitness_function_SXD(phi, chi, omega, params);
    } else if (func_name == "mandi_vary_omega") {
        return fitness_function_MANDIVaryOmega(phi, chi, omega, params);
    } else if (func_name == "snap") {
        return fitness_function_SNAP(phi, chi, omega, params);
    } else if (func_name == "HB3A") {
        return fitness_function_HB3A(phi, chi, omega, params);
    } else {
        throw std::runtime_error("Unknown fitness function!");
    }
}


int angle_fitness(const VectorD& rot_angle_list, const VectorD& ending_vec, const MatrixD& initial_rotation_matrix, VectorD& fitnesses, VectorD& chi_list, VectorD& phi_list, VectorD& omega_list, const std::string& func_name, const VectorD& params) {
    float rot_angle;
    int angle_num;
    int output_index = 0;
    for (angle_num=0;  angle_num < rot_angle_list.size(); angle_num++)
    {
        rot_angle = rot_angle_list(angle_num);
        //printf("angle of %e\\n", rot_angle);
        // --- Make the rotation matrix around the ending_vec ----
        float c = cos(rot_angle);
        float s = sin(rot_angle);
        float x,y,z;
        x = ending_vec(0);
        y = ending_vec(1);
        z = ending_vec(2);

        float extra_rotation_matrix[3][3] = {
        {1 + (1-c)*(x*x-1), -z*s+(1-c)*x*y, y*s+(1-c)*x*z},
        {z*s+(1-c)*x*y, 1 + (1-c)*(y*y-1), -x*s+(1-c)*y*z},
        {-y*s+(1-c)*x*z,  x*s+(1-c)*y*z,  1 + (1-c)*(z*z-1)}
        };

        // Do matrix multiplication
        float total_rot_matrix[3][3];

        int i,j,k;
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
            {
                total_rot_matrix[i][j] = 0;
                for (k=0; k<3; k++)
                {
                    total_rot_matrix[i][j] += extra_rotation_matrix[i][k] * initial_rotation_matrix(k,j);
                }
                // printf("%f, ", total_rot_matrix[i][j]);
            }

        //-------- Now we find angles_from_rotation_matrix() -----------
        float chi, phi, omega;

        //#Let's make 3 vectors describing XYZ after rotations
        float ux = total_rot_matrix[0][0];
        float uy = total_rot_matrix[1][0];
        float uz = total_rot_matrix[2][0];
        float vx = total_rot_matrix[0][1];
        float vy = total_rot_matrix[1][1];
        float vz = total_rot_matrix[2][1];
        float nx = total_rot_matrix[0][2];
        float ny = total_rot_matrix[1][2];
        float nz = total_rot_matrix[2][2];

        //#is v.y vertical?
        if (absolute(vy) < 1e-8)
        {
            //#Chi rotation is 0, so we just have a rotation about y
            chi = 0.0;
            phi = atan2(nx, nz);
            omega = 0.0;
        }
        else if (absolute(vy+1) < 1e-8)
        {
            //#Chi rotation is 180 degrees
            chi = PI;
            phi = -atan2(nx, nz);
            if (phi==-PI) phi=PI;
            omega = 0.0;
        }
        else
        {
            //#General case
            phi = atan2(ny, uy);
            chi = acos(vy);
            omega = atan2(vz, -vx);
        }

        float fitness;
        float old_phi = phi;
        float old_chi = chi;
        float old_omega = omega;

        // Try the original angles
        fitness = fitness_function(phi, chi, omega, func_name, params);
        fitnesses(output_index) = fitness;
        phi_list(output_index) = phi;
        chi_list(output_index) = chi;
        omega_list(output_index) = omega;
        output_index++;

        //Make angles closer to 0
        if (phi > PI) phi -= 2*PI;
        if (chi > PI) chi -= 2*PI;
        if (omega > PI) omega -= 2*PI;
        if (phi < -PI) phi += 2*PI;
        if (chi < -PI) chi += 2*PI;
        if (omega < -PI) omega += 2*PI;
        fitness = fitness_function(phi, chi, omega, func_name, params);
        fitnesses(output_index) = fitness;
        phi_list(output_index) = phi;
        chi_list(output_index) = chi;
        omega_list(output_index) = omega;
        output_index++;

        //(phi-pi, -chi, omega-pi) is always equivalent
        phi = old_phi-PI;
        chi = -old_chi;
        omega = old_omega-PI;
        if (phi > PI) phi -= 2*PI;
        if (chi > PI) chi -= 2*PI;
        if (omega > PI) omega -= 2*PI;
        if (phi < -PI) phi += 2*PI;
        if (chi < -PI) chi += 2*PI;
        if (omega < -PI) omega += 2*PI;
        fitness = fitness_function(phi, chi, omega, func_name, params);
        fitnesses(output_index) = fitness;
        phi_list(output_index) = phi;
        chi_list(output_index) = chi;
        omega_list(output_index) = omega;
        output_index++;
    }
}

double vector_length(const VectorD& vector)
{
    double length = 0;
    for (int i=0; i < vector.size(); i++)
        length += vector(i)*vector(i);
    return sqrt(length);
}


void shrink_q_vector(VectorD& q, double limit)
{
    double length = vector_length(q);
    if (length <= 0)
        return;

    for (int i=0; i < q.size(); i++)
    {
        if (length > limit)
        {q(i) = q(i) * (limit / length); }
        else
        {q(i) = q(i); }
    }
    // return q_out;
}

 xt::pytensor<int, 4> calculate_coverage(const VectorI& xlist, const VectorI& ylist, const MatrixD& azimuthal_angle, const MatrixD& elevation_angle, const MatrixD& rot_matrix, int32_t set_value1, int32_t set_value2, int number_of_ints,  xt::pytensor<int, 4>& coverage, int stride, int max_index, double wl_min, double wl_max, double qlim, double q_resolution) {
    const auto pi = PI;
    //Loop through pixels using the list given before.
    for (int iix = 0; iix < xlist.size(); iix++)
    {
        int ix = xlist[iix];
        for (int iiy = 0; iiy < ylist.size(); iiy++)
        {
            int iy = ylist[iiy];

            //Angles of the detector pixel positions
            double az, elev;
            az = azimuthal_angle(iy, ix);
            elev = elevation_angle(iy, ix);

            //Get the two limits.
            auto q_min = getq(wl_min, az, elev, pi, rot_matrix);
            auto q_max = getq(wl_max, az, elev, pi, rot_matrix);
            //Limit them to the modeled size
            shrink_q_vector(q_min, qlim);
            shrink_q_vector(q_max, qlim);

            //Find out how long of a line that is
            double q_length = (vector_length(q_max) - vector_length(q_min));
            if (q_length<0) q_length = -q_length;

            //How many steps will we take. The multiplication factor here is a fudge to make sure it covers fully.
            int numfrac = (1.25 * q_length) / (q_resolution);

            // if (numfrac == 0) printf("numfrac is %d, q_length %f, qmin and max are %f and %f, abs is %f \\n", numfrac, q_length, vector_length(q_min), vector_length(q_max),   (vector_length(q_max) - vector_length(q_min))  );

            if (numfrac > 0)
            {
                //There is something to measure.

                //Size in q-space of each step
                double dx, dy, dz;
                dx = (double(q_max[0]) - double(q_min[0])) / numfrac;
                dy = (double(q_max[1]) - double(q_min[1])) / numfrac;
                dz = (double(q_max[2]) - double(q_min[2])) / numfrac;


                long index;
                double qx, qy, qz;
                long iqx, iqy, iqz;
                // unsigned int* coverage_int = (unsigned int*) coverage;

                double lim_min = -qlim;
                double lim_max = +qlim;

                double q_min_x = double(q_min[0]);
                double q_min_y = double(q_min[1]);
                double q_min_z = double(q_min[2]);
                //printf("Setvalue1 is %d\\n", set_value1);

                //Okay now we draw the line from q_min to q_max.
                long i;
                for (i=0; i<numfrac; i++)
                {
                    //All of these qx checks might not be necessary anymore...?
                    qx = q_min_x + i*dx;
                    iqx = round((qx - lim_min) / q_resolution);
                    if ((iqx >=0) && (iqx < stride))
                    {
                        qy = q_min_y + i * dy;
                        iqy = round((qy - lim_min) / q_resolution);
                        if ((iqy >=0) && (iqy < stride))
                        {
                            qz = q_min_z + i * dz;
                            iqz = round((qz - lim_min) / q_resolution);
                            if ((iqz >=0) && (iqz < stride))
                            {

                                if (number_of_ints==2)
                                {
                                    coverage(iqx,iqy,iqz,0) = 1;
                                    coverage(iqx,iqy,iqz,1) = 1;
                                }
                                else
                                {
                                    coverage(iqx,iqy,iqz,0) = 1;
                                }

                            }
                        }
                    
                } //for i in numfrac
            } //numfrac > 0 so we can draw the line

        } //for iiy

    }
}
 return coverage;
}



void calculate_coverage_inelastic(double wl_input, const VectorI& xlist, const VectorI& ylist, const MatrixD& azimuthal_angle, const MatrixD& elevation_angle, const MatrixD& rot_matrix, double energy_constant, double ki_squared, double ki,  xt::pytensor<bool, 3>& coverage, int stride, int max_index, double wl_min, double wl_max, double qlim, double q_resolution) {
    double kfz, kf_squared, E;

    double lim_min = -qlim;
    double lim_max = +qlim;

    //Loop through pixels using the list given before.
    for (int iix = 0; iix < xlist.size(); iix++)
    {
        int ix = xlist[iix];
        for (int iiy = 0; iiy < ylist.size(); iiy++)
        {
            int iy = ylist[iiy];

            //Angles of the detector pixel positions
            double az, elev;
            az = azimuthal_angle(iy, ix);
            elev = elevation_angle(iy, ix);

            //Get the two limits.
            auto q_max_both = getq_inelastic(wl_input, wl_min, az, elev, PI, rot_matrix);
            auto q_min_both = getq_inelastic(wl_input, wl_max, az, elev, PI, rot_matrix);

            // Rotated qs
            double q_max[3], q_min[3], q_max_unrot[3], q_min_unrot[3];
            double q_max_length = 0;
            double q_min_length = 0;
            for (int i=0; i<3; i++)
            {
                q_max[i] = q_max_both.first(i);
                q_max_unrot[i] = q_max_both.second(i);
                q_max_length += q_max[i]*q_max[i];

                q_min[i] = q_min_both.first(i);
                q_min_unrot[i] = q_min_both.second(i);
                q_min_length += q_min[i]*q_min[i];
            }
            q_max_length = sqrt(q_max_length);
            q_min_length = sqrt(q_min_length);

            /*
                //Limit them to the modeled size -- SHRINK VECTOR ---
                double reduceby_max = 1.0;
                double reduceby_min = 1.0;
                if (q_max_length > qlim)
                    { reduceby_max = q_max_length/qlim; }
                if (q_min_length > qlim)
                    { reduceby_min = q_min_length/qlim; }
                for (int i=0; i<3; i++)
                {
                    q_max[i] /= reduceby_max;
                    q_min[i] /= reduceby_min;
                    q_max_unrot[i] /= reduceby_max;
                    q_min_unrot[i] /= reduceby_min;
                }
            */

            //Vector difference max-min
            double q_diff[3];
            double q_diff_unrot[3];
            double q_length = 0.0;
            for (int i=0; i<3; i++)
            {
                // We change the sign, because for inelastic, Q = ki-kf, rather than the opposite sign for elastic
                q_diff[i] = q_max[i] - q_min[i];
                q_diff_unrot[i] = q_max_unrot[i] - q_min_unrot[i];
                q_length += q_diff[i]*q_diff[i];
            }

            //Find out how long of a line that is
            q_length = sqrt(q_length);

            //How many steps will we take. The multiplication factor here is a fudge to make sure it covers fully.
            long numfrac = (1.25 * q_length) / (q_resolution);

            if (numfrac > 0)
            {
                //There is something to measure.

                //Size in q-space of each step
                double dx, dy, dz;
                dx = q_diff[0] / numfrac;
                dy = q_diff[1] / numfrac;
                dz = q_diff[2] / numfrac;

                double dx_unrot, dy_unrot, dz_unrot;
                dx_unrot = q_diff_unrot[0] / numfrac;
                dy_unrot = q_diff_unrot[1] / numfrac;
                dz_unrot = q_diff_unrot[2] / numfrac;

                long index;
                double qx, qy, qz;
                double qx_unrot, qy_unrot, qz_unrot;
                long iqx, iqy, iqz;

                //Okay now we draw the line from q_min to q_max.
                double i;
                for (i=0; i<numfrac; i++)
                {
                    //All of these qx checks might not be necessary anymore...?
                    qx = q_min[0] + i*dx;
                    qx_unrot = q_min_unrot[0] + i*dx_unrot;
                    iqx = long(round((qx - lim_min) / q_resolution));
                    if ((iqx >=0) && (iqx < stride))
                    {
                        qy = q_min[1] + i*dy;
                        qy_unrot = q_min_unrot[1] + i*dy_unrot;
                        iqy = long(round((qy - lim_min) / q_resolution));
                        if ((iqy >=0) && (iqy < stride))
                        {
                            qz = q_min[2] + i*dz;
                            qz_unrot = q_min_unrot[2] + i*dz_unrot;
                            iqz = long(round((qz - lim_min) / q_resolution));
                            if ((iqz >=0) && (iqz < stride))
                            {
                                //Calculate the neutron energy gain
                                // But we need to get back the z component of kf
                                // We kept Q = kf - ki
                                kfz = (ki + qz_unrot);

                                // Okay, now calculate kf^2
                                kf_squared = qx_unrot*qx_unrot + qy_unrot*qy_unrot + kfz*kfz;

                                // Get the energy. The constant sets the units (to meV)
                                E = energy_constant * (kf_squared - ki_squared);

                                coverage(iqx,iqy,iqz) = E;
                            }
                        }
                    }

                } //for i in numfrac
            } //numfrac > 0 so we can draw the line

        } //for iiy
    }

}

// Python Module and Docstrings

PYBIND11_MODULE(crystal_plan_c_ext, m)
{
    xt::import_numpy();

    m.doc() = R"pbdoc(
        An example xtensor extension
    )pbdoc";

    m.def("getq", getq, "");
    m.def("get_detector_coordinates", get_detector_coordinates, "");
    m.def("make_volume_symmetry_map", make_volume_symmetry_map, "");
    m.def("apply_volume_symmetry", apply_volume_symmetry, "");
    m.def("calculate_coverage_stats", calculate_coverage_stats, "");
    m.def("angle_fitness", angle_fitness, "");
    m.def("calculate_coverage", calculate_coverage, "");
    m.def("calculate_coverage_inelastic", calculate_coverage_inelastic, "");
    m.def("getq", getq, "");
    m.def("getq_inelastic", getq_inelastic, "");
}
