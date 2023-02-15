/* 4_ode.c
 *
 * Copyright (C) 2023 Song Changli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Original Author: Brian Gough
 */

#include "0_header.h"

// 4th-order Rungekuta algorithm, which calculate a new step for special h
static void rk4 (_complex y[], _complex k[], double h, ode_struct* ode) {
    _complex y0[cst_node];
    _complex temp_y[cst_node];
    memcpy (y0, y, cst_node * sizeof (_complex));
    double c[2] = { global.cst_cal[2], global.cst_cal[3] };
    double h_half = h * 0.5;
    double hc0 = h * c[0];
    double hc1 = h * c[1];

    y[0] += hc0 * k[0];
    y[1] += hc0 * k[1];
    y[2] += hc0 * k[2]; // y1

    temp_y[0] = y0[0] + h_half * k[0];
    temp_y[1] = y0[1] + h_half * k[1];
    temp_y[2] = y0[2] + h_half * k[2];
    equation (k, temp_y, ode); // k2

    y[0] += hc1 * k[0];
    y[1] += hc1 * k[1];
    y[2] += hc1 * k[2]; // y2

    temp_y[0] = y0[0] + h_half * k[0];
    temp_y[1] = y0[1] + h_half * k[1];
    temp_y[2] = y0[2] + h_half * k[2];
    equation (k, temp_y, ode); // k3

    y[0] += hc1 * k[0];
    y[1] += hc1 * k[1];
    y[2] += hc1 * k[2]; // y3

    temp_y[0] = y0[0] + h * k[0];
    temp_y[1] = y0[1] + h * k[1];
    temp_y[2] = y0[2] + h * k[2];
    equation (k, temp_y, ode); // k4

    y[0] += hc0 * k[0];
    y[1] += hc0 * k[1];
    y[2] += hc0 * k[2]; // y4
}

// repeatly try a new step until satisfy error control condition
void step (ode_struct* ode, double* h, double* t) {
    _complex y_1step[cst_node];
    _complex y_2step[cst_node];
    _complex y_err[cst_node];
    _complex k[cst_node];
    _complex* dydt = ode->dydt;
    _complex* y0 = ode->y;
    double r[cst_node];
    double half_h;
    double h0 = *h;
    double h_old;
    double temp;
    double rmax;

    for (;;) {
        { // compare one whole step with two steps that devided into half
            memcpy (y_1step, y0, cst_node * sizeof (_complex));
            memcpy (k, dydt, cst_node * sizeof (_complex));
            rk4 (y_1step, k, h0, ode);

            half_h = 0.5 * h0;
            memcpy (y_2step, y0, cst_node * sizeof (_complex));
            memcpy (k, dydt, cst_node * sizeof (_complex));
            rk4 (y_2step, k, half_h, ode);
            equation (k, y_2step, ode);
            rk4 (y_2step, k, half_h, ode);

            equation (k, y_2step, ode);
            y_err[0] = global.cst_cal[4] * (y_2step[0] - y_1step[0]);
            y_err[1] = global.cst_cal[4] * (y_2step[1] - y_1step[1]);
            y_err[2] = global.cst_cal[4] * (y_2step[2] - y_1step[2]);
        }

        { // adjust h
            rmax = 0;
            h_old = h0;
            r[0] = fabs (y_err[0] / (eps_abs + eps_rel * y_2step[0]));
            r[1] = fabs (y_err[1] / (eps_abs + eps_rel * y_2step[1]));
            r[2] = fabs (y_err[2] / (eps_abs + eps_rel * y_2step[2]));
            rmax = r[0] > r[1] ? (r[0] > r[2] ? r[0] : r[2])
                               : (r[1] > r[2] ? r[1] : r[2]);
            if (rmax > 1.1) {
                temp = 0.9 * pow (rmax, -0.25);
                if (temp < 0.2)
                    temp = 0.2;
                h0 = temp * h_old;
                continue; // try with an updated h0
            } else if (rmax < 0.5) {
                temp = 0.9 * pow (rmax, -0.2);
                if (temp > 4.9)
                    temp = 4.9;
                else if (temp < 1.0)
                    temp = 1.0;
                h0 = temp * h_old; // increase the step with h0 next time and quit
            }

            *t += h_old;
            *h = h0;

            memcpy (y0, y_2step, cst_node * sizeof (_complex));
            memcpy (dydt, k, cst_node * sizeof (_complex));
            break;
        }
    }
}