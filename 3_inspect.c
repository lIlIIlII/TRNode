/* 3_inspect.c
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

#include "0_header.h"

#define thr_second 10
#define thr_stable 1e-13

// judge whether system is stable
static int stable (_complex y[], ode_struct* ode) {
    size_t i, j, para;
    double tmp[cst_node];
    register _complex reg;
    _complex factor[cst_node];

    para = ode->para;
    double (*n)[cst_node] = global.n[para];
    double (*kn)[cst_node] = global.kn[para];
    int (*net_top)[cst_node] = ode->net_top;
    memcpy (factor, ode->cst_top, cst_node * sizeof (_complex));

    for (i = 0; i < cst_node; ++i)
        for (j = 0; j < cst_node; ++j)
            if (net_top[i][j]) {
                reg = cpow (y[i], n[i][j]) * kn[i][j];    // (y/k)^n
                factor[j] *= net_top[i][j] - 1 ? 1 : reg; // node_i activate node_j
                factor[j] *= 1.0 / ++reg;                 // (y/k)^n/(1+(y/k)^n)
            }

    tmp[0] = creal ((factor[0] - y[0]) * (factor[0] - y[0]));
    tmp[1] = creal ((factor[1] - y[1]) * (factor[1] - y[1]));
    tmp[2] = creal ((factor[2] - y[2]) * (factor[2] - y[2]));

    return tmp[0] + tmp[1] + tmp[2] < thr_stable;
}

// calculate adaptation parameters of TRN and record frame information
double inspect (double* result, frame_struct* frame_head, ode_struct* ode) {
    int flag_record = ode->flag_record;
    int flag_outlier;
    int flag_stable;

    double t = 0;
    double y_max = 0;
    double y_min = 10;
    double time_consume;
    double con_eps[2] = { global.cst_cal[5], global.cst_cal[6] };

    time_t t_start;
    _complex* py = ode->y;
    frame_struct* frame = frame_head;

    double* response = &result[0];
    double* peakfold = &result[1];
    double* error = &result[2];
    double* ypeak = &result[3];
    double* ystable[2] = { &result[4], &result[7] };

    { // phase 1
        double h = 1e-6;
        ode->cst_top[0] = global.cst_cal[0];
        equation (ode->dydt, py, ode);

        t_start = time (NULL);
        do {
            step (ode, &h, &t);
            time_consume = difftime (time (NULL), t_start);
            flag_outlier = creal (py[0]) + con_eps[0] < 0 || creal (py[0]) - con_eps[1] > 0
                           || creal (py[1]) + con_eps[0] < 0 || creal (py[1]) - con_eps[1] > 0
                           || creal (py[2]) + con_eps[0] < 0 || creal (py[2]) - con_eps[1] > 0;
            flag_stable = stable (py, ode);

            if (flag_record) {
                frame->next = ( frame_struct* )malloc (sizeof (frame_struct));
                frame = frame->next;
                memcpy (frame->dydt, ode->dydt, cst_node * sizeof (_complex));
                memcpy (frame->y, ode->y, cst_node * sizeof (_complex));
                frame->time_consume = time_consume;
                frame->flag_outlier = flag_outlier;
                frame->flag_stable = flag_stable;
                frame->next = NULL;
                frame->phase = 1;
                frame->t = t;
            }

        } while (!flag_stable && !flag_outlier && time_consume < thr_second);

        if (time_consume >= thr_second)
            return *response = mark_timeout;
        else if (flag_outlier)
            return *response = mark_outlier;
        ystable[0][0] = py[0];
        ystable[0][1] = py[1];
        ystable[0][2] = py[2];
    }

    { // phase 2
        double h = 1e-6;
        ode->cst_top[0] = global.cst_cal[1];
        equation (ode->dydt, py, ode);

        t_start = time (NULL);
        do {
            step (ode, &h, &t);
            time_consume = difftime (time (NULL), t_start);
            y_min = creal (py[2]) < y_min ? creal (py[2]) : y_min;
            y_max = creal (py[2]) > y_max ? creal (py[2]) : y_max;
            flag_outlier = creal (py[0]) + con_eps[0] < 0 || creal (py[0]) - con_eps[1] > 0
                           || creal (py[1]) + con_eps[0] < 0 || creal (py[1]) - con_eps[1] > 0
                           || creal (py[2]) + con_eps[0] < 0 || creal (py[2]) - con_eps[1] > 0;
            flag_stable = stable (py, ode);

            if (flag_record) {
                frame->next = ( frame_struct* )malloc (sizeof (frame_struct));
                frame = frame->next;
                memcpy (frame->dydt, ode->dydt, cst_node * sizeof (_complex));
                memcpy (frame->y, ode->y, cst_node * sizeof (_complex));
                frame->time_consume = time_consume;
                frame->flag_outlier = flag_outlier;
                frame->flag_stable = flag_stable;
                frame->next = NULL;
                frame->phase = 2;
                frame->t = t;
            }

        } while (!flag_stable && !flag_outlier && time_consume < thr_second);

        if (time_consume >= thr_second)
            return *response = mark_timeout;
        else if (flag_outlier)
            return *response = mark_outlier;
        ystable[1][0] = py[0];
        ystable[1][1] = py[1];
        ystable[1][2] = py[2];
    }

    *error = fabs (ystable[1][2] - ystable[0][2]) / fabs (ystable[0][2]);
    *peakfold = fabs (y_min - ystable[1][2]) > fabs (y_max - ystable[1][2])
                    ? ystable[0][2] / (*ypeak = y_min)
                    : (*ypeak = y_max) / ystable[0][2];
    *response = fabs (*ypeak - ystable[0][2]);

    if ((*response > 0.1 || *peakfold > 2) && *error < 1e-2 && ystable[0][2] > 1e-3 && ystable[1][2] > 1e-3)
        if (*response >= 1)
            return *response = mark_outlier; // singularity result
        else
            return *response;
    else
        return *response = mark_ntadapt;
}