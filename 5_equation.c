/* 5_equation.c
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

// calculate TRN dydt by provided ode infomation
void equation (_complex dydt[], _complex y[], ode_struct* ode) {
    size_t i, j, para;
    register _complex reg;
    _complex tmp, retro[cst_node] = { 0, 0, 0 };

    para = ode->para;
    double* t1 = global.t1[para];
    double* eta_list = ode->eta_list;
    double (*n)[cst_node] = global.n[para];
    double (*n2)[cst_node] = global.n2[para];
    double (*kn)[cst_node] = global.kn[para];
    int (*net_top)[cst_node] = ode->net_top;
    memcpy (dydt, ode->cst_top, cst_node * sizeof (_complex));

    for (i = 0; i < cst_node; ++i)
        for (j = 0; j < cst_node; ++j)
            if (net_top[i][j]) {
                reg = cpow (y[i], n[i][j] - 1);         // y^(n-1)
                tmp = reg * n2[i][j] * kn[i][j];        // n^2*y^(n-1)*k^(-n)
                reg *= y[i] * kn[i][j];                 // (y/k)^n
                dydt[j] *= net_top[i][j] - 1 ? 1 : reg; // node_i activate node_j
                reg = 1.0 / ++reg;                      // 1/(1+(y/k)^n)
                dydt[j] *= reg;                         // (y/k)^n/(1+(y/k)^n)
                retro[i] += tmp * reg * reg;            // n^2*y^(n-1)*k^(-n)/(1+(y/k)^n)^2
            }

    dydt[0] = t1[0] * (dydt[0] - y[0]) / (1 + retro[0] * eta_list[0]);
    dydt[1] = t1[1] * (dydt[1] - y[1]) / (1 + retro[1] * eta_list[1]);
    dydt[2] = t1[2] * (dydt[2] - y[2]) / (1 + retro[2] * eta_list[2]);
}