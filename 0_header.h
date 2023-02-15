/* 0_header.h
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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define eps_abs 1e-7
#define eps_rel 1e-7

#define mark_ntadapt -3
#define mark_timeout -2
#define mark_outlier -1

#define cst_node 3
#define cst_para 10000

typedef double complex _complex;

typedef struct {
    double kn[cst_para][cst_node][cst_node]; // k^(-n)
    double n2[cst_para][cst_node][cst_node]; // n^2
    double n[cst_para][cst_node][cst_node];
    double t1[cst_para][cst_node]; // tau^(-1)

    char path_env[FILENAME_MAX];
    double cst_cal[7]; // constant: input, rk4, eps
    time_t clock;
    int total;
} global_struct;

typedef struct case_struct {
    int case_id;
    int eta_qty;
    int* content;  // single: net, para; loop: net_start, net_end; discret: uid
    int case_type; // 0: single; 1: loop; 2: discrete
    int content_qty;
    int progress[2]; // loop: net_i, para_i; discrete: case_i
    double** eta_all;
    struct case_struct* next;
} case_struct;

typedef struct frame_struct {
    _complex dydt[cst_node];
    _complex y[cst_node];
    double t;
    int phase;
    int flag_stable;
    int flag_outlier;
    int time_consume;
    struct frame_struct* next;
} frame_struct;

typedef struct {
    _complex y[cst_node];
    _complex dydt[cst_node];
    _complex cst_top[cst_node];
    double eta_list[cst_node];
    int net_top[cst_node][cst_node];
    int flag_record;
    size_t para;
} ode_struct;

extern global_struct global;

double inspect (double*, frame_struct*, ode_struct*);
void equation (_complex*, _complex*, ode_struct*);
void step (ode_struct*, double*, double*);
void init (case_struct*, int, int);