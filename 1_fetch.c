/* 1_fetch.c
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
#include <fcntl.h>
#include <sys/stat.h>

global_struct global;

case_struct* CASE;

// check legality of program argument
static int check (int argc, char** argv, char* path) {
    FILE* fp;
    getcwd (global.path_env, FILENAME_MAX);
    char error_msg[FILENAME_MAX] = "\033[1m\033[31merror:\033[0m ";

    if (argc < 2 || argc > 6 || argc == 5) {
        printf ("%sargument incorrect.\n", error_msg);
        return 1;
    }
    sprintf (path, "%s/data/para/const", global.path_env);
    if (access (path, 0) || !(fp = fopen (path, "r"))) {
        printf ("%scan not open const file or path incorrect: %s.\n", error_msg, path);
        return 1;
    }
    fclose (fp);
    sprintf (path, "%s/data/para/LHS", global.path_env);
    if (access (path, 0) || !(fp = fopen (path, "r"))) {
        printf ("%scan not open LHS file or path incorrect: %s.\n", error_msg, path);
        return 1;
    }
    fclose (fp);

    switch (argc) {
    case 2:
        if (argv[1][0] == '/')
            sprintf (path, "%s", argv[1]);
        else
            sprintf (path, "%s/%s", global.path_env, argv[1]);
        if (access (path, 0) || !(fp = fopen (path, "r"))) {
            printf ("%scan not access config file or path incorrect: %s\n", error_msg, path);
            return 1;
        }
        fclose (fp);
        break;
    case 3:
        if (atoi (argv[1]) < 0 || atoi (argv[1]) >= 19683 || atoi (argv[2]) < 0 || atoi (argv[2]) >= 19683 || atoi (argv[1]) > atoi (argv[2])) {
            printf ("%sargument incorrect.\n", error_msg);
            return 1;
        }
        sprintf (path, "%s/data/para/eta", global.path_env);
        if (access (path, 0) || !(fp = fopen (path, "r"))) {
            printf ("%scan not access eta file or path incorrect: %s\n", error_msg, path);
            return 1;
        }
        fclose (fp);
        break;
    default:
        if (atoi (argv[1]) < 0 || atoi (argv[1]) >= 19683 || atoi (argv[2]) < 0 || atoi (argv[2]) >= cst_para) {
            printf ("%sargument incorrect.\n", error_msg);
            return 1;
        }
        break;
    }
    return 0;
}

// fetch 'const' and 'LHS' data then pre-proceed
static int data () {
    int up;
    FILE* fp;
    size_t i, j, k;
    char path[FILENAME_MAX];
    double kl, kh, nl, nh, taul, tauh, base, fold, tmp;

    sprintf (path, "%s/data/para/const", global.path_env);
    fp = fopen (path, "r");
    fscanf (fp, "%lf%lf%lf%lf%lf%lf%lf%lf%d", &kl, &kh, &nl, &nh, &taul, &tauh, &base, &fold, &up);
    fclose (fp);

    sprintf (path, "%s/data/para/LHS", global.path_env);
    fp = fopen (path, "r");
    for (i = 0; i < cst_para; ++i) {
        for (j = 0; j < cst_node; ++j)
            for (k = 0; k < cst_node; ++k) {
                fscanf (fp, "%lf", &tmp);
                global.kn[i][j][k] = kl * pow (10, tmp * log10 (kh / kl));
            }
        for (j = 0; j < cst_node; ++j)
            for (k = 0; k < cst_node; ++k) {
                fscanf (fp, "%lf", &tmp);
                global.n[i][j][k] = nl + tmp * (nh - nl);
                global.n2[i][j][k] = pow (global.n[i][j][k], 2);
                global.kn[i][j][k] = pow (1.0 / global.kn[i][j][k], global.n[i][j][k]);
            }
        for (j = 0; j < cst_node; ++j) {
            fscanf (fp, "%lf", &tmp);
            global.t1[i][j] = 1.0 / (taul * pow (10, tmp * log10 (tauh / taul)));
        }
    }
    fclose (fp);

    global.cst_cal[2] = 1.0 / 6.0;
    global.cst_cal[3] = 1.0 / 3.0;
    global.cst_cal[4] = 4.0 / 15.0;
    global.cst_cal[5] = 2 * eps_abs;
    global.cst_cal[6] = 2 * eps_abs + 1;
    global.cst_cal[0] = base / 0.4 / (1 + base / 0.4);
    base *= up * 10;
    global.cst_cal[1] = base / 0.4 / (1 + base / 0.4);

    return 0;
}

// fetch config file and generate case link list
static int cases (int argc, char** argv, int* eta_max, char* path) {
    case_struct* cur;
    size_t i, j, k;
    int consistent;
    int case_qty;
    double eta;
    FILE* fp;

    CASE = ( case_struct* )malloc (sizeof (case_struct));
    global.total = 0;
    *eta_max = 0;
    cur = CASE;
    switch (argc) {
    case 2: // config
        fp = fopen (path, "r");
        fscanf (fp, "%d", &case_qty);
        for (i = 0; i < case_qty; ++i) {
            cur->case_id = i;
            fscanf (fp, "%d %d", &consistent, &cur->eta_qty);
            *eta_max = cur->eta_qty > *eta_max ? cur->eta_qty : *eta_max;
            cur->eta_all = ( double** )malloc (cur->eta_qty * sizeof (double*));
            for (j = 0; j < cur->eta_qty; ++j) {
                cur->eta_all[j] = ( double* )malloc (cst_node * sizeof (double));
                if (consistent) {
                    fscanf (fp, "%lf", &eta);
                    for (k = 0; k < cst_node; ++k)
                        cur->eta_all[j][k] = eta;
                } else
                    fscanf (fp, "%lf %lf %lf", &cur->eta_all[j][0], &cur->eta_all[j][1], &cur->eta_all[j][2]);
            }
            fscanf (fp, "%d", &cur->content_qty);
            if (cur->content_qty == -1) { // config-loop
                cur->case_type = 1;
                cur->content_qty = 2;
                cur->content = ( int* )malloc (cur->content_qty * sizeof (int));
                fscanf (fp, "%d %d", &cur->content[0], &cur->content[1]);
                cur->progress[0] = cur->content[0];
                cur->progress[1] = 0;
                global.total += cst_para * (cur->content[1] - cur->content[0] + 1);
            } else { // discrete
                cur->case_type = 2;
                cur->content = ( int* )malloc (cur->content_qty * sizeof (int));
                for (j = 0; j < cur->content_qty; ++j)
                    fscanf (fp, "%d", &cur->content[j]);
                cur->progress[0] = 0;
                global.total += cur->content_qty;
            }
            if (i != case_qty - 1) {
                cur->next = ( case_struct* )malloc (sizeof (case_struct));
                cur = cur->next;
            }
            cur->next = NULL;
        }
        fclose (fp);
        return case_qty;
    case 3: // loop
        fp = fopen (path, "r");
        fscanf (fp, "%d %d", &consistent, &cur->eta_qty);
        *eta_max = cur->eta_qty > *eta_max ? cur->eta_qty : *eta_max;
        cur->eta_all = ( double** )malloc (cur->eta_qty * sizeof (double*));
        for (i = 0; i < cur->eta_qty; ++i) {
            cur->eta_all[i] = ( double* )malloc (cst_node * sizeof (double));
            if (consistent) {
                fscanf (fp, "%lf", &eta);
                for (j = 0; j < cst_node; ++j)
                    cur->eta_all[i][j] = eta;
            } else
                fscanf (fp, "%lf %lf %lf", &cur->eta_all[i][0], &cur->eta_all[i][1], &cur->eta_all[i][2]);
        }
        cur->case_id = 0;
        cur->next = NULL;
        cur->case_type = 1;
        cur->content_qty = 2;
        cur->content = ( int* )malloc (2 * sizeof (int));
        cur->content[1] = atoi (argv[2]);
        cur->progress[0] = cur->content[0] = atoi (argv[1]);
        global.total = cst_para * (cur->content[1] - cur->content[0] + 1);
        cur->progress[1] = 0;
        fclose (fp);
        return 1;
    default: // single
        global.total = 1;
        cur->case_id = 0;
        cur->next = NULL;
        cur->case_type = 0;
        cur->content_qty = 2;
        cur->eta_qty = argc - 3;
        *eta_max = cur->eta_qty;
        cur->content = ( int* )malloc (2 * sizeof (int));
        cur->eta_all = ( double** )malloc (cur->eta_qty * sizeof (double*));
        for (i = 0; i < cur->eta_qty; ++i) {
            cur->eta_all[i] = ( double* )malloc (cst_node * sizeof (double));
            for (j = 0; j < cst_node; ++j)
                cur->eta_all[i][j] = atof (argv[argc == 4 ? 3 : 3 + j]);
        }
        cur->content[0] = atoi (argv[1]);
        cur->content[1] = atoi (argv[2]);
        return 0;
    }
}

// create result folders and files base on current time and case quantity
static int files (int case_qty, char** argv) {
    char file[FILENAME_MAX], folder[FILENAME_MAX];
    struct tm* timer;
    size_t i;

    global.clock = time (NULL);
    timer = gmtime (&global.clock);
    sprintf (global.path_env, "%s/data/%d%02d%02d%02d%02d%02d", global.path_env, (1900 + timer->tm_year), (1 + timer->tm_mon), timer->tm_mday, timer->tm_hour, timer->tm_min, timer->tm_sec);
    if (access (global.path_env, 0))
        mkdir (global.path_env, 0777);

    if (case_qty)
        for (i = 0; i < case_qty; ++i) {
            sprintf (folder, "%s/%ld", global.path_env, i);
            if (access (folder, 0))
                mkdir (folder, 0777);
            sprintf (file, "%s/response", folder);
            if (access (file, 0))
                creat (file, 0777);
            sprintf (file, "%s/timeout", folder);
            if (access (file, 0))
                creat (file, 0777);
            sprintf (file, "%s/detail", folder);
            if (access (file, 0))
                creat (file, 0777);
        }
    else {
        sprintf (file, "%s/%09d", global.path_env, atoi (argv[1]) * 10000 + atoi (argv[2]));
        if (access (file, 0))
            creat (file, 777);
    }
    return 0;
}

// main procedure
int main (int argc, char* argv[]) {
    int eta_max;
    char path[FILENAME_MAX];
    if (check (argc, argv, path) || data () || files (cases (argc, argv, &eta_max, path), argv))
        exit (-1);

    init (CASE, eta_max, argc > 3);

    return 0;
}