/* 2_queue.c
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
#include <pthread.h>

#define QUEUE_MAX_LENGTH 10000

typedef struct {
    ode_struct ode;
    frame_struct* frame;
    double** eta_all;
    double** result;
    int case_type;
    int case_id;
    int eta_qty;
    int status; // 0.unsolve 1.solved 2.done
    int uid;
} slot_struct;

static struct {
    pthread_mutex_t lock;
    slot_struct* slot;
    size_t front, rear;
} QUEUE;

static case_struct* CASE;

// generate and judge whether net topology can possibly achieve adaptation
static int topology (int (*net_top)[cst_node], int flag_hill[], int net) {
    size_t i, j;
    int tmp = net * 3;
    for (i = 0; i < cst_node; ++i)
        for (j = 0; j < cst_node; ++j)
            net_top[i][j] = (tmp /= 3) % 3 - 1;
    if (fabs (net_top[0][1]) + fabs (net_top[0][2]) + fabs (net_top[1][0]) + fabs (net_top[1][2]) + fabs (net_top[2][0]) + fabs (net_top[2][1]) <= 1 || !(net_top[0][2] || (net_top[0][1] && net_top[1][2])))
        return 0;
    for (i = 1; i < cst_node; ++i) {
        flag_hill[i - 1] = 0;
        for (j = 0; j < cst_node; ++j)
            flag_hill[i - 1] = net_top[j][i] ? 1 : flag_hill[i - 1];
    }
    return 1;
}

// output completed task to result file
static int depart () {
    char status[4][FILENAME_MAX] = { "outlier", "timeout", "ntadapt", "\033[1m\033[32madapt\033[0m" };
    char path[2][FILENAME_MAX];
    double **datas, *data;
    int flag_output;
    size_t i, j;
    frame_struct *frame, *frame_pre;
    slot_struct* slot = &QUEUE.slot[QUEUE.rear];
    FILE* fp[2];

    if (slot->status != 1)
        return slot->status;

    if (slot->ode.flag_record) {
        data = slot->result[0];
        sprintf (path[0], "%s/%09d", global.path_env, slot->uid);
        fp[0] = fopen (path[0], "w+");
        frame = slot->frame->next;
        slot->frame->next = NULL;
        while (frame) {
            fprintf (fp[0], " %12lf %12lf %12lf %12lf\n", frame->t, creal (frame->y[0]), creal (frame->y[1]), creal (frame->y[2]));
            frame_pre = frame;
            frame = frame->next;
            free (frame_pre);
        }

        printf ("> %s\n", data[0] < 0 ? status[( int )fabs (data[0]) - 1] : status[3]);
        if (data[0] > 0 || data[0] == mark_ntadapt)
            printf ("response\t%12lf\npeakfold\t%12lf\nerror\t\t%12lf\nypeak\t\t%12lf\nystable\t\t%12lf\t%12lf\t%12lf\n\t\t%12lf\t%12lf\t%12lf\n", data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9]);
    } else {
        datas = slot->result;
        if (datas[0][0] == mark_timeout) {
            sprintf (path[0], "%s/%d/timeout", global.path_env, slot->case_id);
            fp[0] = fopen (path[0], "a+");
            fprintf (fp[0], "%09d\n", slot->uid);
            fclose (fp[0]);
        } else {
            flag_output = 0;
            for (i = 0; !flag_output && i < slot->eta_qty; ++i)
                flag_output = datas[i][0] > 0;
            if (flag_output) {
                sprintf (path[0], "%s/%d/response", global.path_env, slot->case_id);
                sprintf (path[1], "%s/%d/detail", global.path_env, slot->case_id);
                fp[0] = fopen (path[0], "a+");
                fp[1] = fopen (path[1], "a+");

                fprintf (fp[0], "%09d", slot->uid);
                for (i = 0; i < slot->eta_qty; ++i)
                    if (datas[i][0] > 0)
                        fprintf (fp[0], "%10.6lf", datas[i][0]);
                    else
                        fprintf (fp[0], "%10d", ( int )datas[i][0]);
                fprintf (fp[0], "\n");

                for (i = 0; i < slot->eta_qty; ++i) {
                    fprintf (fp[1], "%09d %6.2lf %6.2lf %6.2lf", slot->uid, slot->eta_all[i][0], slot->eta_all[i][1], slot->eta_all[i][2]);
                    if (datas[i][0] > 0)
                        fprintf (fp[1], "%10.6lf", datas[i][0]);
                    else
                        fprintf (fp[1], "%10d", ( int )datas[i][0]);
                    for (j = 1; j < 10; ++j)
                        fprintf (fp[1], " %10.6lf", datas[i][j]);
                    fprintf (fp[1], "\n");
                }
                fclose (fp[0]);
                fclose (fp[1]);
            }
        }
    }

    return 1;
}

// assign new task then drop into empty slot
static void enter () {
    int i;
    int is_pass = 0;
    int flag_hill[cst_node - 1];
    case_struct* case_pre;
    _complex y0[cst_node] = { 0.2, 0.2, 0.2 };
    slot_struct* slot = &QUEUE.slot[QUEUE.rear];

    slot->status = 2;
    while (CASE && slot->status == 2) {
        switch (CASE->case_type) {
        case 0:
            if (!topology (slot->ode.net_top, flag_hill, CASE->content[0])) {
                free (CASE->content);
                CASE = CASE->next;
                --global.total;
                break;
            }

            slot->status = 0;
            slot->ode.flag_record = 1;
            slot->case_id = CASE->case_id;
            slot->eta_qty = CASE->eta_qty;
            slot->case_type = CASE->case_type;
            slot->ode.para = CASE->content[1];
            slot->ode.cst_top[1] = ( _complex )flag_hill[0];
            slot->ode.cst_top[2] = ( _complex )flag_hill[1];
            slot->uid = CASE->content[0] * 10000 + (CASE->content[1]);
            memcpy (slot->ode.y, y0, cst_node * sizeof (_complex));
            for (i = 0; i < slot->eta_qty; ++i)
                slot->eta_all[i] = CASE->eta_all[i];
            slot->frame->next = NULL;
            CASE = CASE->next;
            break;
        case 1:
            if (CASE->progress[0] > CASE->content[1]) {
                free (CASE->content);
                CASE = CASE->next;
                break;
            }
            if (CASE->progress[1] >= cst_para || (is_pass = !topology (slot->ode.net_top, flag_hill, CASE->progress[0]))) {
                ++CASE->progress[0];
                CASE->progress[1] = 0;
                global.total -= is_pass ? cst_para : 0;
                break;
            }

            slot->status = 0;
            slot->ode.flag_record = 0;
            slot->case_id = CASE->case_id;
            slot->eta_qty = CASE->eta_qty;
            slot->case_type = CASE->case_type;
            slot->ode.para = CASE->progress[1];
            slot->ode.cst_top[1] = ( _complex )flag_hill[0];
            slot->ode.cst_top[2] = ( _complex )flag_hill[1];
            slot->uid = CASE->progress[0] * 10000 + (CASE->progress[1]);
            memcpy (slot->ode.y, y0, cst_node * sizeof (_complex));
            for (i = 0; i < slot->eta_qty; ++i)
                slot->eta_all[i] = CASE->eta_all[i];
            slot->frame->next = NULL;
            ++CASE->progress[1];
            break;
        case 2:
            if (CASE->progress[0] >= CASE->content_qty) {
                free (CASE->content);
                CASE = CASE->next;
            }
            if (!topology (slot->ode.net_top, flag_hill, CASE->content[CASE->progress[0]] / 10000)) {
                ++CASE->progress[0];
                --global.total;
                break;
            }

            slot->status = 0;
            slot->ode.flag_record = 0;
            slot->case_id = CASE->case_id;
            slot->eta_qty = CASE->eta_qty;
            slot->case_type = CASE->case_type;
            slot->uid = CASE->content[CASE->progress[0]];
            slot->ode.cst_top[1] = ( _complex )flag_hill[0];
            slot->ode.cst_top[2] = ( _complex )flag_hill[1];
            memcpy (slot->ode.y, y0, cst_node * sizeof (_complex));
            slot->ode.para = CASE->content[CASE->progress[0]] % 10000;
            for (i = 0; i < slot->eta_qty; ++i)
                slot->eta_all[i] = CASE->eta_all[i];
            slot->frame->next = NULL;
            ++CASE->progress[0];
            break;
        }
    }
}

// monitor task loop queue, output estimate remain time
static void monitor (int flag_record) {
    time_t clock = global.clock;
    int done = 0;

    for (;;) {
        switch (depart ()) {
        case 2:
            return;
        case 1:
            ++done;
            enter ();
            QUEUE.rear = (QUEUE.rear + 1) % QUEUE_MAX_LENGTH;
            if (!flag_record && global.total)
                printf ("\rprogress: %3d %%, estimate time remain: %.0f sec ", done * 100 / global.total, (global.total - done) * difftime (time (NULL), clock) / done);
        default:
            break;
        }
    }
}

// processor for calculating model
static void* processor () {
    size_t i, fetch, *front = &QUEUE.front;
    slot_struct* slot;

    for (;;) {
        fetch = -1;
        pthread_mutex_lock (&QUEUE.lock);
        fetch = *front;
        *front = (*front + 1) % QUEUE_MAX_LENGTH;
        pthread_mutex_unlock (&QUEUE.lock);
        if (!(fetch + 1))
            continue;
        slot = &QUEUE.slot[fetch];
        while (slot->status == 1)
            ;
        if (slot->status == 2)
            break;
        for (i = 0; i < slot->eta_qty; ++i) {
            memcpy (slot->ode.eta_list, slot->eta_all[i], cst_node * sizeof (double));
            if (inspect (slot->result[i], slot->frame, &slot->ode) == mark_timeout) {
                slot->result[0][0] = mark_timeout;
                break;
            }
        }
        slot->status = 1;
    }
    return ( void* )slot;
}

// initialize task loop queue
void init (case_struct* head, int eta_max, int flag_record) {
    int thread_qty = sysconf (_SC_NPROCESSORS_ONLN);
    pthread_t pid[thread_qty];
    case_struct* case_pre;
    size_t i, j;

    CASE = head;
    QUEUE.front = QUEUE.rear = 0;
    pthread_mutex_init (&QUEUE.lock, NULL);
    QUEUE.slot = ( slot_struct* )malloc (QUEUE_MAX_LENGTH * sizeof (slot_struct));
    for (i = 0; i < QUEUE_MAX_LENGTH; ++i) {
        slot_struct* slot = &QUEUE.slot[i];
        slot->status = 2;
        slot->frame = ( frame_struct* )malloc (sizeof (frame_struct));
        slot->result = ( double** )malloc (eta_max * sizeof (double*));
        slot->eta_all = ( double** )malloc (eta_max * sizeof (double*));
        for (j = 0; j < eta_max; ++j)
            slot->result[j] = ( double* )malloc (10 * sizeof (double));
        slot->uid = slot->eta_qty = slot->case_id = slot->case_type = slot->ode.para = slot->ode.flag_record = 0;
        CASE ? enter () : 1;
        QUEUE.rear = (QUEUE.rear + 1) % QUEUE_MAX_LENGTH;
    }

    QUEUE.rear = 0;
    for (i = 0; i < thread_qty; ++i)
        pthread_create (&pid[i], NULL, processor, NULL);
    monitor (flag_record);

    for (i = 0; i < thread_qty; ++i)
        pthread_join (pid[i], NULL);

    // free malloc space
    while (head) {
        for (i = 0; i < eta_max; ++i)
            free (head->eta_all[i]);
        free (head->eta_all);
        case_pre = head;
        head = head->next;
        free (case_pre);
    }
    for (i = 0; i < QUEUE_MAX_LENGTH; ++i) {
        for (j = 0; j < eta_max; ++j)
            free (QUEUE.slot[i].result[j]);
        free (QUEUE.slot[i].result);
        free (QUEUE.slot[i].eta_all);
    }
    free (QUEUE.slot);
    printf ("\r%s", flag_record ? "" : "\ndone.\n");
}