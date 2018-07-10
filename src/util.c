// Copyright (c) 2015-2016, Massachusetts Institute of Technology
// Copyright (c) 2018, University of Michigan
//
// This file is part of the C3 for Stochastic Optimal Control (C3SC) toolbox
// Author: Alex A. Gorodetsky 
// Contact: goroda@umich.edu
// Website: https://www.alexgorodetsky.com

// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//Code



#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "cdyn/simulate.h"

/////////////////////////////////////////////
/// Packing a double
/// from 
/// http://stackoverflow.com/questions/3418319/serialize-double-and-float-with-c#3418560
////////////////////////////////////////////
#define FRAC_MAX 9223372036854775807LL /* 2**63 - 1 */

struct dbl_packed
{
    int exp;
    long long frac;
};

void pack(double x, struct dbl_packed *r)
{

    double xf = fabs(frexp(x, &(r->exp))) - 0.5;

    if (xf < 0.0)
    {
        r->frac = 0;
        return;
    }
    
    r->frac = 1 + (long long)(xf * 2.0 * (FRAC_MAX - 1));

    if (x < 0.0){
        r->frac = -r->frac;
    }
    printf("pack x = %G, exp=%d, xf=%G \n",x,r->exp,xf);
}

double unpack(const struct dbl_packed *p)
{
    double xf, x;

    if (p->frac == 0)
        return 0.0;

    /* xf = ((double)(llabs(p->frac) - 1) / (FRAC_MAX - 1)) / 2.0; */
    /* x = ldexp(xf + 0.5, p->exp); */

    xf = ((double)(llabs(p->frac) - 1) / (FRAC_MAX - 1)) / 2.0 + 0.5;
    /* x = xf * (1 << p->exp); */
    x = ldexp(xf + 0.5, p->exp);

    if (p->frac < 0)
        x = -x;

    return x;
}

/////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////

struct ApproxArgs
{
    double cross_tol;
    double round_tol;
    size_t kickrank;
    size_t startrank;
    size_t maxrank;
    int adapt;
    enum function_class fc;
};

struct ApproxArgs * approx_args_init(void)
{
    struct ApproxArgs * aargs = malloc(sizeof(struct ApproxArgs));
    if (aargs == NULL){
        fprintf(stderr,"Failure allocating approximation arguments\n");
        exit(1);
    }
    
    aargs->cross_tol = 1e-10;
    aargs->round_tol = 1e-10;
    aargs->kickrank = 10;
    aargs->startrank = 5;
    aargs->maxrank = 40;
    aargs->adapt = 1;
    aargs->fc = LINELM;
    return aargs;
}

void approx_args_free(struct ApproxArgs * aargs)
{
    if (aargs != NULL){
        free(aargs); aargs = NULL;
    }
}

void approx_args_set_function_class(struct ApproxArgs * aargs, enum function_class fc)
{
    assert (aargs != NULL);
    aargs->fc = fc;
}

enum function_class approx_args_get_function_class(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->fc;
}

void approx_args_set_cross_tol(struct ApproxArgs * aargs, double cross_tol)
{
    assert (aargs != NULL);
    aargs->cross_tol = cross_tol;
}

double approx_args_get_cross_tol(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->cross_tol;

}
void approx_args_set_round_tol(struct ApproxArgs * aargs, double round_tol)
{
    assert (aargs != NULL);
    aargs->round_tol = round_tol;
}

double approx_args_get_round_tol(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->round_tol;
}

void approx_args_set_kickrank(struct ApproxArgs * aargs, size_t kickrank)
{
    assert (aargs != NULL);
    aargs->kickrank = kickrank;
}

size_t approx_args_get_kickrank(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->kickrank;
}

void approx_args_set_maxrank(struct ApproxArgs * aargs, size_t maxrank)
{
    assert (aargs != NULL);
    aargs->maxrank = maxrank;
}

size_t approx_args_get_maxrank(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->maxrank;
}

void approx_args_set_startrank(struct ApproxArgs * aargs, size_t startrank)
{
    assert (aargs != NULL);
    aargs->startrank = startrank;
}

size_t approx_args_get_startrank(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->startrank;
}

void approx_args_set_adapt(struct ApproxArgs * aargs, int adapt)
{
    assert (aargs != NULL);
    aargs->adapt = adapt;
}

int approx_args_get_adapt(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->adapt;
}

int c3sc_check_bounds(size_t dx, double * lbx,
                      double * ubx,const  double * x)
{
    if ((lbx == NULL) || (ubx == NULL)){
        return 0;
    }
    for (size_t ii = 0; ii < dx; ii++){
        if (x[ii] < lbx[ii]){
//            printf("x-lbx = %G",x[ii]-lbx[ii]);
            return - (ii+1);
        }
        else if (x[ii] > ubx[ii]){
            return ii+1;
        }
    }
    return 0;
}

static int compareDouble (const void * a, const void * b)
{
    if (*(const double*)a < *(const double*)b){
        return -1;
    }
    else if (*(const double*)a > *(const double*)b){
        return 1;
    }
    return 0;
}

double * c3sc_combine_and_sort(size_t Nx, double * x, size_t Ny, double * y,
                               size_t * Ntot)
{
    double * temp = malloc( (Nx+Ny)*sizeof(double));
    assert (temp != NULL);
    memmove(temp,x,Nx*sizeof(double));
    memmove(temp+Nx,y,Ny*sizeof(double));
    qsort(temp,Nx+Ny,sizeof(double),compareDouble);

    double * out = calloc(Nx+Ny,sizeof(double));
    *Ntot = 1;
    out[0] = temp[0];
    for (size_t ii = 1; ii < Nx+Ny; ii++){
        if (fabs(temp[ii]-out[*Ntot-1]) > 1e-15){
            out[*Ntot] = temp[ii];
            *Ntot = *Ntot + 1;
        }
    }
    free(temp); temp = NULL;
    return out;
}


struct c3sc_SortCouple
{
    size_t ind;
    double val;
};
    
int c3sc_compare (const void * a, const void * b)
{
    const struct c3sc_SortCouple *as = (const struct c3sc_SortCouple *)a;
    const struct c3sc_SortCouple *bs = (const struct c3sc_SortCouple *)b;
    if ( as->val < bs->val ){
      return -1;  
    } 
    else if ( as->val > bs->val ){
      return 1;  
    }
    else{
        return 0;
    }
}

// probs is overwritten
size_t c3sc_sample_discrete_rv(size_t n, double * probs,
                               double sample)
{
    struct c3sc_SortCouple sc[1000];        
    if (n > 1000){
        fprintf(stderr,"Not enough memory allocate in ");
        fprintf(stderr,"discrete_sample\n");
        exit(1);
    }

    for (size_t ii = 0; ii < n; ii++){
        sc[ii].ind = ii;
        sc[ii].val = probs[ii];
//        perm[ii] = ii+1;
    }

    // sort in ascending order
    qsort (sc,n,sizeof(struct c3sc_SortCouple),c3sc_compare);

    probs[0] = sc[0].val;
    for (size_t ii = 1; ii < n; ii++){
        probs[ii] = sc[ii].val + probs[ii-1];
    }
    for (size_t jj = 0; jj < n; jj++){
        if (sample <= probs[jj]){
            return sc[jj].ind;
        }
    }
    
    fprintf(stderr,"Warning: problem with gen. sample\n");
    fprintf(stderr,"Uniform sample is %G\n",sample);
    return 0;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////// Hashing doubles /////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/***********************************************************//**
    Simple hash function

    \param[in] size - hashtable size
    \param[in] str  - string to hash

    \return hashval

    \note
    http://www.sparknotes.com/cs/searching/hashtables/section3/page/2/
***************************************************************/
static size_t hashd(size_t size, char * str)
{   
    size_t hashval = 0;
    
    for (; *str != '\0'; str++){
        hashval = (size_t) *str + (hashval << 5) - hashval;
        /* printf("\t hashval = %zu ", hashval ); */
    }
    /* printf("\n"); */
    
    return hashval % size;
}

struct HashGridList{

//    double x;
    char * key;
    size_t ind;
    struct HashGridList * next;
};

void hash_grid_list_push(struct HashGridList ** head, size_t ind,char * string)//double xIn, size_t ind)
{
    struct HashGridList * newNode = malloc(sizeof(struct HashGridList));
    //  newNode->x = xIn;
    newNode->ind = ind;
    newNode->next = *head;
    size_t N1 = strlen(string);
    newNode->key = malloc((N1+1)*sizeof(char));
    strncpy(newNode->key,string,N1);
    newNode->key[N1] = '\0';
    *head = newNode;
}

void hash_grid_list_free(struct HashGridList ** head){

    // This is a bit inefficient and doesnt free the last node
    struct HashGridList * current = *head;
    struct HashGridList * next;
    while (current != NULL){
        next = current->next; current->next = NULL;
        free(current->key); current->key = NULL;
        /* free(current->x); */
        free(current); current = NULL;
        current = next;
    }
    *head = NULL;
}

struct HashGrid
{
    size_t size;
    struct HashGridList ** table;
};

/***********************************************************//**
    Allocate memory for a new hashtable of cpairs

    \param[in] size - size table
    
    \return new_table 
***************************************************************/
struct HashGrid * hash_grid_create(size_t size)
{

    struct HashGrid * new_table = NULL;
    if (size < 1) return NULL;
    
    if (NULL == (new_table = malloc(sizeof(struct HashGrid)))){
        fprintf(stderr, "failed to allocate memory for hashtable.\n");
        exit(1);
    }
    if (NULL == (new_table->table = malloc(size * sizeof(struct HashGridList *)))){
        fprintf(stderr, "failed to allocate memory for hashtable.\n");
        exit(1);
    }

    new_table->size = size;
    size_t ii;
    for (ii = 0; ii < size; ii++){
        new_table->table[ii] = NULL;
    }

    return new_table;
}

/***********************************************************//**
    Allocate memory for and has a grid

    \param[in] size - size table
    \param[in] grid - grid to has
    
    \return new_table 
***************************************************************/
struct HashGrid * hash_grid_create_grid(size_t size, const struct c3Vector * grid)
{
    struct HashGrid * hg = hash_grid_create(size);
    for (size_t ii = 0; ii < grid->size; ii++){
        hash_grid_add_element(hg,ii,grid->elem[ii]);
    }
    return hg;
}

/***********************************************************//**
    Allocate memory for and hash a grid

    \param[in] size - size table
    \param[in] d    - number of dimensions
    \param[in] grid - grids to hash
    
    \return new_table 
***************************************************************/
struct HashGrid ** hash_grid_create_ndgrid(size_t size, size_t d, struct c3Vector ** grid)
{
    struct HashGrid ** hg = malloc(d * sizeof(struct HashGrid *));
    for (size_t jj = 0; jj < d; jj++){
        hg[jj] = hash_grid_create(size);
        for (size_t ii = 0; ii < grid[jj]->size; ii++){
            hash_grid_add_element(hg[jj],ii,grid[jj]->elem[ii]);
        }
    }
    return hg;
}

/***********************************************************//**
    Free an array of hash grids
***************************************************************/
void hash_grid_free_ndgrid(size_t d, struct HashGrid ** hg)
{
    if (hg != NULL){
        for (size_t ii = 0; ii < d; ii++){
            hash_grid_free(hg[ii]); hg[ii] = NULL;
        }
        free(hg); hg = NULL;
    }
}

/***********************************************************//**
    Lookup a key in the hashtable

    \param[in]     ht     - hashtable
    \param[in]     key    - key to lookup
    \param[in,out] exists - returns 1 if exists 0 if doesnt
    
    \return out - either NULL or the second element in the pair stored under the key 
***************************************************************/
static size_t lookup_keyd(struct HashGrid * ht, char * key, int * exists)
{
    struct HashGridList * pl = NULL;
    size_t val = hashd(ht->size,key);
    
    /* printf("val ind = %zu\n",val); */
    *exists = 0;
    size_t out = 0;
    for (pl = ht->table[val]; pl != NULL; pl = pl->next){
        if (strcmp(key,pl->key) == 0){
            out = pl->ind;
//            out[N]='\0';
            *exists = 1;
            return out;
        }
    }
    return out;
}

/***********************************************************//**
    Add an index and value to the has table

    \param[in] ht  - hashtable
    \param[in] ind - index
    \param[in] val - double to add
    
    \return
        0 if good, 1 if some err, 2 if exists
***************************************************************/
int hash_grid_add_element(struct HashGrid * ht,size_t ind,double val)
{

    /* struct dbl_packed dbl; */
    /* pack(val,&dbl); */
    /* printf("exp = %d, frac=%lld \n",dbl.exp,dbl.frac); */


    /* char last = '\0'; */
    /* char * key = malloc(nvals); */
    /* assert (key != NULL); */
    /* memmove(key,&(dbl.exp),sizeof(int)); */
    /* memmove(key+sizeof(int),(char *)&(dbl.frac),sizeof(long long)); */
    /* memmove(key+sizeof(int)+sizeof(long long),&last,sizeof(char)); */
    char * key = serialize_double_to_text(val);
    /* printf("%s\n",key); */
    /* key[nvals-1] = '\0'; */
    /* char  * key = "asaksjdhaskdh\0"; */
    /* for (size_t ii = 0; ii < nvals; ii++){ */
    /*     printf("%c",key[ii]); */
    /* } */
    /* printf("\n"); */
    
    // check if key already exists
    int exists = 0;
    lookup_keyd(ht,key,&exists);
    if (exists == 1) {
        free(key); key = NULL;
        return 2;
    }

    size_t hashval = hashd(ht->size,key);
    /* printf("adding hashval = %zu\n",hashval); */
    hash_grid_list_push(&(ht->table[hashval]),ind,key);
    free(key); key = NULL;

    return 0;
}

void hash_grid_print(struct HashGrid * ht, FILE *fp)
{
    for (size_t ii = 0; ii < ht->size; ii++){
        struct HashGridList * hgl = ht->table[ii];
        if (hgl != NULL){
            while (hgl != NULL){
                double val = deserialize_double_from_text(hgl->key);
                fprintf(fp,"ind=%zu,val = %3.15G ",ii,val);
                hgl = hgl->next;
            }
            fprintf(fp,"\n");
        }
    }
}

/***********************************************************//**
    Get the index associated with a particular value

    \param[in] ht  - hashtable
    \param[in] val - double to add
    
    \return index
***************************************************************/
size_t hash_grid_get_ind(struct HashGrid * ht,double val,int *exists)
{

    char * key = serialize_double_to_text(val);
    /* struct dbl_packed dbl; */
    /* pack(val,&dbl); */
    /* size_t nvals = sizeof(int) + sizeof(long long)+1; */
    /* char * key = malloc(nvals); */
    /* assert (key != NULL); */
    /* memmove(key,&(dbl.exp),sizeof(int)); */
    /* memmove(key+sizeof(int),&(dbl.frac),sizeof(long long)); */
    /* key[nvals-1] = '\0'; */
    
    // check if key already exists
    int one = 1;
    *exists = one;
    size_t ind = lookup_keyd(ht,key,exists);
    if (*exists == 0) {
        /* fprintf(stderr,"Value doesn't exists %3.15G\n",val); */
        /* assert (1 == 0); */
        free(key);
        return 0;
        /* free(key); key = NULL; */
        /* return 2; */
    }
    free(key); key = NULL;
    return ind;
}

/***********************************************************//**
    get index associated with an array of grids

    \return 0 if got index, 1 if didnt
***************************************************************/
int hash_grid_ndgrid_get_ind(struct HashGrid ** ht, size_t dim, const double * x, size_t * out)
{
    int success = 0;
    int exists;
    for (size_t ii = 0; ii < dim; ii++){
        out[ii] = hash_grid_get_ind(ht[ii],x[ii],&exists);
        if (exists == 0){
            //success = 1;
            /* fprintf(stderr,"value does not exist in array\n"); */
            /* dprint(dim,x); */
            /* exit(1); */
            success = 1;
            break;
        }
    }
    return success;
}

/***********************************************************//**
    Free memory allocated to the hashtable

    \param[in,out] ht - hashtable
***************************************************************/
void hash_grid_free(struct HashGrid * ht)
{
    if (ht != NULL){
        size_t ii;
        for (ii = 0; ii < ht->size; ii++){
            hash_grid_list_free(&(ht->table[ii]));
        }
        free(ht->table); ht->table = NULL;
        free(ht); ht = NULL;
    }
}

/* unsigned char * serialize_int(unsigned char *buffer, int value) */
/* { */
/*     buffer[0] = value >> 24; */
/*     buffer[1] = value >> 16; */
/*     buffer[2] = value >> 8; */
/*     buffer[3] = value; */
/*     return buffer + 4; */
/* } */



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//                       Workspace manager                                 //

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

struct Workspace
{
    size_t dx;
    size_t du;
    size_t dw;
    size_t N; // copies of memory
    size_t index[11]; // indexs into memory mem
    double ** mem; // memory for node based variables

    // separate memory for fiber-based
    double * costs_pol;
    int    * absorbed_pol;
    size_t * ind_to_serialize;
    size_t * absorbed_no;
    size_t * absorbed_yes;


    // Hash tables for various things
    struct HTable * vi_htable;
    size_t vi_iter;

    struct HTable * pi_htable;
    struct HTable * pi_prob_htable;
    size_t pi_iter;
    size_t pi_subiter;

    char ** saved_keys;
    char ** saved_keys2;
};

struct Workspace * workspace_alloc(size_t dx,size_t du,size_t dw, size_t N)
{
    struct Workspace * work = malloc(sizeof(struct Workspace));
    if (work == NULL){
        fprintf(stderr,"Workspace cannot be allocated\n");
        exit(1);
    }
    
    work->dx = dx;
    work->du = du;
    work->dw = dw;
    work->N = N;
    work->mem = malloc(N * sizeof(double *));
    if (work->mem == NULL){
        fprintf(stderr,"Workspace cannot be allocated\n");
        exit(1);
    }
    
    // integers denote space in memory required by specified variable
    work->index[0] = dx;                             // drift
    work->index[1] = dx * du + work->index[0];       // grad_drift
    work->index[2] = dx * dw + work->index[1];       // diff
    work->index[3] = dx * dw * du + work->index[2];  // grad_diff
    work->index[4] = 1 + work->index[3];             // dt
    work->index[5] = du + work->index[4];            // grad_dt
    work->index[6] = 2 * dx + 1 + work->index[5];    // prob
    work->index[7] = du * (2*dx+1) + work->index[6]; // grad_prob
    work->index[8] = du + work->index[7];            // grad_stage
    work->index[9] = du + work->index[8];            // extra control-sized mem
    work->index[10] = du + work->index[9];           // u

    for (size_t ii = 0; ii < N; ii++){
        work->mem[ii] = calloc_double(work->index[10]);
    }

    work->costs_pol = calloc_double(N*(2*dx+1));
    work->absorbed_pol = calloc_int(N);
    work->ind_to_serialize = calloc_size_t(dx+3);
    work->absorbed_no = calloc_size_t(N);
    work->absorbed_yes = calloc_size_t(N);

    size_t ncreate = 1000000; // original
    /* size_t ncreate = 10000; // original     */
    work->vi_htable = htable_create(ncreate);
    work->vi_iter = 0;

    work->pi_htable = htable_create(ncreate);
    work->pi_prob_htable = htable_create(ncreate);
    work->pi_iter = 0;
    work->pi_subiter = 0;

    work->saved_keys = malloc(N * sizeof(char*));
    work->saved_keys2 = malloc(N * sizeof(char*));
    for (size_t ii = 0; ii < N; ii++){
        work->saved_keys[ii] = malloc(256*sizeof(char));
        work->saved_keys2[ii] = malloc(256*sizeof(char));
    }
    /* assert (work->saved_keys != NULL); */
    
    return work;
}

void workspace_reset_pi_prob_htable(struct Workspace * w)
{
    size_t ncreate = 1000000; // original
    htable_destroy(w->pi_prob_htable);
    w->pi_prob_htable = htable_create(ncreate);
}

void workspace_reset_pi_htable(struct Workspace * w)
{
    size_t ncreate = 1000000; // original
    htable_destroy(w->pi_htable);
    w->pi_htable = htable_create(ncreate);
}

void workspace_reset_vi_htable(struct Workspace * w)
{
    size_t ncreate = 1000000; // original
    htable_destroy(w->vi_htable);
    w->vi_htable = htable_create(ncreate);
}

void workspace_free(struct Workspace * w)
{
    if (w != NULL){
        for (size_t ii = 0; ii <  w->N; ii++){
            free(w->mem[ii]); w->mem[ii] = NULL;
        }
        free(w->mem); w->mem = NULL;
        free(w->costs_pol); w->costs_pol = NULL;
        free(w->absorbed_pol); w->absorbed_pol = NULL;
        free(w->ind_to_serialize); w->ind_to_serialize = NULL;
        free(w->absorbed_no); w->absorbed_no = NULL;
        free(w->absorbed_yes); w->absorbed_yes = NULL;

        htable_destroy(w->vi_htable); w->vi_htable = NULL;
        htable_destroy(w->pi_htable); w->pi_htable = NULL;
        htable_destroy(w->pi_prob_htable); w->pi_prob_htable = NULL;

        for (size_t ii = 0; ii < w->N; ii++){
            free(w->saved_keys[ii]); w->saved_keys[ii] = NULL;
            free(w->saved_keys2[ii]); w->saved_keys2[ii] = NULL;
        }
        free(w->saved_keys); w->saved_keys = NULL;
        free(w->saved_keys2); w->saved_keys2 = NULL;
        free(w); w = NULL;
    }
}

void workspace_increment_vi_iter(struct Workspace * w)
{
    w->vi_iter++;
}

size_t workspace_get_vi_iter(const struct Workspace * w)
{
    return w->vi_iter;
}

void workspace_increment_pi_iter(struct Workspace * w)
{
    w->pi_iter++;
}

size_t workspace_get_pi_iter(const struct Workspace * w)
{
    return w->pi_iter;
}

void workspace_increment_pi_subiter(struct Workspace * w)
{
    w->pi_subiter++;
}

size_t workspace_get_pi_subiter(const struct Workspace * w)
{
    return w->pi_subiter;
}


struct HTable * workspace_get_vi_htable(const struct Workspace * w)
{
    return w->vi_htable;
}

struct HTable * workspace_get_pi_htable(const struct Workspace * w)
{
    return w->pi_htable;
}

struct HTable * workspace_get_pi_prob_htable(const struct Workspace * w)
{
    return w->pi_prob_htable;
}


double * workspace_get_drift(struct Workspace * w, size_t node)
{
    return w->mem[node];
}

double * workspace_get_grad_drift(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[0];
}

double * workspace_get_diff(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[1];
}

double * workspace_get_grad_diff(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[2];
}

double * workspace_get_dt(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[3];
}

double * workspace_get_grad_dt(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[4];
}

double * workspace_get_prob(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[5];
}

double * workspace_get_grad_prob(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[6];
}

double * workspace_get_grad_stage(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[7];
}

double * workspace_get_control_size_extra(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[8];
}

double * workspace_get_u(struct Workspace * w, size_t node)
{
    return w->mem[node] + w->index[9];
}

double * workspace_get_costs(struct Workspace * w, size_t node)
{
    return w->costs_pol + node * (2 * w->dx + 1);
}

int * workspace_get_absorbed(struct Workspace * w, size_t node)
{
    return w->absorbed_pol + node;
}

size_t * workspace_get_ind_to_serialize(struct Workspace * w)
{
    return w->ind_to_serialize;
}

size_t * workspace_get_absorbed_no(struct Workspace * w)
{
    return w->absorbed_no;
}

size_t * workspace_get_absorbed_yes(struct Workspace * w)
{
    return w->absorbed_yes;
}

char ** workspace_get_saved_keys(struct Workspace * w)
{
    return w->saved_keys;
}

char ** workspace_get_saved_keys2(struct Workspace * w)
{
    return w->saved_keys2;
}
    

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//                       Other useful stuff                                //

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


/**********************************************************//**
    Compute the stride length between indices to choose
    *M* indices almost uniformly from *N* options

    \note
    N must be greater than or equal to M (N>=M)
**************************************************************/
size_t uniform_stride(size_t N, size_t M)
{
    assert (N >= M);
    size_t stride = 1;
    size_t M1 = M-1;
    size_t N1 = N-1;
    while (stride * M1 < N1){
        stride++;
    }
    stride--;
    return stride;
}
