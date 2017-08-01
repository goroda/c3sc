// Copyright (c) 2015-2016, Massachusetts Institute of Technology
//
// This file is part of the C3 for Stochastic Optimal Control (C3SC) toolbox
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

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

#include "c3/stringmanip.h"


char * size_t_a_to_char(size_t * arr, size_t n, char * buffer)
{

    /* char * buffer = malloc(256 * sizeof(char)); */

    int cx = snprintf(buffer,256,"%zu ",arr[0]);
    for (size_t ii = 1; ii < n; ii++){
        int nx = snprintf(buffer+cx,256-cx,"%zu ",arr[ii]);
        cx += nx;
    }
    /* printf("%s",buffer); */
    return buffer;
}


/***********************************************************//**
    Simple hash function

    \param[in] size   - hashtable size
    \param[in] str_in - string to hash

    \return hashval

    \note
    http://www.sparknotes.com/cs/searching/hashtables/section3/page/2/
***************************************************************/
static size_t hashchar(size_t size, char * str_in)
{   
    size_t hashval = 0;
    char * str = str_in;
    for (; *str != '\0'; str++){
        hashval = (size_t) *str + (hashval << 5) - hashval;
        /* hashval =  (size_t) *str + (hashval << 6) + (hashval << 16) - hashval; */
        /* printf("\t hashval = %zu ", hashval ); */
    }
    /* printf("\n"); */
    
    return hashval % size;
}

struct HList{

//    double x;
    char key[256];
    double * data;
    size_t N;
    
    struct HList * next;
};

/* void hlist_push(struct HList ** head, char * string, void * data, size_t nbytes)//double xIn, size_t ind) */
void hlist_push(struct HList ** head, char * string, double * data, size_t N)//double xIn, size_t ind)    
{
    struct HList * newNode = malloc(sizeof(struct HList));
    /* newNode->key = NULL; */
    newNode->data = NULL;

    /* printf("push-key = %s\n",string); */
    newNode->next = *head;
    /* size_t N1 = strlen(string); */
    // create key
    /* newNode->key = malloc(256 * sizeof(char)); */
    strcpy(newNode->key,string);
    /* newNode->key = string; */

    // create data
    newNode->N = N;
    newNode->data = calloc(N,sizeof(double));
    assert (newNode->data != NULL);
    /* memmove((char *)(newNode->data), (char *) data, nbytes); */
    memmove(newNode->data, data, N * sizeof(double));
    
    *head = newNode;
}

void hlist_free(struct HList ** head){

    if (*head != NULL){
        // This is a bit inefficient and doesnt free the last node
        struct HList * current = *head;
        struct HList * next;
        while (current != NULL){
            next = current->next; current->next = NULL;
            /* printf("free key\n"); */
            /* free(current->key); current->key = NULL; */
            /* printf("freed data\n"); */
            free(current->data); current->data = NULL;
            /* printf("freed keys and data\n"); */
            /* free(current->x); */
            free(current); current = NULL;
            current = next;
        }
        *head = NULL;
    }
}

struct HTable
{
    size_t size;
    struct HList ** table;
};

struct HTable * htable_create(size_t size)
{

    struct HTable * new_table = NULL;
    if (size < 1) return NULL;
    
    if (NULL == (new_table = malloc(sizeof(struct HTable)))){
        fprintf(stderr, "failed to allocate memory for hashtable.\n");
        exit(1);
    }
    if (NULL == (new_table->table = malloc(size * sizeof(struct HList *)))){
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
    Free memory allocated to the hashtable

    \param[in,out] ht - hashtable
***************************************************************/
void htable_destroy(struct HTable * ht)
{
    if (ht != NULL){
        size_t ii;
        /* printf("destroy the little lists\n"); */
        for (ii = 0; ii < ht->size; ii++){
            hlist_free(&(ht->table[ii])); ht->table[ii] = NULL;
        }
        /* printf("destroyed all little lists\n"); */
        free(ht->table); ht->table = NULL;
        /* printf("destroyed table\n"); */
        free(ht); ht = NULL;
    }
}


/***********************************************************//**
    Lookup a key in the hashtable

    \param[in]     ht  - hashtable
    \param[in]     key - key to lookup
    \param[in,out] N   - number of values
    
    \return out - either NULL or a reference to the data
***************************************************************/
static double * lookup_key(struct HTable * ht, char * key, size_t * N)
{
    struct HList * pl = NULL;
    /* printf("lookup key = %s\n",key); */
    size_t val = hashchar(ht->size,key);
    /* printf("lookup key2 = %s\n",key); */

    /* printf("\t val = %zu\n",val); */
    /* printf("looking up key ind is = %zu\n",val); */
    *N = 0;
    double * out = NULL;

    /* for (pl = ht->table[val]; pl != NULL; pl = pl->next){ */
    /*     printf("key=%s, pl->key=%s\n",key,pl->key); */
    /*     out = pl->data; */
    /*     for (size_t jj = 0; jj < 10; jj++){ */
    /*         printf("%G ",((double*)out)[jj]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    for (pl = ht->table[val]; pl != NULL; pl = pl->next){
        if (strcmp(key,pl->key) == 0){
        /* if (memcmp(key,pl->key,4*sizeof(unsigned char)) == 0){             */
            /* printf("key=%s, pl->key=%s\n",key,pl->key); */
            out = pl->data;
            /* printf("pl-> N == %zu\n",pl->N); */
            *N = pl->N;
            return out;
        }
        /* printf("going to next list!\n"); */
    }
    /* printf("lookup key3 = %s\n",key); */
    return out;
}


/***********************************************************//**
    Add an index and value to the has table

    \param[in] ht   - hashtable
    \param[in] key  - must be a string! 
    \param[in] data - pointer to data
    \param[in] N    - number of values
    
    \return
        0 if good, 1 if some err,
***************************************************************/
int htable_add_element(struct HTable * ht, char * key, double * data, size_t N)
{
    size_t hashval = hashchar(ht->size,key);
    /* printf("add-key = %s\n",key);     */
    /* printf("adding hashval = %zu\n",hashval); */
    hlist_push(&(ht->table[hashval]),key,data,N);
    /* free(key); key = NULL; */

    return 0;
}

/***********************************************************//**
    Get the index associated with a particular value

    \param[in]     ht  - hashtable
    \param[in]     key - key to retrieve
    \param[in,out] N   - number of values
    
    \return pointer to data (or null)
***************************************************************/
double * htable_get_element(struct HTable * ht, char * key, size_t * N)
{

    // check if key already exists
    *N = 0;
    double * data = lookup_key(ht,key,N);
    return data;
}
