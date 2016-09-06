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

char * size_t_to_char(size_t n)
{

    char * out = malloc(256*sizeof(char));
    sprintf(out,"%zu",n);

    /* printf("n = %zu\n",n); */
    /* printf("sizeof(size_t) = %zu\n",sizeof(size_t)); */
    /* printf("strlen out = %zu\n",strlen(out)); */
    assert (strlen(out) > 0);

    return out;
}

char * size_t_a_to_char(size_t * arr, size_t n)
{

    char * out = malloc(256*sizeof(char));
    int offset = 0;
    for (size_t ii = 0; ii < n; ii++){
        int nchars = sprintf(out+offset,"%zu",arr[ii]);
        offset += nchars;
    }
    /* printf("offset = %d\n",offset); */

    /* printf("n = %zu\n",n); */
    /* printf("sizeof(size_t) = %zu\n",sizeof(size_t)); */
    /* printf("strlen out = %zu\n",strlen(out)); */
    /* assert (strlen(out) > 0); */
    /* exit(1); */
    return out;
}


/***********************************************************//**
    Simple hash function

    \param[in] size - hashtable size
    \param[in] str  - string to hash

    \return hashval

    \note
    http://www.sparknotes.com/cs/searching/hashtables/section3/page/2/
***************************************************************/
static size_t hashchar(size_t size, char * str)
{   
    size_t hashval = 0;
    
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
    char * key;
    void * data;
    size_t nbytes;
    
    struct HList * next;
};

void hlist_push(struct HList ** head, char * string, void * data, size_t nbytes)//double xIn, size_t ind)
{
    struct HList * newNode = malloc(sizeof(struct HList));
    newNode->key = NULL;
    newNode->data = NULL;
    
    newNode->next = *head;
    /* size_t N1 = strlen(string); */
    // create key
    newNode->key = malloc(256 * sizeof(char));
    sprintf(newNode->key,"%s",string);
    /* memmove(newNode->key,string,(N1+1)*sizeof(char)); */
    /* strncpy(newNode->key,string,N1); */
    /* newNode->key[N1] = '\0'; */

    // create data
    newNode->nbytes = nbytes;
    newNode->data = malloc(nbytes);
    memmove((char *)(newNode->data), (char *) data, nbytes);
    
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
            free(current->key); current->key = NULL;
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

    \param[in]     ht     - hashtable
    \param[in]     key    - key to lookup
    \param[in,out] nbytes - number of bytes of return
    
    \return out - either NULL or a reference to the data
***************************************************************/
static void * lookup_key(struct HTable * ht, char * key, size_t * nbytes)
{
    struct HList * pl = NULL;
    size_t val = hashchar(ht->size,key);
    
    /* printf("looking up key ind is = %zu\n",val); */
    *nbytes = 0;
    void * out = NULL;

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
            /* printf("key=%s, pl->key=%s\n",key,pl->key); */
            out = pl->data;
            *nbytes = pl->nbytes;
            return out;
        }
        /* printf("going to next list!\n"); */
    }
    return out;
}


/***********************************************************//**
    Add an index and value to the has table

    \param[in] ht     - hashtable
    \param[in] key    - must be a string! 
    \param[in] data   - pointer to data
    \param[in] nbytes - number of bytes of data (will copy!)
    
    \return
        0 if good, 1 if some err,
***************************************************************/
int htable_add_element(struct HTable * ht, char * key, void * data, size_t nbytes)
{

    size_t hashval = hashchar(ht->size,key);
    /* printf("adding hashval = %zu\n",hashval); */
    hlist_push(&(ht->table[hashval]),key,data,nbytes);
    free(key); key = NULL;

    return 0;
}

/***********************************************************//**
    Get the index associated with a particular value

    \param[in]     ht     - hashtable
    \param[in]     key    - double to add
    \param[in,out] nbytes - size in bytes of data
    
    \return pointer to data (or null)
***************************************************************/
void * htable_get_element(struct HTable * ht, char * key, size_t * nbytes)
{

    // check if key already exists
    *nbytes = 0;
    void * data = lookup_key(ht,key,nbytes);
    return data;
}
