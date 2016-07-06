#ifndef HASHTABLE_H
#define HASHTABLE_H

char * size_t_to_char(size_t);
char * size_t_a_to_char(size_t *, size_t);

struct HTable;
struct HTable * htable_create(size_t);
void htable_destroy(struct HTable *);
int htable_add_element(struct HTable *, char *,void*,size_t);
void * htable_get_element(struct HTable *, char *, size_t *);

#endif
