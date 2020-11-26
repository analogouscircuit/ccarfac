#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ------------------------- STRUCTURES ------------------------- */

/** Node containing a piece of data (and associated metadeta) that one wishes
 * to transfer to Python.  This is the main component of the linked list
 * structure used by the library */
typedef struct pydata_node_s {
    char * name;
    char type_py;
    size_t size;
    unsigned int len;
    int * dims;
    int num_dims;
    void * data;
    struct pydata_node_s * next;
} pydata_node;

/** Structure for the actual linked list (stringing together the nodes) */
typedef struct pydata_list_s {
    pydata_node * head;
} pydata_list;



/* ------------------------- FUNCTION DECLARATIONS ------------------------- */

pydata_list * pydata_list_init();

void pydata_list_push(void * data, char * type, 
                      int num_dims, int * dims, 
                      char * name, pydata_list * l);
    
pydata_node * pydata_list_pop(pydata_list * l);

void pydata_node_free(pydata_node * node);

void pydata_list_free(pydata_list * l);

void pydata_write(pydata_list * l, char * file_name);
