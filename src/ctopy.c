/** \file
 * A simple library for transferring numerical data from C into NumPy arrays in
 * Python. The main usage case: plotting with Matplotlib. */

#include "ctopy.h"

/** \mainpage
 * This library provides a simple way to wrap up numerical data in C and save it
 * to a pair of files (one text, one binary) that can then be read by a Python
 * function.  This Python function then generates a dictionary of NumPy arrays.
 * This is especially useful when hacking together numerical algorithms in C,
 * the output of which one wishes to plot in Matplotlib.
 *
 * The main user-facing functions are: \n\n
 * \ref pydata_list_init, for initializing a list of data \n 
 * \ref pydata_list_push, for adding data to the list \n 
 * \ref pydata_write, for writing the data to files \n 
 * \ref pydata_list_free, for freeing memory when done \n\n
 *
 * The order listed above is the standard order for use. */


/** Initialize the data structure (list) containing data to pass to Python
 *
 * \return A pydata_list, ready for accepting data */
pydata_list * pydata_list_init()
{
	pydata_list * l = malloc(sizeof(pydata_list));
	l->head = NULL;
	return l;
}


/** Add data to the list of data that will be trasnferred to Python 
 *
 * \param data A pointer to the actual data that will be transferred.  Should
 * be of type * float, * double, or * int. (Note user retains responsibility for
 * freeing this data when done.)
 * 
 * \param type A string inidicating the type. Should be "float", "double", or 
 * "int"
 *
 * \param num_dims An integer indicating the number of dimensions of the
 * vector/matrix/tensor
 *
 * \param dims An array of integers (which must be the same length as
 * num_dims!) indicating the length of each dimension.
 *
 * \param name A string that determines the variable name in the Python
 * environment
 *
 * \param pydata_list A pointer to a pydata_list 
 *
 * \return void */
void pydata_list_push(void * data, char * type, 
					  int num_dims, int * dims,
					  char * name, pydata_list * l)
{
	size_t size;
	char type_py;
	pydata_node * node = malloc(sizeof(pydata_node));

	// Determine the size of the data type and record label for Python
	if(strcmp(type, "float") == 0) {
		size = sizeof(float);
		type_py = 'f';
	}
	else if(strcmp(type, "double") == 0) {
		size = sizeof(double);
		type_py = 'd';
	}
	else if(strcmp(type, "int") == 0) {
		size = sizeof(int);
		type_py = 'i';
	}
	else {
		printf("Invalid type given to ctopy. Type must be float, double or int. Exiting");
		exit(0);
	}

	// Resolve the dimensions and length variables
	unsigned int len = 1;
	int * node_dims = malloc(sizeof(int)*num_dims);
	for(int k=0; k < num_dims; k++) {
		node_dims[k] = dims[k];
		len *= dims[k];
	}

	// Set the values of the new node
	node->size = size;
	node->type_py = type_py;
	int name_len = strlen(name);
	node->name = malloc(sizeof(char)*(1+name_len));
	node->name[name_len] = '\0';
	strcpy(node->name, name); 
	node->len = len;
	node->data = data;
	node->num_dims = num_dims;
	node->dims = node_dims;

	// Linked list book keeping
	node->next = l->head;
	l->head = node;
}

/** Function for popping off nodes from the list of data. Not typically called by user. 
 *
 * \param l The list from which one wishes to pop the top element
 * \return The node that is popped off the list */
pydata_node * pydata_list_pop(pydata_list * l)
{
	pydata_node * node_out;
	if(l->head != NULL) {
		node_out = l->head;
		l->head = l->head->next;
		return node_out;
	}
	return NULL;
}

/** Function for freeing a node (and all data to which it points!)
 *
 * \param node The node that one wishes to free
 * \return void */
void pydata_node_free(pydata_node * node)
{
	free(node->name);
	// free(node->data); // make user free their own data
	free(node->dims);
	free(node);
}

/** Function for freeing all the data one wishes to pass to python. Frees all
 * nodes and all data that the nodes point to. Should be called before program
 * termination.
 *
 * \param l List that of data that one wishes to free 
 * \return void */
void pydata_list_free(pydata_list * l)
{
	pydata_node * node;
	while(l->head != NULL) {
		node = pydata_list_pop(l);
		pydata_node_free(node);
	}
	free(l);
}

/** Function for actually writing the data to two files (one text metadata, one
 * binary), to be read by Python. 
 *
 * \param l List of data that one wishes to transfer
 *
 * \param file_name Prefix for both file names 
 *
 * \return void */

void pydata_write(pydata_list * l, char * file_name)
{
	FILE * file_txt;
	FILE * file_bin;
	char * file_name_txt;
	char * file_name_bin;
	int name_len = strlen(file_name);

	// Set up the file names
	file_name_txt = malloc(sizeof(char)*(5+name_len));
	file_name_bin = malloc(sizeof(char)*(5+name_len));
	file_name_txt[name_len+4] = '\0';
	file_name_bin[name_len+4] = '\0';
	strncpy(file_name_txt, file_name, name_len);
	strncpy(file_name_bin, file_name, name_len);
	strcpy(file_name_txt+name_len, ".txt");
	strcpy(file_name_bin+name_len, ".bin");

	// Write the data 
	file_txt = fopen(file_name_txt, "w");
	file_bin = fopen(file_name_bin, "w");

	pydata_node * node = pydata_list_pop(l);
	while(node != NULL) {
		fprintf(file_txt, "%s\t%c\t%d\t", node->name, node->type_py, node->num_dims);
		for(int k = 0; k < node->num_dims; k++) {
			fprintf(file_txt,"%d\t", node->dims[k]);
		}
		fprintf(file_txt, "\n");
		fwrite(node->data, node->size, node->len, file_bin);
		pydata_node_free(node);
		node = pydata_list_pop(l);
	}
	
	// Clean up
	fclose(file_txt);
	fclose(file_bin);
	free(file_name_txt);
	free(file_name_bin);
}

