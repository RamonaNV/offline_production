#ifndef _SPLINE_TABLE_H
#define _SPLINE_TABLE_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct splinetable {
	int ndim;
	int *order;

	double **knots;
	long *nknots;

	double **extents;

	double *periods;

	float *coefficients;
	long *naxes;
	unsigned long *strides;

	int naux;
	char ***aux;
};

struct splinetable_buffer {
	void *data;
	size_t size;
	void *(*mem_alloc)(size_t newsize);
	void *(*mem_realloc)(void *p, size_t newsize);
};

typedef enum {
	SPLINETABLE_INT,
	SPLINETABLE_DOUBLE
} splinetable_dtype;

int readsplinefitstable(const char *path, struct splinetable *table);
int readsplinefitstable_mem(struct splinetable_buffer *buffer,
    struct splinetable *table);
int writesplinefitstable(const char *path, const struct splinetable *table);
int writesplinefitstable_mem(struct splinetable_buffer *buffer,
    const struct splinetable *table);
void splinetable_free(struct splinetable *table);
void splinetable_permute(struct splinetable *table, int *permutation);
char * splinetable_get_key(const struct splinetable *table, const char *key);
int splinetable_read_key(const struct splinetable *table, splinetable_dtype type,
    const char *key, void *result);

#ifdef __cplusplus
}
#endif

#endif /* _SPLINE_TABLE_H */

