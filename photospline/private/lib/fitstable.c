/*
 * fitstable.c: Reads spline coefficient tables from on-disk FITS tables.
 *
 * These tables are laid out in the following way:
 *  Header Keys:
 *   TYPE: "Spline Coefficient Table"
 *   ORDERn: Order of B-splines on the n-th axis
 *   PERIODn: Periodicity of the n-th axis, or 0 if not a periodic basis
 *   BIAS: Logarithm offset (optional)
 *   GEOMETRY: Photonics geometry (optional)
 *  Images:
 *   Primary: N-D array of spline coefficients
 *   KNOTSn: Vector of knot locations on axis n
 *   EXTENTS: 2-D array of table boundaries (optional)
 *   
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "photospline/splinetable.h"

static int parsefitstable(fitsfile *fits, struct splinetable *table);
static int fillfitstable(fitsfile *fits, const struct splinetable *table);

int
readsplinefitstable(const char *path, struct splinetable *table)
{
	fitsfile *fits;
	int error = 0;

	memset(table, 0, sizeof(struct splinetable));

	fits_open_file(&fits, path, READONLY, &error);
	if (error != 0)
		return (error);

	error = parsefitstable(fits, table);
	fits_close_file(fits, &error);
	fits_report_error(stderr, error);

	return (error);
}

int
readsplinefitstable_mem(struct splinetable_buffer *buf,
    struct splinetable *table)
{
	fitsfile *fits;
	int error = 0;

	memset(table, 0, sizeof(struct splinetable));

	fits_open_memfile(&fits, "", READONLY, &buf->data, &buf->size,
	    0, NULL, &error);
	if (error != 0)
		return (error);

	error = parsefitstable(fits, table);
	fits_close_file(fits, &error);
	fits_report_error(stderr, error);

	return (error);
}

int
writesplinefitstable(const char *path, const struct splinetable *table)
{
	fitsfile *fits;
	int error = 0;
	
	fits_create_diskfile(&fits, path, &error);
	if (error != 0)
		return (error);
	
	error = fillfitstable(fits, table);
	fits_close_file(fits, &error);
	fits_report_error(stderr, error);
	
	return (error);
}

int
writesplinefitstable_mem(struct splinetable_buffer *buf,
    const struct splinetable *table)
{
	fitsfile *fits;
	int error = 0;
	
	buf->size = 2880;
	buf->data = buf->mem_alloc(buf->size);
	if (!buf->data)
		return (ENOMEM);
	fits_create_memfile(&fits, &buf->data, &buf->size, 0u,
	    buf->mem_realloc, &error);
	if (error != 0)
		return (error);
	
	error = fillfitstable(fits, table);
	fits_close_file(fits, &error);
	fits_report_error(stderr, error);
	
	return (error);
}

static void
splinetable_free_aux(struct splinetable *table)
{
	int i;

	for (i = 0; i < table->naux; i++) {
		free(table->aux[i][0]);
		free(table->aux[i][1]);
		free(table->aux[i]);
	}

	free(table->aux);
	table->naux = 0;
	table->aux = NULL;
}

void splinetable_free(struct splinetable *table)
{
	int i;

	splinetable_free_aux(table);
	if (table->knots) {
		for (i = 0; i < table->ndim; i++)
			if (table->knots[i])
				free(table->knots[i]-table->order[i]);
		free(table->knots);
	}
	free(table->nknots);
	free(table->naxes);
	free(table->coefficients);
	free(table->periods);
	free(table->strides);
	
	if (table->extents) {
		free(table->extents[0]);
		free(table->extents);
	}
}

/* Re-order the dimensions of a spline table */
void
splinetable_permute(struct splinetable *table, int *permutation)
{
	int i, j;
	int iperm[table->ndim];
	unsigned long pos, npos;
	unsigned long total_size;
	int *order;
	unsigned long *strides;
	long *naxes, *nknots;
	double **extents, **knots;
	float *coefficients;
	
	assert(table->ndim >= 1);
	order   = malloc(sizeof(table->order[0])*table->ndim);
	naxes   = malloc(sizeof(table->naxes[0])*table->ndim);
	strides = malloc(sizeof(table->strides[0])*table->ndim);
	nknots  = malloc(sizeof(table->nknots[0])*table->ndim);
	knots   = malloc(sizeof(table->knots[0])*table->ndim);
	extents = malloc(sizeof(table->extents[0])*table->ndim);
	extents[0] = malloc(sizeof(table->extents[0][0])*2*table->ndim);
	for (i = 1; i < table->ndim; i++) {
		extents[i] = &extents[0][2*i];
	}
	
	/* Permute various per-axis properties */
	for (i = 0; i < table->ndim; i++) {
		j = permutation[i];
		iperm[j] = i;
		order[i] = table->order[j];
		naxes[i] = table->naxes[j];
		nknots[i] = table->nknots[j];
		knots[i] = table->knots[j];
		extents[i][0] = table->extents[j][0];
		extents[i][1] = table->extents[j][1];
	}
	
	/* Compute new strides */
	total_size = 1;
	for (i = table->ndim-1; i >= 0; i--) {
		strides[i] = total_size;
		total_size *= naxes[i];
	}
	
	/* Re-order coefficient array */
	coefficients = malloc(sizeof(table->coefficients[0])*total_size);
	for (pos = 0; pos < total_size; pos++) {
		npos = 0;
		/* Multiply index of point in old shape by new stride */
		for (i = 0; i < table->ndim; i++)
			npos += (pos / table->strides[i] % table->naxes[i])*strides[iperm[i]];
		coefficients[npos] = table->coefficients[pos];
	}
	
	free(table->coefficients);
	free(table->order);
	free(table->naxes);
	free(table->strides);
	free(table->nknots);
	free(table->knots);
	free(table->extents[0]);
	free(table->extents);
	
	table->coefficients = coefficients;
	table->order = order;
	table->naxes = naxes;
	table->strides = strides;
	table->nknots = nknots;
	table->knots = knots;
	table->extents = extents;
}

static int
fillfitstable(fitsfile *fits, const struct splinetable *table)
{
	int error = 0;
	int i;
	
	/*
	 * Create the coefficient array with transposed axis
	 * counts, like PyFITS does.
	 */
	{
		long *naxes = malloc(sizeof(long)*table->ndim);
		for (i = 0; i < table->ndim; i++)
			naxes[i] = table->naxes[table->ndim - i - 1];
		
		fits_create_img(fits, FLOAT_IMG, table->ndim, naxes, &error);
		free(naxes);
		
		if (error != 0)
			return (error);
	}
	
	/*
	 * Write coefficient array
	 */
	{
		long *fpixel = malloc(sizeof(long)*table->ndim);
		long arraysize = 1;
		for (i = 0; i < table->ndim; i++) {
			fpixel[i] = 1;
			arraysize *= table->naxes[i];
		}
		
		fits_write_pix(fits, TFLOAT, fpixel, arraysize,
		    table->coefficients, &error);
		free(fpixel);
		
		if (error != 0)
			return (error);
	}
	
	/*
	 * Write required splinetable keywords
	 */
	for (i = 0; i < table->ndim; i++) {
		char name[255];
		sprintf(name,"ORDER%d",i);
		fits_write_key(fits, TINT, name, &table->order[i],
		    NULL, &error);
		if (error != 0)
			return (error);
		sprintf(name,"PERIOD%d",i);
		fits_write_key(fits, TINT, name, &table->periods[i],
		    NULL, &error);
		if (error != 0)
			return (error);
	}
	
	/*
	 * Auxiliary keywords have no imposed type, and were read in as raw records.
	 * Write them back the same way.
	 */
	for (i = 0; i < table->naux; i++) {
		char card[FLEN_CARD];
		snprintf(card, FLEN_CARD, "%-8s= %s", table->aux[i][0], table->aux[i][1]);
		fits_write_record(fits, card, &error);
		if (error != 0)
			return (error);
	}
	
	/*
	 * Write each of the knot vectors in an extension HDU
	 */
	for (i = 0; i < table->ndim; i++) {
		char name[255];
		long fpixel = 1;
		sprintf(name,"KNOTS%d",i);
		fits_create_img(fits, DOUBLE_IMG, 1, &table->nknots[i], &error);
		if (error != 0)
			return (error);
		fits_write_key(fits, TSTRING, "EXTNAME", name,
		    NULL, &error);
		if (error != 0)
			return (error);
		fits_write_pix(fits, TDOUBLE, &fpixel, table->nknots[i],
		    table->knots[i], &error);
		if (error != 0)
			return (error);
	}
	
	/*
	 * Write the boundaries of support to an extension HDU
	 */
	{
		long fpixel = 1;
		long nextents = 2*table->ndim;
		fits_create_img(fits, DOUBLE_IMG, 1, &nextents, &error);
		if (error != 0)
			return (error);
		fits_write_key(fits, TSTRING, "EXTNAME", "EXTENTS", NULL, &error);
		if (error != 0)
			return (error);
		fits_write_pix(fits, TDOUBLE, &fpixel, nextents,
		    table->extents[0], &error);
		if (error != 0)
			return (error);
	}
	
	return (error);
}

static int
parsefitstable(fitsfile *fits, struct splinetable *table)
{
	int error = 0;
	int hdus, type, i, nkeys;
	size_t arraysize;
	long *fpixel;

	fits_get_num_hdus(fits, &hdus, &error);
	fits_movabs_hdu(fits, 1, &type, &error);
	if (error != 0)
		return (error);

	if (type != IMAGE_HDU)
		return (ENOENT);

	/*
	 * Read header information
	 */

	fits_get_img_dim(fits, &table->ndim, &error);
	if (error != 0)
		return (error);
	assert(table->ndim >= 1);

	/*
	 * Read in any auxiliary keywords.
	 */
	nkeys = 0;
	fits_get_hdrspace(fits, &nkeys, NULL, &error);
	if (nkeys > 0) {
		char key[FLEN_KEYWORD], value[FLEN_VALUE];
		int keylen, valuelen;
		table->aux = calloc(sizeof(char**), nkeys);
		i = 0;
		int j = 1;
		for ( ; (i < nkeys) && (j-1 < nkeys); j++) {
			error = 0;
			fits_read_keyn(fits, j, key, value, NULL, &error);
			if (error != 0)
				continue;
			if (strncmp("TYPE", key, 4) == 0 ||
			    strncmp("ORDER", key, 5) == 0 || 
			    strncmp("NAXIS", key, 5) == 0 ||
			    strncmp("BITPIX", key, 6) == 0 ||
			    strncmp("SIMPLE", key, 6) == 0 ||
			    strncmp("PERIOD", key, 6) == 0 ||
			    strncmp("EXTEND", key, 6) == 0 ||
			    strncmp("COMMENT", key, 7) == 0)
				continue;

			keylen = strlen(key) + 1;
			valuelen = strlen(value) + 1;
			table->aux[i] = calloc(sizeof(char*), 2);
			table->aux[i][0] = calloc(sizeof(char), keylen);
			table->aux[i][1] = calloc(sizeof(char), valuelen);
			memcpy(table->aux[i][0], key, keylen);
			memcpy(table->aux[i][1], value, valuelen);
			i++;
		}
		table->aux = realloc(table->aux, i*sizeof(char**));
		table->naux = i;
	} else {
		table->aux = NULL;
		table->naux = 0;
	}

	table->order = malloc(sizeof(table->order[i])*table->ndim);
	fits_read_key(fits, TINT, "ORDER", &table->order[0], NULL, &error);
	if (error != 0) {
		error = 0;
		
		for (i = 0; i < table->ndim; i++) {
			char name[255];
			sprintf(name,"ORDER%d",i);
			fits_read_key(fits, TINT, name, &table->order[i],
			    NULL, &error);
			assert(table->order[i] > 0);
		}
	} else {
		assert(table->order[0] > 0);
		for (i = 1; i < table->ndim; i++)
			table->order[i] = table->order[0];
	}

	if (error != 0)
		return (error);

	error = 0;
	table->periods = malloc(sizeof(table->periods[i])*table->ndim);
	for (i = 0; i < table->ndim; i++) {
		char name[255];
		sprintf(name,"PERIOD%d",i);
		fits_read_key(fits, TDOUBLE, name, &table->periods[i], NULL, &error);
		/*
		 * If the PERIOD keys cannot be read, just interpret this as
		 * a non-periodic table.
		 */
		if (error != 0) {
			table->periods[i] = 0;
			error = 0;
		}
	}

	/*
	 * Read the coefficient table
	 */

	table->naxes = malloc(sizeof(long)*table->ndim);
	fits_get_img_size(fits, table->ndim, table->naxes, &error);

	/*
	 * FITS multidimensional arrays are stored as FORTRAN arrays,
	 * not C arrays, so we need to swizzle the matrix into being
	 * a C array. Or we should. Instead, PyFITS, which writes these
	 * files, writes a C array, but with the axis counts transposed.
	 * Fix it.
	 */
	{
		long *naxestmp = malloc(sizeof(long)*table->ndim);
		for (i = 0; i < table->ndim; i++)
			naxestmp[i] = table->naxes[table->ndim - i - 1];

		free(table->naxes);
		table->naxes = naxestmp;
	}

	/* Compute the total array size and the strides into each dimension */
	table->strides = malloc(sizeof(unsigned long)*table->ndim);
	table->strides[table->ndim - 1] = arraysize = 1;
	for (i = table->ndim-1; i >= 0; i--) {
		arraysize *= table->naxes[i];
		if (i > 0)
			table->strides[i-1] = arraysize;
	}
	table->coefficients = malloc(sizeof(float)*arraysize);

	fpixel = malloc(sizeof(long)*table->ndim);
	for (i = 0; i < table->ndim; i++)
		fpixel[i] = 1;

	fits_read_pix(fits, TFLOAT, fpixel, arraysize, NULL,
	    table->coefficients, NULL, &error);

	free(fpixel);

	if (error != 0) {
		fprintf(stderr, "Error reading table coefficients\n");
		splinetable_free(table);
		return (error);
	}

	/*
	 * Read the knot vectors, which are stored one each in extension
	 * HDUs
	 */

	table->knots = malloc(sizeof(table->knots[0])*table->ndim);
	for (i = 0; i < table->ndim; i++)
		table->knots[i] = NULL;
	table->nknots = malloc(sizeof(table->nknots[0])*table->ndim);

	for (i = 0; i < table->ndim; i++) {
		char hduname[255];
		long fpix = 1;
		double *knot_scratch;
		sprintf(hduname,"KNOTS%d",i);

		fits_movnam_hdu(fits, IMAGE_HDU, hduname, 0, &error);
		fits_get_img_size(fits, 1, &table->nknots[i], &error);
		if (error != 0) {
			fprintf(stderr, "Error reading knot vector %d\n", i);
			break;
		}

		/* 
		 * Allow spline evaluations to run off the ends of the
		 * knot field without segfaulting.
		 */
		knot_scratch = calloc(table->nknots[i]
		    +2*table->order[i], sizeof(double));
		table->knots[i] = knot_scratch + table->order[i];
		fits_read_pix(fits, TDOUBLE, &fpix, table->nknots[i], NULL,
		    table->knots[i], NULL, &error);
	}
	if (error != 0) {
		splinetable_free(table);
		return (error);
	}

	/*
	 * Read the axes extents, stored in a single extension HDU.
	 */

	table->extents = malloc(sizeof(double*) * table->ndim);
	table->extents[0] = malloc(sizeof(double) * 2 * table->ndim);

	for (i = 1; i < table->ndim; i++) {
		table->extents[i] = &table->extents[0][2*i];
	}

	long n_extents = 0;
	long fpix = 1;
	int ext_error = 0;
	fits_movnam_hdu(fits, IMAGE_HDU, "EXTENTS", 0, &ext_error);
	fits_get_img_size(fits, 1, &n_extents, &ext_error);
	if (n_extents != 2*table->ndim)
		ext_error = 1;

	if (ext_error != 0) { /* No extents. Make up some reasonable ones. */
		for (i = 0; i < table->ndim; i++) {
			table->extents[i][0] = table->knots[i][table->order[i]];
			table->extents[i][1] = table->knots[i][table->nknots[i]
			    - table->order[i] - 1];
		}
	} else {
		fits_read_pix(fits, TDOUBLE, &fpix, n_extents, NULL,
		     table->extents[0], NULL, &ext_error);
	}

	return (error);
}

char *
splinetable_get_key(const struct splinetable *table, const char *key)
{
	int i = 0;
	char *value = NULL;
	
	for ( ; i < table->naux; i++) {
		if (strcmp(key, table->aux[i][0]) == 0) {
			value = table->aux[i][1];
		}
	}

	return (value);
}

int
splinetable_read_key(const struct splinetable *table, splinetable_dtype type,
    const char *key, void *result)
{
	int error = 0;
	char *value = splinetable_get_key(table, key);

	if (!value)
		return (-1);

	switch (type) {
		case SPLINETABLE_INT:
			ffc2i(value, (long*)result, &error);
			break;
		case SPLINETABLE_DOUBLE:
			ffc2d(value, (double*)result, &error);
			break;
		default:
			error = BAD_DATATYPE;
	}

	if (error != 0)
		return (-1);
	else
		return (0);
	
}
