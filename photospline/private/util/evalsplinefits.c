#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

static void usage() {
	fprintf(stderr,"evalsplinefits <path> x1 [x2_i-x2_f] x3 ...\n");
	exit(1);
}

#define TIMING
#define GRADIENTS
#define SAMPLES 1000

int main(int argc, char **argv) {
	struct splinetable table;
	int i, j, iterdim;
	double *x, x_i, x_f, value, *eval1, *eval2;
	double *randomseq;
	int *centers;
	struct timeval tp1, tp2;

	if (argc < 2)
		usage();

	gettimeofday(&tp1, NULL);
	readsplinefitstable(argv[1], &table);
	gettimeofday(&tp2, NULL);

    #ifdef TIMING
	tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
	printf("Time to open table: %ld.%06ld seconds\n", 
	    (long)(tp2.tv_sec - tp1.tv_sec), (long)(tp2.tv_usec - tp1.tv_usec));
    #endif

	if (argc < 2+table.ndim)
		usage();

	x = malloc(table.ndim*sizeof(double));
	centers = malloc(table.ndim*sizeof(int));
	iterdim = -1;
	for (i = 0; i < table.ndim; i++) { 
		sscanf(argv[i+2],"%lf",&x[i]);
		if (strchr(argv[i+2], '-') != NULL) {
			sscanf(strchr(argv[i+2], '-') + 1, "%lf", &x_f);
			iterdim = i;
		}
	}
	gettimeofday(&tp1, NULL);
	if (tablesearchcenters(&table, x, centers)) {
		printf("Table search centers failed for coordinates: [");
		for (i = 0; i < table.ndim; i++)
			printf("%lf, ", x[i]);
		printf("\b\b]\n");
		free(x);
		free(centers);
		splinetable_free(&table);
		return (-1);
	}
	value = ndsplineeval(&table, x, centers, 0);
	gettimeofday(&tp2, NULL);

    #ifdef TIMING
	tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
	printf("Cold cache evaluation time: %ld microseconds\n",
	    (long)(tp2.tv_usec - tp1.tv_usec));

	gettimeofday(&tp1, NULL);
	for (i = 0; i < SAMPLES; i++)
		value = ndsplineeval(&table, x, centers, 0);
	gettimeofday(&tp2, NULL);

	tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
	printf("Warm cache evaluation time: %ld.%02ld microseconds\n",
	    (long)((tp2.tv_usec - tp1.tv_usec)/SAMPLES),
	    (long)((((tp2.tv_usec - tp1.tv_usec)%SAMPLES) * 100) /SAMPLES));
    #endif

	printf("NDim: %d\n",table.ndim);
	printf("Order: %d\n",table.order[0]);

	printf("Value: %lf\n",value);
	printf("e^(value): %e\n",exp(value));

    #ifdef TIMING
	value = ndsplineeval_linalg(&table, x, centers, 0);
	gettimeofday(&tp1, NULL);
	for (i = 0; i < SAMPLES; i++)
		value = ndsplineeval_linalg(&table, x, centers, 0);
	gettimeofday(&tp2, NULL);

	tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
	printf("Warm cache evaluation time (BLAS): %ld.%02ld microseconds\n",
	    (long)((tp2.tv_usec - tp1.tv_usec)/SAMPLES),
	    (long)((((tp2.tv_usec - tp1.tv_usec)%SAMPLES) * 100) /SAMPLES));
	printf("Value (BLAS): %lf\n",value);
    #endif
	
	
    #ifdef GRADIENTS
	eval1 = calloc(table.ndim + 1, sizeof(double));
	eval2 = calloc(table.ndim + 1, sizeof(double));
	
	eval1[0] = ndsplineeval(&table, x, centers, 0);
	for (j = 0; j < table.ndim; j++)
		eval1[1+j] = ndsplineeval(&table, x, centers,
		    (1 << j));
	ndsplineeval_gradient(&table, x, centers, eval2);
	
	printf("Sequential gradient:");
	for (j = 0; j < table.ndim+1; j++)
		printf(" %- .2e", eval1[j]);
	printf("\n");
	printf("Vectorized gradient:");
	for (j = 0; j < table.ndim+1; j++)
		printf(" %- .2e", eval2[j]);
	printf("\n");
	
    #endif

	if (iterdim >= 0) {
		double min, max; 
		randomseq = calloc(SAMPLES, sizeof(double));
		x_i = x[iterdim];

		min = x_f - x_i; max = 0;
		for (i = 0; i < SAMPLES; i++) {
			randomseq[i] = ((double)rand())/(double)RAND_MAX * (x_f - x_i);
			if (randomseq[i] < min) min = randomseq[i];
			if (randomseq[i] > max) max = randomseq[i];
		}
		printf("Evaluating %d samples from %lf to %lf on axis %d\n",
		    SAMPLES, x_i + min, x_i + max, iterdim);
		gettimeofday(&tp1, NULL);
		for (i = 0; i < SAMPLES; i++) {
			x[iterdim] = x_i + randomseq[i];
			if (tablesearchcenters(&table, x, centers)) {
				printf("Search table centers failed for %lf\n", x[iterdim]);
				break;
			}
			value = ndsplineeval(&table, x, centers, 0);
			if (!isfinite(value))
				fprintf(stderr, "Encountered a NaN!\n");
		}
		gettimeofday(&tp2, NULL);
	    #ifdef TIMING
		tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
		printf("Multiple evaluation time: %ld.%02ld microseconds\n",
		    (long)((tp2.tv_usec - tp1.tv_usec)/SAMPLES),
		    (long)((((tp2.tv_usec - tp1.tv_usec)%SAMPLES) * 100) /SAMPLES));
	    #endif
	
	    #ifdef GRADIENTS
		
		gettimeofday(&tp1, NULL);
		for (i = 0; i < SAMPLES; i++) {
			x[iterdim] = x_i + randomseq[i];
			if (tablesearchcenters(&table, x, centers)) {
				printf("Search table centers failed for %lf\n", x[iterdim]);
				break;
			}
			eval1[0] = ndsplineeval(&table, x, centers, 0);
			for (j = 0; j < table.ndim; j++)
				eval1[1+j] = ndsplineeval(&table, x, centers,
				    (1 << j));
		}
		gettimeofday(&tp2, NULL);
		tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
		printf("Sequential gradient evaluation time: %ld.%02ld microseconds\n",
		    (long)((tp2.tv_usec - tp1.tv_usec)/SAMPLES),
		    (long)((((tp2.tv_usec - tp1.tv_usec)%SAMPLES) * 100) /SAMPLES));
		
		
		gettimeofday(&tp1, NULL);
		for (i = 0; i < SAMPLES; i++) {
			x[iterdim] = x_i + randomseq[i];
			if (tablesearchcenters(&table, x, centers)) {
				printf("Search table centers failed for %lf\n", x[iterdim]);
				break;
			}
			ndsplineeval_gradient(&table, x, centers, eval2);
		}
		gettimeofday(&tp2, NULL);
		tp2.tv_usec += (tp2.tv_sec - tp1.tv_sec)*1e6;
		printf("Combined gradient evaluation time:   %ld.%02ld microseconds\n",
		    (long)((tp2.tv_usec - tp1.tv_usec)/SAMPLES),
		    (long)((((tp2.tv_usec - tp1.tv_usec)%SAMPLES) * 100) /SAMPLES));
		
	    #endif
	
		free(randomseq);
	
	}
	
    #ifdef GRADIENTS
	free(eval1);
	free(eval2);
    #endif

	return 0;
}
