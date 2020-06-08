#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <photospline/splinetable.h>

static void usage() {
	fprintf(stderr,"dumpsplinefits <path> [coeff]\n");
	exit(1);
}

int main(int argc, char **argv) {
	struct splinetable table;
	int i, j;

	if (argc < 2)
		usage();

	readsplinefitstable(argv[1], &table);

	printf("NDim: %d\n",table.ndim);

	for (i = 0; i < table.ndim; i++) {
		printf("Dimension %d\n",i);
		printf("\tOrder: %d\n",table.order[i]);
		printf("\tN Knots: %d",(int)table.nknots[i]);
		for (j = 0; j < table.nknots[i]; j++) {
			if (j % 6 == 0)
				printf("\n\t\t");
			printf("%lf ",table.knots[i][j]);
		}
		printf("\n");
		printf("\tPeriod: %lf\n",table.periods[i]);
		printf("\tN Splines: %d\n",(int)table.naxes[i]);
	}
	
	if (argv[2] != NULL && strcmp(argv[2],"coeff") == 0) {
		j = table.naxes[0];
		for (i = 1; i < table.ndim; i++)
			j *= table.naxes[i];

		printf("\nCoefficients:");
		for (i = 0; i < j; i++) {
			if (i % 6 == 0)
				printf("\n");
			printf("%e ",table.coefficients[i]);
		}
		printf("\n");
	}

	return 0;
}

