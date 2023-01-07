#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main( int argc, char **argv){
	char binfilename[512];

	FILE *pFile;
	double value;

	int pos = 2;
	int actual_pos_x;
	int actual_pos_y;
	int length_x = 1000;
	
	// PARAMETERS
	char c;
	while ((c = getopt (argc, argv, "X:")) != -1)
	         switch (c)
	           {
	           case 'X':
	        	   length_x = atoi(optarg);
	        	   //fprintf( stderr, "Set Remodeling Iterations to %i\n", maxit);
	        	   break;
	           }

	// READ FILES
	for( int f=optind; f<argc; f++){
		pFile = fopen( argv[f],"rb");
		
		actual_pos_x = 0;
		actual_pos_y = 0;
		
		while( fread (&value , 1 , sizeof(value) , pFile ))
		//while( fread (&value , 1 , sizeof(value) , pFile ) && actual_pos<pos)
		{
			printf( "%i %i %i %lf\n", f, actual_pos_x, actual_pos_y, value);
			actual_pos_x++;
			if( actual_pos_x == length_x){
				actual_pos_x=0;
				actual_pos_y++;
			}

		}
		printf( "\n");
		//printf( "%lf\n", value);

		fclose(pFile);
	}
}
