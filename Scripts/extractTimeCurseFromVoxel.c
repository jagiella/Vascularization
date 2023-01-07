#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main( int argc, char **argv)
{
	// options
	/*char name[512];
	strcpy( name, argv[1]);
	fprintf( stderr, "open file %s\n", name);*/

	int X,Y,Z;
	X=Y=Z=0;
	bool OUTPUT = false;

	char c;
	while ((c = getopt(argc, argv, "x:y:z:O")) != -1)
		switch (c) {
		case 'x':
			X=atoi(optarg);
			break;
		case 'y':
			Y=atoi(optarg);
			break;
		case 'z':
			Z=atoi(optarg);
			break;
		case 'O':
			OUTPUT=true;
		}
	fprintf( stderr, "Read voxel (%i,%i,%i)\n", X,Y,Z);

	for( int i=optind; i<argc; i++){
		//Open file
		FILE *file = fopen( argv[i], "rb");
		if (!file)
		{
			fprintf(stderr, "Unable to open file %s", argv[i]);
			return 0;
		}

		// read dimensions x,y,z
		int x=0;
		int y=0;
		int z=0;
		if(!fread( &x, sizeof(int), 1, file)) {fprintf(stderr, "Unexpected end of file\n");}
		if(!fread( &y, sizeof(int), 1, file)) {fprintf(stderr, "Unexpected end of file\n");}
		if(!fread( &z, sizeof(int), 1, file)) {fprintf(stderr, "Unexpected end of file\n");}
		fprintf( stderr, "Read file %s (data size: %ix%ix%i)\n",argv[i],x,y,z);

		// read values
		float c=0;
		int ix,iy,iz;
		ix=iy=iz=0;
		while( fread( &c, sizeof(float), 1, file) != 0){
			if(ix==X&&iy==Y&&iz==Z)
				fprintf( stdout, "%f\n", c);
			if(OUTPUT)
				fprintf( stderr, "%f\n", c);
			ix++;
			if(ix==x){
				iy++;
				if(iy==y)
					iz++;
			}
		}

		fclose(file);
	}
}
