// #include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "output_vtu_foreach.h"
scalar f[];

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(ferr, "Usage: %s t_real tsnap output_dir\n", argv[0]);
        return 1;
    }
    double t_real = atof(argv[1]);
    double tsnap = atof(argv[2]);
    const char *output_dir = argv[3];

    char nameOut[80];
    sprintf(nameOut, "intermediate/snapshot-%5.4f", t_real);
    
    if (restore(file=nameOut)) {
        char vtkOut[80];
        sprintf(vtkOut, "%s/out-%5.4f-CC", output_dir, t_real/tsnap*10);
        output_vtu((scalar *){f,p}, (vector *){u}, vtkOut); 
        fprintf(ferr, "snapshot-%5.4f have been converted to out-%5.4f-CC\n", t_real, t_real/tsnap*10);
        return 0;
    } else {
        // fprintf(ferr, "Failed to restore %s\n", nameOut);
        return 1;
    }
}