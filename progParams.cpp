/// A class that will carry all my parameters

#include "progParams.h"


/*class programParams
 * {
 * public:
 * programParams();
 * int initialize(char *filename);
 * int matrixsize;
 * int stepnumber;
 * char *finfile;
 * char *outfile;
 * int initialconfig;
 * int measure;
 * char *inifile;
 * double gD2;
 * double gD22;
 * double gD4;
 * double wmoveA;
 * double wmoveM;
 * int Type;
 *
 * }; */

programParams::programParams()
{
    // some hardcoded default values, for lazy people
    matrixsize = 2;
    wmoveA     = 1. / 6. / pow(matrixsize, 1.5);
    stepnumber = 1000;
    outfile    = (char *) calloc(150, sizeof(char));
    strcpy(outfile, "default_output.txt");
    initialconfig = 1;
    measure       = 0;
    Type    = 13;
    gD2     = 1;
    gD22    = 0;
    gD4     = 0;
    moveT   = 0;
    inifile = (char *) calloc(150, sizeof(char));
    finfile = (char *) calloc(150, sizeof(char));
    strcpy(finfile, "default_finalmatrix.txt");
}

int
programParams::initialize(char * filename)
{
    FILE * fin;
    int j, iniI = 0, buf;
    char * value;
    char buffer[150];
    double p;


    value = (char *) calloc(150, sizeof(char));
    //	buffer= (char*) calloc (30,sizeof(char));

    fin = fopen(filename, "r");
    if (!fin) {
        printf("Sorry but your input file does not exist");
        return 1;
    }


    // while there are more lines we keep going
    while (j != -1) {
        j = fscanf(fin, "%s  %s \n", value, buffer);

        if (strcmp(value, "matrixsize") == 0) {
            buf        = atoi(buffer);
            matrixsize = buf;
            wmoveA     = 1. / 6. / pow(matrixsize, 1.5);


            if (DEBUG) printf(" %d %s \n", buf, value);
        } else if (strcmp(value, "Type") == 0) { /// I am assuming that the type is (p,q) written as pq.
            /// this only works as long as p,q<10 but I think it's a safe bet
            buf  = atoi(buffer);
            Type = buf;

            if (DEBUG) printf(" %d %s \n", buf, value);
        } else if (strcmp(value, "steps") == 0) {
            buf        = atoi(buffer);
            stepnumber = buf;

            if (DEBUG) printf(" %d %s \n", buf, value);
        } else if (strcmp(value, "outputfile") == 0) {
            strcpy(outfile, buffer);

            if (DEBUG) printf(" %s %s \n", outfile, value);
        } else if (strcmp(value, "initialconfig") == 0) {
            buf = atoi(buffer);
            initialconfig = buf;

            if (DEBUG) printf(" %d %s \n", buf, value);
        } else if (strcmp(value, "measurement") == 0) {
            buf     = atoi(buffer);
            measure = buf;

            if (DEBUG) printf(" %d %s \n", buf, value);
        } else if (strcmp(value, "moveT") == 0) {
            buf   = atoi(buffer);
            moveT = buf;

            if (DEBUG) printf(" %d %s \n", buf, value);
        } else if (strcmp(value, "initalfile") == 0) {
            strcpy(inifile, buffer);
            iniI = 1;

            if (DEBUG) printf(" %s %s \n", inifile, value);
        } else if (strcmp(value, "movesradiusA") == 0) {
            p      = atof(buffer);
            wmoveA = p;

            if (DEBUG) printf(" %e %s \n", p, value);
        } else if (strcmp(value, "couplingD4") == 0) {
            p   = atof(buffer);
            gD4 = p;

            if (DEBUG) printf(" %e %s \n", p, value);
        } else if (strcmp(value, "couplingD2") == 0) {
            p   = atof(buffer);
            gD2 = p;

            if (DEBUG) printf(" %e %s \n", p, value);
        } else if (strcmp(value, "couplingD22") == 0) {
            p    = atof(buffer);
            gD22 = p;

            if (DEBUG) printf(" %e %s \n", p, value);
        } else if (strcmp(value, "finalmatrixfile") == 0) {
            strcpy(finfile, buffer);
            iniI = 1;

            if (DEBUG) printf(" %s %s \n", finfile, value);
        } else {
            printf("Are you sure %s is a valid option? \n", value);
        }
    }

    if (iniI == 0 && initialconfig == 4) {
        printf("You need to tell me an initial file.");
        return 1;
    }


    fclose(fin);
    free(value);
    // free(buffer);
    return 0;
} // programParams::initialize

void
programParams::announce()
{
    std::cout << "You are simulating a " << matrixsize << "x" << matrixsize << " matrix" << std::endl;
    std::cout << "Of type (p,q)=pq " << Type << " " << std::endl;
    std::cout << "There will be " << stepnumber << " sweeps of " << matrixsize * 4 << " attempted MC moves"
              << std::endl;
    std::cout << "The  D^2 coupling is  " << gD2 << "The  (Tr(D^2)^2 coupling is  " << gD22
              << "and the D^4 coupling is  " << gD4 << std::endl;
    if (measure == 0) {
        std::cout << "You are measuring set 0, that means only the action, kind of obsolete see action_monitor.txt"
                  << std::endl;
    } else if (measure == 1) {
        std::cout << "You are measuring set 1, that means the action and the determinants for the mi" << std::endl;
    } else if (measure == 2) {
        std::cout << "You are measuring set 2, that means the eigenvalues of D" << std::endl;
    } else if (measure == 3) {
        std::cout << "You are measuring set 3, that means all elements of D" << std::endl;
    }
    if (wmoveA == 0.) std::cout << " The additive move distance will be determined dynamically" << std::endl;
    else std::cout << "And a additive move distance of " << wmoveA << std::endl;

    std::cout << "And the resulting data will be gathered in " << outfile << std::endl;
}
