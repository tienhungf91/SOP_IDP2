#include <algorithm>
#include <numeric>
//#include <cmath>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <sstream>
#include <getopt.h>
#include "dcd.h"
//#include "IO.h"
//#include "cluster.h"
//#include "COM.h"
#include "Rg_shape.h"

using namespace std;
const double PI = 3.1415926535;

////////////////////////////////////////////////////////////
void usage () {
    printf ("   SAXS_2bead_full_res    -t    traj\n");
    printf ("                          -o    SAXS output\n");
    printf ("                          -f    number of frame skipped\n");
    printf ("                          -s    processed every n steps\n");
    printf ("                          -q    input file for coarse-grained amino acid form factor\n");
    printf ("                          -S    sequence\n");
}

////////////////////////////////////////////////
vector< vector<double> > read_form_factor (char * file_q) {
    vector< vector<double> > result;
    std::ifstream f (file_q);
    if (f.is_open()) {
        string line;
        int count = 0;
        while (getline(f, line)) {
            count++;
            if (count <= 3) continue;

            stringstream stream (line);
            double tmp;
            stream >> tmp;
            vector<double> form_factor;
            form_factor.resize(21);
            for (int i = 0; i < 21; i++)
                stream >> form_factor[i];
            result.push_back(form_factor);
        }
    } else {
        printf ("Not able to read amino acid form factor !!!\n");
        exit(0);
    }
    return result;
}

///////////////////////////////////////////////
int get_index_form_factor (int i, char ch) {
    if (i % 2 == 0) return 0;
    else switch (ch) {
        case 'A': return 1;
        case 'C': return 2;
        case 'D': return 3;
        case 'E': return 4;
        case 'F': return 5;
        case 'G': return 6;
        case 'H': return 7;
        case 'I': return 8;
        case 'K': return 9;
        case 'L': return 10;
        case 'M': return 11;
        case 'N': return 12;
        case 'P': return 13;
        case 'Q': return 14;
        case 'R': return 15;
        case 'S': return 16;
        case 'T': return 17;
        case 'V': return 18;
        case 'W': return 19;
        case 'Y': return 20;
        default:
            printf ("Something's wrong here %s!!\n", ch);
            exit(0);
    }
}

///////////////////////////////////////////////////////////
vector<double> calc_saxs (const vector<coordDS> &coord,
                          const vector< vector<double> > &form_factor,
                          char * sequence) {
    vector<double> result;
    result.resize(form_factor.size(), 0);

    for (int i = 0; i < coord.size(); i++) {
        if (i % 2 == 1 and sequence[i/2] == 'G') continue;
        int id1 = get_index_form_factor (i, sequence[i/2]);

        for (int j = 0; j < coord.size(); j++) {
            if (j % 2 == 1 and sequence[j/2] == 'G') continue;
            int id2 = get_index_form_factor (j, sequence[j/2]);
            double r = 0;
            if (i != j) {
                double dx = coord[i].x - coord[j].x;
                double dy = coord[i].y - coord[j].y;
                double dz = coord[i].z - coord[j].z;
                r = sqrt(dx*dx + dy*dy + dz*dz);
            }
            for (int q = 0; q < form_factor.size(); q++) {
                double f1 = form_factor[q][id1];
                if (sequence[i/2] == 'G')
                    f1 += form_factor[q][6];

                double f2 = form_factor[q][id2];
                if (sequence[j/2] == 'G')
                    f2 += form_factor[q][6];

                double F = f1*f2;
                if (q == 0 or i == j) result[q] += F;
                else {
                    double qr = q*0.01*r;
                    result[q] += F * sin(qr) / qr;
                }
            }
        }
    }
    return result;
}

///////////////////////////////////////////////////////
void print_saxs (const vector<double> &saxs,
                 char * file_saxs) {
    std::ofstream f (file_saxs);
    if (f.is_open()) {
        f << "# q (A^-1)   I(q)\n";
        for (int i = 0; i < saxs.size(); i++) {
            f << std::setw(10) << std::setprecision(2) << std::fixed << 0.01*i;
            f << std::setw(16) << std::setprecision(2) << std::fixed << saxs[i] << std::endl;
        }
    }
}

/////////////////////////////////////////////////////
////////            MAIN CODE        ///////////////
/////////////////////////////////////////////////////
int main (int argc, char *argv[]) {
    int opt;
    char * file_dcd = NULL;
    char * file_saxs = "SAXS.dat";
    char * file_q = "aa_ff_2bead.dat";
    char * sequence = NULL;

    long frame_skip = 0;
    int step = 1;

    if (argc == 1) {
        printf ("Need arguments !!!\n");
        usage();
        exit (0);
    }
    do {
        struct option long_options[] =
        {
            {"traj_file",     required_argument,  0,  't'},
            {"frame_skip",    required_argument,  0,  'f'},
            {"step",          required_argument,  0,  's'},
            {"SAXS",          required_argument,  0,  'o'},
            {"form_factor",   required_argument,  0,  'q'},
            {"sequence",      required_argument,  0,  'S'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        opt = getopt_long (argc, argv, "t:f:s:o:q:S:", long_options, &option_index);
        // followed by 1 colon - require an argument; no colon - not argument required

        switch (opt) {
            case 't':
                file_dcd = optarg;
                break;
            case 'f':
                frame_skip = atol (optarg);
                break;
            case 's':
                step = atoi (optarg);
                break;
            case 'o':
                file_saxs = optarg;
                break;
            case 'q':
                file_q = optarg;
                break;
            case 'S':
                sequence = optarg;
                break;
        }
    } while (opt != -1);

    if (sequence == NULL) {
        printf ("Must provide the sequence !!!\n");
        exit(0);
    }

    int i, natoms;
    printf ("Reading amino acid form factor ...\n");
    vector< vector<double> > form_factor = read_form_factor (file_q);

    char * fdx = strtok (file_dcd, " ");
    long total_frame = 0;
    int Ndcd = 0;

    vector<double> saxs;
    saxs.resize (form_factor.size());
    double box[3];

    while (fdx != NULL) {
        ++argv;
        Ndcd++;
        natoms = 0;
        void *v;
        dcdhandle *dcd;

        printf("Processing file %s\n", fdx);

        v = open_dcd_read (fdx, "dcd", &natoms);
        if (!v) {
            fprintf (stderr, "open_dcd_read failed for file %s\n", fdx);
            return 1;
        } else if (2*strlen(sequence) != natoms) {
            printf ("Sequence does not match number of atoms!!\n");
            exit (0);
        }
        dcd = (dcdhandle *)v;

        for (i = 0; i < dcd->nsets; i++) {
            float unitcell[6];
            unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
            unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;

            int rc = read_next_timestep (v, natoms, unitcell);
            if (!rc) {
                if ((i % step != 0) or (Ndcd == 1 && i < frame_skip)) continue;
                if (i % 5000 == 0) printf("Processing frame %d\n", i);

                box[0] = unitcell[0];
                box[1] = unitcell[2];
                box[2] = unitcell[5];
                //printf ("unit cell  %f %f %f\n", box[0], box[1], box[2]);

                vector<coordDS> coord;
                for (int j = 0; j < natoms; j++) {
                    coordDS tmp;
                    tmp.x = dcd->x[j];
                    tmp.y = dcd->y[j];
                    tmp.z = dcd->z[j];
                    coord.push_back (tmp);
                }

                unwrap_mol (coord, box);
                vector<double> saxs_tmp = calc_saxs (coord, form_factor, sequence);

                for (int j = 0; j < saxs.size(); j++)
                    saxs[j] += saxs_tmp[j];

                total_frame += step;
            } else {
                fprintf (stderr, "error in read_next_timestep on frame %d\n", i);
                return 1;
            }
        }

        close_file_read(v);
        fdx = strtok (NULL, " ");
    }

    printf ("Total frames read  %ld\n", total_frame);
    for (int i = 0; i < saxs.size(); i++)
        saxs[i] /= total_frame;

    print_saxs (saxs, file_saxs);
    return 0;
}
