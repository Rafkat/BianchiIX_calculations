//
// Created by rafkat on 8/27/24.
//

#ifndef ODECALCULATIONS_OUTPUT_WRITER_H
#define ODECALCULATIONS_OUTPUT_WRITER_H

#include <fstream>
#include "../nr/nr3.h"
#include "../nr/odeint.h"


inline void write_reverse_anisotropic(Output &out, ofstream &file) {
    for (Int i = out.count - 1; i > 0; i--)
        file << out.xsave[i] << " " << out.ysave[0][i]
//             << " " << out.ysave[1][i] << " " << out.ysave[2][i] << " " << out.ysave[3][i] << endl;
                 << " " << out.ysave[1][i] << " " << out.ysave[2][i] << " " << out.ysave[3][i]
             << " " << out.ysave[4][i] << " " << out.ysave[5][i] << " " << out.ysave[6][i] << endl;
}

inline void write_anisotropic(Output &out, ofstream &file) {
    for (Int i = 0; i < out.count; i++)
        file << out.xsave[i] << " " << out.ysave[0][i]
//             << " " << out.ysave[1][i] << " " << out.ysave[2][i] << " " << out.ysave[3][i] << endl;
                    << " " << out.ysave[1][i] << " " << out.ysave[2][i] << " " << out.ysave[3][i]
                << " " << out.ysave[4][i] << " " << out.ysave[5][i] << " " << out.ysave[6][i] << endl;
}

inline void write_reverse_isotropic(Output &out, ofstream &file) {
    for (Int i = out.count - 1; i > 0; i--)
        file << out.xsave[i] << " " << out.ysave[0][i]
             << " " << out.ysave[1][i] << " " << out.ysave[2][i] << endl;
}

inline void write_isotropic(Output &out, ofstream &file) {
    for (Int i = 0; i < out.count; i++)
        file << out.xsave[i] << " " << out.ysave[0][i]
             << " " << out.ysave[1][i] << " " << out.ysave[2][i] << endl;
}

inline void write_to_file(Output &out, ofstream &file, bool reverse, bool anisotropic = false) {
    if (reverse) {
        if (anisotropic) {
            write_reverse_anisotropic(out, file);
        } else {
            write_reverse_isotropic(out, file);
        }
    } else {
        if (anisotropic) {
            write_anisotropic(out, file);
        } else {
            write_isotropic(out, file);
        }
    }

}


#endif //ODECALCULATIONS_OUTPUT_WRITER_H
