/*
This project executes the image reconstruction for SPCI measurement data.
To do so, two different algorithms can be used:
- maximum likelihood expectation maximization (ML-EM),
- origin ensemble (OE).

Author: Dominik Kornek <dominik.kornek@gmail.com>
*/


/*
TO-DO:

1) create projection matrix:
    1.1 read measurement data (histograms)
    1.2 get N_dbc (Number of events of detector pair d/c in bin b)
    1.3 model real activity distribution provided by the measurements A_v
    1.3 calculate matrix elements p_dcbv

2) implement ML-EM algorithm:
    2.1 create a-priori activity distribution (e.g. homogeneous distribution)
    2.2 implement formula

3) implement OE algorithm:
    3.1 ?

4) analysis:
    4.1 create plots for every iteration
    4.2 update graph of maximization (for ML-EM (log?))
    4.3 generate file with image properties (e.g. SNR, resolution ...)
*/


int main(){
    return 0;
}
