# Spatially-Coupled-QLDPC-Codes
Software that optimizes SC-HGP codes in [https://arxiv.org/abs/2305.00137](https://arxiv.org/pdf/2305.00137v4)

# Code Construction of SC-HGP codes
Run `sc_hgp_grade_ao_construction.cpp` to construct SC-HGP codes of desired parameters. Please change the corresponding parameters in line 14-23:
```
    int r1=3;
    int r2=3; 
    int n1=7;
    int n2=7; // (r1, n1) and (r2, n2) represent the dimensions of the component matrices that form the underlying hypergraph product codes
    int m1=2;
    int m2=2; // m1 and m2 are memories of 2D-SC codes
    double w=25; // this is the weight of cycles 8, assuming cycles 6 are of weight 1
    int L1=10; 
    int L2=10; // L1 and L2 are coupling lengths of 2D-SC codes
    double step=0.02; // this is the step size of the gradient descent algorithm used in GRADE
```

Below is an example of output of AO optimization of partitioning matrices initialized by distribution optimized by GRADE:
```
P_init: // this is the partitioning matrices pair initialized by distribution optimizaed by GRADE
0,8,3,0,5,2,2,6,0,1,6,7,4,2,1,8,6,8,6,0,7,
3,0,7,6,1,1,5,2,8,8,0,6,0,8,0,4,7,6,6,2,2,
3,38,
4,65,
compute the total start number
n_cur start: // these are numbers of cycles 4, cycles 6, cycles 8 of the random code initialized by distribution optimized by GRADE (e.g. Code 4 in the paper)
50,70,10688,
N_cur start=12438
P_end:
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,
finished P_gd ao
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,

P_init:
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,
0,38,
0,32,
compute the total start number
n_cur start:
0,0,1596,
N_cur start=1596
P_end:
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,
check P_gd ao
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,

P_init:
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,
0,38,
0,32,
compute the total start number
n_cur start: // these are numbers of cycles 4, cycles 6, cycles 8 of the code optimized by GRADE-AO
0,0,1596, 
N_cur start=1596
P_end:
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,
check P1_gd ao // this is the code optimized by GRADE-AO (e.g., Codes 1, 3, 5-7 in the paper)
3,0,6,2,2,8,2,8,0,5,6,7,2,3,2,8,1,8,6,0,7,
8,0,1,6,2,4,3,2,3,8,6,1,0,8,0,5,3,1,6,8,2,
```

Below is an example of output of AO optimization of partitioning matrices initialized by uniform distribution:
```
P_init: // this is the partitioning matrices pair initialized by uniform distribution
1,6,2,3,0,0,8,4,7,3,1,5,1,6,4,8,2,7,5,7,4,
1,5,6,7,4,2,5,8,8,0,4,2,4,0,7,1,6,1,7,3,3,
3,56,
8,59,
compute the total start number
n_cur start: // these are numbers of cycles 4, cycles 6, cycles 8 of the random code initialized by uniform distribution (e.g., Code 2 in the paper)
40,110,19216,
N_cur start=21966
P_end:
2,5,8,8,1,3,3,0,6,3,2,5,1,6,4,8,2,3,7,7,4,
8,3,2,1,3,4,1,0,8,0,4,2,4,0,7,2,6,1,7,2,3,
finished P_uni ao
2,5,8,8,1,3,3,0,6,3,2,5,1,6,4,8,2,3,7,7,4,
8,3,2,1,3,4,1,0,8,0,4,2,4,0,7,2,6,1,7,2,3,

P_init:
2,5,8,8,1,3,3,0,6,3,2,5,1,6,4,8,2,3,7,7,4,
8,3,2,1,3,4,1,0,8,0,4,2,4,0,7,2,6,1,7,2,3,
0,57,
1,58,
compute the total start number
n_cur start:
0,10,6388,
N_cur start=6638
P_end:
2,0,8,8,1,3,3,0,6,3,2,5,1,6,6,8,2,6,7,2,4,
4,3,2,1,3,4,5,0,8,0,4,2,6,0,7,0,6,2,7,2,4,
Check P_uni ao
2,0,8,8,1,3,3,0,6,3,2,5,1,6,6,8,2,6,7,2,4,
4,3,2,1,3,4,5,0,8,0,4,2,6,0,7,0,6,2,7,2,4,

P_init:
2,0,8,8,1,3,3,0,6,3,2,5,1,6,6,8,2,6,7,2,4,
4,3,2,1,3,4,5,0,8,0,4,2,6,0,7,0,6,2,7,2,4,
0,46,
0,52,
compute the total start number
n_cur start: // these are numbers of cycles 4, cycles 6, cycles 8 of the code optimized by AO only 
0,0,2388,
N_cur start=2388
P_end:
2,0,8,8,1,3,3,0,6,3,2,5,1,6,6,8,2,6,7,2,4,
4,3,2,3,6,6,5,0,8,0,4,2,6,0,7,0,6,2,7,2,4,
Check P1_uni ao // this is the code optimized by AO only (e.g., Code 2 in the paper)
2,0,8,8,1,3,3,0,6,3,2,5,1,6,6,8,2,6,7,2,4,
4,3,2,3,6,6,5,0,8,0,4,2,6,0,7,0,6,2,7,2,4,
```

# Write codes to .txt file
If you have a pair of partitioning matrices obtained from `sc_hgp_grade_ao_construction.cpp`, run `sc_generate_code.cpp` to generate the 2D SC-HGP code and write it to .txt file.

Please replace line 34-44 with the paritioning matrices and the corresponding parameters. For example, if you obtained the following matrix from `sc_hgp_grade_ao_construction.cpp`: 
```
8,3,1,7,7,2,1,2,4,8,1,6,0,4,7,5,5,4,3,6,0,
5,0,4,8,1,6,1,6,8,1,4,5,4,3,3,7,2,7,7,0,2,
```
You just need to assign this matrix to P in `sc_generate_code.cpp` as follows:
```
    vector<vector<int> > Code2_3_7{{8,3,1,7,7,2,1,2,4,8,1,6,0,4,7,5,5,4,3,6,0},{5,0,4,8,1,6,1,6,8,1,4,5,4,3,3,7,2,7,7,0,2}};//30,150,17804
    vector<vector<int> > P=Code2_3_7;
    string filename = "Code2_3_7.txt";# this is the .txt file name of your code

    int r1=3;
    int r2=3;
    int n1=7;
    int n2=7;
    int m1=2;
    int m2=2;
    int L1=10;
    int L2=10;
```
# Attribution
If you use this software in your research, please cite it the original paper as follows:
```
@article{yang2023quantum,
  title={Quantum spatially-coupled codes},
  author={Yang, Siyi and Calderbank, Robert},
  journal={arXiv preprint arXiv:2305.00137},
  year={2023}
}
```
