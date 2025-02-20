# Spatially-Coupled-QLDPC-Codes
Software that optimizes SC-HGP codes in https://arxiv.org/abs/2305.00137

# Write codes to .txt file
If you have a pair of partitioning matrices at hand, run sc_generate_code.cpp to generate the 2D SC-HGP code and write it to .txt file.

Please replace line 34-44 with the paritioning matrices and the corresponding parameters, the example is as follows:
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

# Attributions
If you use this software in your research please cite as follows:
```
@software{SpatiallyCoupledQLDPCCodes,
  author = Siyi Yang,
  title = Spatially Coupled QLDPC Codes,
  year = {2025},
  url = {https://github.com/SiyiPriscaYang/Spatially-Coupled-QLDPC-Codes},
}
```
