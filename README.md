This is part of the code used to implement several computation models for V1 and V2 developed by Haruo Hosoya and Aapo Hyvarinen [1,2].
  * Hierarchical sparse coding model for V2 [1]
  * Pooling model for V1 complex cells [2]

The code relies on overcomplete ICA algorithm based on score matching [3].  This needs to be downloaded separatedly.
  
To use the code, run first initpath.m.  

The entry scripts are:
  * v1_example.m for V1 models (pooling based on PCA dimension reduction [1])
  * v2_example.m for V2 models (hierarchical sparse coding [2])

If you publish a paper based on this code, please cite [1], [2], or any following conference/journal publication.

[1] Haruo Hosoya, Aapo Hyv?rinen. Learning Visual Spatial Pooling by Strong PCA Dimension Reduction. Neural Computation, 82:1-16, 2016.

[2] Hosoya H, Hyv?rinen A. A Hierarchical Statistical Model of Natural Images Explains Tuning Properties in V2. Journal of Neuroscience, 35:10412-10428, 2015.

[3] https://github.com/HaruoHosoya/smica
