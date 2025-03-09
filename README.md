# Solution for Nintendo HireMe Challenge (0.07ms~ per solution)

The solve process was relatively straightforward. 

1. **Linear Transformation (L)**:  
   First, I recognized the linear transformation `L`, which uses the diffusion matrix.

2. **Non-linear and Non-invertible Transformations (H and S)**:  
   Then, there are the non-linear, non-invertible transformations `H` and `S`.

3. **Forward Transformation**:  
   The complete forward transformation can be written as:  
   `I(S * L)^256 * H`

4. **Inversion of `L`**:  
   `L^-1` can be implemented simply by inverting the matrix. Since it is a small 32x32 matrix, it is fast to compute.

5. **Inversion of `H` and `S`**:  
   Inverting `H^-1` and `S^-1` is trickier because, when there are multiple possibilities, we need to branch out to find all possible values.

6. **Handling S Compatibility**:  
   Another key aspect was recognizing that `S` cannot produce all possible values. Therefore, we need to handle cases where we choose a value incompatible with `S`.
