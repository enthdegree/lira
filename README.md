# lira
Lovasz Integer Relation Algorithm: find a basis for the lattice of integer vectors orthogonal to the one specified  

##

The Lovasz Integer Relation Algorithm from *Finding Integer Relations in Polynomial Time* by Hastad, Just, Lagrias and Schnorr.

Find a basis for the collection of integer vectors orthogonal to `v_x`.  Paraphrasing the paper:
> On input `v_x` in `Z^n`, this algorithm finds a basis of the `n-1`-dimensional integer lattice `L` whose complement is the integer-span of `v_x`. The algorithm and output have nice properties by Theorem 6.1.


## Input 

  `v_x` = nonzero integer `n`-vector

## Output

 `mtx_C` = `n`*`(n-1)` integer matrix, each column an integer vector orthogonal
         to `v_x`, with properties according to 

## Usage Example
    
    >> v_x = [1;4;2;1];
    >> mtx_C = lira(v_x)
    
    ans =
  
        -1    -1     1
         1     0     0
        -1     1     0
        -1    -1    -1
    
    >> v_x'*mtx_C 
    
    ans = 

        0     0     0

