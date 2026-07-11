### POD example for InP

We will fit a potential to the `InP` training data in the `XYZ` directory, which houses `.xyz` files 

Please download the training data from [the repo](https://github.com/cesmix-mit/pod-examples/tree/main/JCP2023_InP/XYZ)

Fit POD with

    lmp -in in.fitpod

This creates `InP_coefficients.pod` for the linear model, which we can use to run MD with

    lmp -in in.pod


 
