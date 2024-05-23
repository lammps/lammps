### POD example for Ta

We will fit a potential to the `Ta` training data in the `XYZ` directory, which houses `.xyz` files 
of the training data taken from [the FitSNAP repo](https://github.com/FitSNAP/FitSNAP/tree/master/examples/Ta_XYZ/XYZ)

Fit POD with

    lmp -in in.fitpod

This creates `Ta_coefficients.pod` for the linear model, which we can use to run MD with

    lmp -in in.pod


 
