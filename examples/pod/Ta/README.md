### POD example for Ta

First, obtain the training data (the `XYZ` directory) from [here](https://github.com/FitSNAP/FitSNAP/tree/master/examples/Ta_XYZ/XYZ)

Fit POD with

    lmp < in.podfit

This creates `coefficients.txt`, which we can use to run MD with

    lmp < in.podrun


 