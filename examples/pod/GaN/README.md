### POD example for GaN

First, obtain the training data (the `XYZ` directory) from [here](https://github.com/cesmix-mit/TrainingData/tree/main/data/GaN)

Fit POD with

    lmp < in.podfit

This creates `coefficients.txt`, which we can use to run MD with

    lmp < in.run


