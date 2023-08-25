.. _hooke_model

For the *hooke* model, the normal, elastic component of force acting
on particle *i* due to contact with particle *j* is given by:

.. math::

   \mathbf{F}_{ne, Hooke} = k_n \delta_{ij} \mathbf{n}

Where 
:math:`\delta_{ij} = R_i + R_j - \|\mathbf{r}_{ij}\|` is the particle overlap, 
:math:`R_i, R_j` are the particle radii, 
:math:`\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j` is the vector separating the two
particle centers (note the i-j ordering so that :math:`\mathbf{F}_{ne}` is
positive for repulsion), and 
:math:`\mathbf{n} = \frac{\mathbf{r}_{ij}}{\|\mathbf{r}_{ij}\|}`.  

The units of the spring constant :math:`k_n` are
*force*\ /\ *distance*, or equivalently *mass*\ /*time\^2*.

.. _hertz_model
For the *hertz* model, the normal component of force is given by:

.. math::

   \mathbf{F}_{ne, Hertz} = k_n R_{eff}^{1/2}\delta_{ij}^{3/2} \mathbf{n}

Where 
:math:`\delta_{ij} = R_i + R_j - \|\mathbf{r}_{ij}\|` is the particle overlap, 
:math:`R_i, R_j` are the particle radii, 
:math:`R_{eff} = \frac{R_i R_j}{R_i + R_j}` is the effective radius
:math:`\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j` is the vector separating the two
particle centers (note the i-j ordering so that :math:`\mathbf{F}_{ne}` is
positive for repulsion), and 
:math:`\mathbf{n} = \frac{\mathbf{r}_{ij}}{\|\mathbf{r}_{ij}\|}`. 

The units of the spring constant :math:`k_n` are *force*\ /\ *length*\ \^2, or
equivalently *pressure*\ .

