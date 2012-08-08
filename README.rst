====================
G'DAY
====================

GDAY simulates carbon, nitrogen and water cycling between the plant and the soil. The model is structured into three plant pools (foliage, wood and fine roots), four litter pools (above/below metabolic and structural litter) and three soil organic matter (SOM) pools with varying turnover rates (active, slow
and passive).

.. contents:: :local:


Installation
=============

G'DAY model depends on very few other packages, in fact only one I think??

http://www.voidspace.org.uk/python/configobj.html

And replacing this with a small bit of code to read user input files would be
trivial.

Once you have downloaded the source code, or cloned the repository there is 
a simple makefile, e.g. ::

    sudo make install

or ::

    sudo python setup.py install

I need to add some simple scripts (to do list!). But... ::

    from gday import gday as model
    G = model.Gday(cfg_fname)
    G.run_sim()

(Also spin up option, adjusting parameter file...add!).

References
=============
1). Comins, H. N. and McMurtrie, R. E. (1993) Ecological Applications, 3, 666-681.


Link to me
============
https://sites.google.com/site/mdekauwe/