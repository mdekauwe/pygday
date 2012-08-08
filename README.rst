====================
G'DAY
====================

GDAY simulates carbon, nitrogen and water cycling between the plant and the soil. The model is structured into three plant pools (foliage, wood and fine roots), four litter pools (above/below metabolic and structural litter) and three soil organic matter (SOM) pools with varying turnover rates (active, slow
and passive).

The model is coded entirely in 
.. _Python: http://www.python.org/

Key Reference
=============
1). Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO"2 Enrichment; Equilibrium Behavior of Plant-Soil Models. Ecological Applications, 3, 666-681.

.. contents:: :local:

Installation
=============

The G'DAY model depends on very few non-standard python packages, in fact only one I think??!

http://www.voidspace.org.uk/python/configobj.html

And replacing this with a small bit of code to read user input files would be
trivial. I might write that...

Once you have downloaded the source code, or clone the repository (go on...) there is a simple makefile, e.g. ::

    sudo make install

or ::

    sudo python setup.py install

Running the model
=================

I need to add some simple scripts (on the todo list!). But... ::
    
    from gday import gday as model
    G = model.Gday(cfg_fname, spin_up=True)
    G.spin_up_pools()

will spin the model up. Spin up expects a met forcing file with a 1000 yrs of data, how you recycle this is up to you. The model automatically stops once the soil, plant and litter pools have reached equilibrium (check code for finer details).

Changing the model default parameters for user defined ones is trivial and utilises a python dictionary, e.g. ::

    from gday import adjust_gday_param_file as ad
    replace_dict = { "albedo": "0.123" }
    ad.adjust_param_file(cfg_fname, replace_dict)

And finally running the model... ::

    from gday import gday as model
    G = model.Gday(cfg_fname)
    G.run_sim()

In all cases ``cfg_fname`` is simple a string with a link to the users parameter file, it can be named anything you please.
    
Link to me
============
`My home page -> https://sites.google.com/site/mdekauwe`