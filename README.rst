=============================================
GDAY model
=============================================

** NOTE discontinued, see the C version **

GDAY (Generic Decomposition And Yield) is a simple, daily time step ecosystem model that represents carbon, nitrogen, and water dynamics at the stand scale. The model is coded entirely in `Python <http://www.python.org/>`_ without any dependancies. 


Key References
==============
1. Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO2 Enrichment; Equilibrium Behavior of Plant-Soil Models. *Ecological Applications*, 3, 666-681.
2. Medlyn, B. E., McMurtrie, R. E., Dewar, R. C. and Jeffreys, M. P. (2000), Soil processes dominate the long-term response of forest net primary productivity to increased temperature and atmospheric CO2 concentration, *Canadian Journal of Forest Research*, 30, 873–888.

**Note** there are many subtle changes from those original papers included in the code.



.. contents:: :local:

Installation
=============
Setting up python on your system is very easy. For window or mac users the `Enthought <http://www.enthought.com/>`_ or `Anaconda <http://continuum.io/downloads>`_ python packages are perhaps your simplest avenue. On a Linux machine it is simply as case of using whatever your default package manager is, e.g. sudo apt-get install python2.7. If you are on a mac and don't want to use `Enthought <http://www.enthought.com/>`_ or `Anaconda <http://continuum.io/downloads>`_ then python comes as standard with your system, so infact you don't need to do anything! However, in my personal experience I've found that it is easier to set up your own separate working copy using a package manager such as `Macports <http://www.macports.org/>`_ or `Homebrew <http://brew.sh/>`_. I read that all the cool kids are now using the later, but personally I've had no issues with Macports. 

Once you have downloaded the source code, or clone the repository (go on...) there is a simple makefile, e.g. ::

    sudo make install

or the standard python approach ::

    sudo python setup.py install

Running the model
=================

I need to add some simple scripts (on the todo list!). But within a python script you would require the following lines... ::
    
    from gday import gday as model
    
    G = model.Gday(cfg_fname, spin_up=True)
    G.spin_up_pools()

which will spin the model up. Spin up expects a met forcing file with a 50 yrs of data, how you recycle this is up to you. The model automatically stops once the soil, plant and litter pools have reached equilibrium (check code for finer details).

Changing the model default parameters for user defined ones is trivial and utilises a python dictionary, e.g. ::

    from gday import adjust_gday_param_file as ad
    
    replace_dict = { "albedo": "0.123" }
    ad.adjust_param_file(cfg_fname, replace_dict)

And finally running the model... ::

    from gday import gday as model
    
    G = model.Gday(cfg_fname)
    G.run_sim()

In all cases ``cfg_fname`` is simple a string with a link to the users parameter file, it can be named anything you please.
    
Contacts
========
* `Martin De Kauwe <http://mdekauwe.github.io/>`_  (mdekauwe at gmail.com)
* `Belinda Medlyn <http://bio.mq.edu.au/people/person.php?user=bmedlyn>`_ (bmedlyn at bio.mq.edu.au).
