# simple gui for filename loading


from gday import gday as model
from Tkinter import *
from tkFileDialog import askopenfilename
import os

def main():
    
    filename = askopenfilename(filetypes=[("G'DAY config files","*.cfg"), 
                                ("All files","*")])
    
    fname = os.path.basename(filename).split(".")[0]
    fdir = os.path.dirname(filename)
    
    G = model.Gday(fname=fname, default_dir=fdir)
    G.run_sim()
    






if __name__ == "__main__":
    
    main()
    