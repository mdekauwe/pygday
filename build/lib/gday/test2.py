from gday import gday as model

def main1():
    fname = "/Users/mdekauwe/src/python/pygday/params/duke_testing.cfg"
    G = model.Gday(cfg_fname)
    G.run_sim()
    
def main2():
    fname = "/Users/mdekauwe/src/python/pygday/params/duke_testing.cfg"
    G = model.Gday(cfg_fname)
    G.run_sim()
    
    
if __name__ == "__main__":
    main1()
    main2()