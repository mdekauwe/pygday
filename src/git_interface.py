"""
git revision information. Provides short revision code, tahs and URL
for any generated files/model output
"""

__author__  = "Douglas Kelley"
__version__ = "is what this returns"
__email__   = "douglas.kelley@mq.edu.au"

import subprocess

def info(x):
    gitArg = ["git"]
    gitArg.extend(x)
    git_info = subprocess.check_output(gitArg)
    return(git_info.split('\n'))

def URL_info():
    return(info(("remote", "show", "origin")))


if __name__ == "__main__":

    import git_interface as git

    revision_code =     git.info(["rev-parse","HEAD"])[0]
    URL_Fetch = git.URL_info()[1]
    URL_push = git.URL_info()[2]
    branch = git.URL_info()[3]
    remote_branch = git.URL_info()[4]
    #tag = git.info(["tag"])[-2]

    print revision_code
    print URL_Fetch
    print URL_push
    print branch
    print remote_branch
    #print tag
