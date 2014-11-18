import subprocess

def info(x):
	gitArg=["git"]
	gitArg.extend(x)
	git_info=subprocess.check_output(gitArg)
	return(git_info.split('\n'))

def URL_info(): return(info(("remote","show","origin")))