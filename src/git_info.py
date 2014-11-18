"""
git revision information. Provides short revision code, tahs and URL
for any generated files/model output
"""

__author__  = "Douglas Kelley"
__version__ = "is what this returns"
__email__   = "douglas.kelley@mq.edu.au"

revision_code	=	subprocess.check_output(["git", "describe"])
URL_Fetch		=	URL_info.split('\n')[1]
URL_pull		= 	URL_git_info()[2]
branch			=	URL_git_info()[3]
remote_branch	= 	URL_git_info()[4]
tag				= 	subprocess.check_output(["git", "tag"]).split('\n')[-2]


def URL_git_info():
    URL_info=subprocess.check_output(["git", "remote","show","origin"])
    return(URL_info.split('\n'))