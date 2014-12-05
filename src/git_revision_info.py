"""
git revision information. Provides short revision code, tahs and URL
for any generated files/model output
"""

import git

__author__  = "Douglas Kelley"
__version__ = "is what this returns"
__email__   = "douglas.kelley@mq.edu.au"

revision_code	=	git.info(["rev-parse","HEAD"])[0]
URL_Fetch		=	git.URL_info()[1]
URL_push		= 	git.URL_info()[2]
branch			=	git.URL_info()[3]
remote_branch	= 	git.URL_info()[4]
tag				= 	git.info(["tag"])[-2]