### django related packages
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.urls import reverse

from .config import *

#
from . import utils

# Create your views here.
def form(request):
    '''
    html form for the source to be analyzed
    '''
    instrumentForm = utils._loadinstruments_()
    return render(request, 'multilambda/form.html', {'instrumentForm': instrumentForm})

def submit(request):
    ### dev version to show the submit page ###
    # request_post = utils._parserequest_(request)
    # request_post = request.POST
    # return render(request, 'multilambda/submit.html', {'request_post':request_post})
    ### normal version for submission
    request_info = utils._parserequest_(request)
    taskname = utils._assigntask_(request_info['coord'])
    utils._savejson_(request_info, taskname, 'request.info')
    ## open a backgroud thread ##
    ## warning!!! this may cause DoS attack ##
    utils._openthread_(
        utils._sourceanalysis_,
        args = (request_info['coord'], taskname, request_info['todownload'])
    )
    return HttpResponseRedirect(reverse('multilambda:status', args=(taskname,)))

def status(request, taskname):
    statusjson = utils._checkdownloadstatus_(taskname)
    return render(request, 'multilambda/status.html', {'taskname':taskname, 'statusjson':statusjson})
    
def webpage(request, taskname):
    requestjson = utils._loadrequestjson_(taskname, filenameconvert=True)
    ra, dec = requestjson['coord']
    overlaypath = utils._getmulticutoutpath_(taskname)
    isvast = utils._checkmeasure_(taskname)
    return render(
        request, 'multilambda/webpage.html', 
        {
            'taskname':taskname,
            'ra':ra, 'dec':dec,
            'overlaypath': overlaypath,
            'isvast': isvast,
        }
    )