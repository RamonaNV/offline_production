#!/usr/bin/env python

#############################################################################
#
#  	General Description: 	Monitors SysChk Jobs submitted to simulation production sites
#                            - Periodically submits jobs in 'chunks' to sites based on total number of active jobs
#                            - Gathers information from all sites at the same time once a week
#                            - initiated by cron job

# Copyright: (C) 2014 The IceCube collaboration
#
# @file    $Id$
# @version $Revision$

# @date    02/10/2014
# @author  Oladipo Fadiran <ofadiran@icecube.wisc.edu>
#
#############################################################################

import os, sys
import datetime
import logging

sys.path.append("/data/user/i3filter/SQLServers_n_Clients/")
try:
    import SQLClient_Simdbs4 as dbs4
    dbs4_ = dbs4.MySQL()
    
except Exception as err:
    raise Exception("Error: %s "%str(err))


def UpdateQueues(dataset_id,TotalRunning):
    try:
        logging.info("Checking submission for %s"%dataset_id)
        
        QInfo = dbs4_.fetchall(""" select * from job j where dataset_id=%s"""%dataset_id,UseDict=True)
    
        idle = [i["queue_id"] for i in QInfo if i['status'] in ["IDLE"]]
        idle.sort()
        active = [i["queue_id"] for i in QInfo if i['status'] not in ["IDLE","OK","ERROR","EVICTED","SUSPENDED","FAILED"]]
        
        if len(idle) and len(active)<TotalRunning:
            ToSubmit = TotalRunning - len(active)
            
            if len(idle) < ToSubmit:
                ToSubmit = len(idle)
            
            logging.info("""active: %s RequiredQLength: %s idle: %s
                     %s new jobs will be queued"""%(len(active),TotalRunning,len(idle),ToSubmit))
            

            SubmitQId = [str(idle[s]) for s in range(ToSubmit)]
            SubmitQId = ",".join(SubmitQId)
            
            dbs4_.execute(""" update job j set j.status="WAITING" where dataset_id=%s and queue_id in (%s)"""%(dataset_id,SubmitQId))

            
        else:
            logging.info("""active: %s RequiredQLength: %s idle: %s
                     no new jobs will be queued"""%(len(active),TotalRunning,len(idle)))
            
    
    except Exception as err:
        logging.error("failed to update queue for dataset %s"%dataset_id)
        logging.error("Error: %s "%str(err))
        #raise Exception("Error: %s "%str(err))


def GatherStats(dataset_id):
    try:
        logging.info("\nGathering stats for dataset %s"%dataset_id)
        
        Site = dbs4_.fetchall("""SELECT reqs FROM task_def t
                         where dataset_id=%s and name="MainTask" """%dataset_id,UseDict=True)
        
        Site = Site[0]['reqs'].split("'")[1]
        
        Jobs = dbs4_.fetchall("""SELECT  * FROM job j
                              where dataset_id=%s """%dataset_id,UseDict=True)
        
        ErrorJobs = [j for j in Jobs if j['status'] in ["ERROR","EVICTED","SUSPENDED","FAILED"] ]
        CompletedJobs = [j for j in Jobs if j['status'] in ["OK","ERROR","EVICTED","SUSPENDED","FAILED"]]
        SuccesfulJobs = [j for j in Jobs if j['status'] in ["OK"]]
        
        
        # rate of all jobs that reached a terminal state
        CompletionRate = 100.*float(len(CompletedJobs))/float(len(Jobs))
        
        # rate of all jobs that terminated succesfully
        SuccesfulCompletionRate = 100.*float(len(SuccesfulJobs))/float(len(Jobs))
        
        # ratio of jobs terminating unsuccesfully to those that terminated 
        if len(CompletedJobs):
            ErrorRate = 100.*float(len(ErrorJobs))/float(len(CompletedJobs))
        else:
            ErrorRate = 0
            
        
        ProblemNodeID = ""
        ProblemNodeName = ""
        MaxFailures = 0
        if ErrorRate > 0:
            MaxFailures = 1
            ProblemNode = dbs4_.fetchall("""SELECT t.host,ns.name,ns.domain,count(t.host) as count_problem_hosts
                                           FROM task t join job j on t.job_id=j.job_id
                                           left join node_statistics ns on t.host=ns.host_id
                                           where j.dataset_id=%s
                                           and t.status in ("ERROR","EVICTED","SUSPENDED","FAILED")
                                           group by t.host order by count_problem_hosts desc
                                         """%(dataset_id),UseDict=True)
            
            NodeCount = [n['count_problem_hosts'] for n in ProblemNode] 
            
            if not len(NodeCount):
                pass
            else:
                if(len(NodeCount) == 1 or max(NodeCount) > min(NodeCount)):
                    ProblemNodeID = ProblemNode[0]['host']
                    ProblemNodeName = ProblemNode[0]['name']
                    MaxFailures = ProblemNode[0]['count_problem_hosts']

        
        
        UniqueNodes = 0
        UniqueNodes = dbs4_.fetchall(""" SELECT count(distinct(t.host)) as count_hosts
                                     FROM task t join job j on t.job_id=j.job_id
                                     where j.dataset_id=%s 
                                    """%(dataset_id),UseDict=True)
        
        UniqueNodes = UniqueNodes[0]['count_hosts']

        SiteSpeed = dbs4_.fetchall(""" SELECT avg(js.value) FROM job_statistics js
                                       join job j on js.dataset_id=j.dataset_id
                                       where js.dataset_id=%s and js.name="cpu_speed"
                                       and j.dataset_id=%s and j.status="OK"
                                       """%(dataset_id,dataset_id),UseDict=True)
        
        AvgCPUSpeed = 0
        if SiteSpeed[0]['avg(js.value)'] is not None : AvgCPUSpeed = SiteSpeed[0]['avg(js.value)']
        
        return [Site, CompletionRate,SuccesfulCompletionRate,ErrorRate,AvgCPUSpeed,UniqueNodes,ProblemNodeID,ProblemNodeName,MaxFailures]
        
    except Exception as err:
        raise Exception("Error: %s "%str(err))

def UpdateSummaryTable(CheckTime,dataset,SummaryInfo):
    logging.info("\nupdating summary DB stats for site:%s, dataset_id:%s"%(SummaryInfo[0],dataset))
   
    try:
        dbs4_.execute( """ insert into syschk
                      (site,syschk_time,dataset_id,completion_rate,success_rate,error_rate,
                       avg_cpu_speed,unique_nodes,problem_node_id,problem_node_name,max_failures)
                       values("%s","%s",%s,%s,%s,%s,%s,%s,"%s","%s",%s)
                         """%(SummaryInfo[0],CheckTime,dataset,SummaryInfo[1],
                              SummaryInfo[2],SummaryInfo[3],SummaryInfo[4],
                               SummaryInfo[5],SummaryInfo[6],SummaryInfo[7],SummaryInfo[8]))
    except Exception as err:
        raise Exception("Error: %s "%str(err))

def ClearAll(dataset_id):
    
    try:
        logging.info("Cleaning %s"%dataset_id)
        
        dbs4_.execute("""delete t FROM job j  join task t on j.job_id=t.job_id
                         where j.dataset_id=%s"""%dataset_id)
        
        dbs4_.execute(""" update job j set j.status="IDLE",j.failures=0,j.errormessage=NULL
                          where dataset_id=%s """%dataset_id)
        
        dbs4_.execute(""" update grid_statistics g set suspend=0 where dataset_id=%s """%dataset_id)

    except Exception as err:
        raise Exception("Error: %s "%str(err))


TotalRunning = 100

SysChkDatasets = [10551,10552,10553,10554,10555,10556,10558,10566]
#SysChkDatasets = [10551,10552,10553,10554,10555,10556,10557,10558]

CheckTime = datetime.datetime.now()
logging.info("\n=========== %s ==========="%CheckTime)

# gather stats and reset jobs between 5 and 6am on Friday
if datetime.date.today().isoweekday() == 5:
    if CheckTime.hour >=5 and CheckTime.hour < 6:
        logging.info("\n ==== Will attempt to gather current site status, update DB, and reset jobs ====")
        SummaryInfo = {}
        for dataset in SysChkDatasets:
            SummaryInfo[dataset] = GatherStats(dataset)
            
        for dataset in SysChkDatasets:
            UpdateSummaryTable(CheckTime,dataset,SummaryInfo[dataset])
    
    SysChkDatasets = [10551,10552,10553,10554,10555,10556,10558,10566]
    #SysChkDatasets = [10551,10552,10553,10554,10555,10556,10557,10558]
    
    CheckTime = datetime.datetime.now()
    logging.info("\n=========== %s ==========="%CheckTime)
    
    # gather stats and reset jobs between 5 and 6am on Friday
    #if datetime.date.today().isoweekday() == 5:
    #if True:        # now every day
    if not CheckTime.day%2:     # check every even day
        #if True:
        if CheckTime.hour >= 7 and CheckTime.hour < 8:
            logging.info("\n ==== Will attempt to gather current site status, update DB, and reset jobs ====")
            SummaryInfo = {}
            for dataset in SysChkDatasets:
                SummaryInfo[dataset] = GatherStats(dataset)
                
            for dataset in SysChkDatasets:
                UpdateSummaryTable(CheckTime,dataset,SummaryInfo[dataset])
        
            for dataset in SysChkDatasets:
                ClearAll(dataset)
    
    
    # was needed when we needed to submit 1000 jobs in 100 chunks
    # no longer needed since only 100 jobs need to be submitted
    for dataset in SysChkDatasets:  
        UpdateQueues(dataset,TotalRunning)
    



