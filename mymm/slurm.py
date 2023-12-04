import re, os
import mymm
import subprocess

class Slurm:
    def __init__(self, account=None,nodes=1, ntasks=1, cpusPerTask=1,partition="normal",time="01:00:00",job_name_base=""):
        self.nodes  = str(nodes)
        self.ntasks = str(ntasks)
        self.cpusPerTask = str(cpus_per_task)
        self.partition = partition
        self.time = time
        self.account = account
        self.job_name_base = job_name_base
        self.job_list = []

    def add_job(self, name="", directory="", output_fname="", error_fname="", command=""):
        self.job_list.append(SlurmJob())

    def write_all_jobs(self, name=""):
        self.partition= 1
