#!/bin/bash
#set a job name  
#SBATCH --job-name=run1
#################  
#a file for job output, you can check job progress
#SBATCH --output=run1.out
#################
# a file for errors from the job
#SBATCH --error=run1.err
#################
#time you think you need; default is one day
#in minutes in this case, hh:mm:ss
#SBATCH --time=30:00:00
#################
#number of tasks you are requesting
#SBATCH -n 32
#SBATCH --exclusive
#################
#partition to use
#SBATCH --partition=ht-10g
#################
#number of nodes to distribute n tasks across
#SBATCH -N 32
#################

work=/home/li.yi3/081816_doubletaper

cd $work

mpirun -np 1 R --no-save < doubletaper_exp_simu.R
