import os
import luigi
from NGS_process import *

#Execute cutadapt with luigi
#Dependencies : os, luigi, Cutadapt(from NGS_process)
class LCutadapt(luigi.Task):

	#args :
	sample = luigi.Parameter() #Name of .fastq file : use by self.run()
	adapter = luigi.Parameter() #Sequence given to cutadapt who must be cut for each read

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.cutadapt.fastq" % self.sample)

	#method : give the name of the input file use by cutadapt  
	def sample_input(self):
		#if cutadapt was already use for the self.sample name, the command will be execute for the self.output().path instead of self.sample+".fastq" 
		cutadapted = ".cutadapt" if os.path.isfile(self.output().path) else ""
		return luigi.LocalTarget("%s%s.fastq" % (self.sample,cutadapted))
	#method : running by defaut when LCutadapt is called by luigi.
	#Call cutadapt_command from Cutadapt class
	def run(self):
		file_input = self.sample+".fastq"
		file_output = self.sample+".cutadapt_prev.fastq"
		if os.path.isfile(file_output):
			file_input = file_output
			file_output = self.sample+".cutadapt.fastq"

		tmp = Cutadapt.cutadapt_command(self.adapter,file_output,file_input)

#Execute bwa index with luigi
#Dependencies : os, luigi, BwaIndex(from NGS_process)
class LBwaIndex(luigi.Task,BwaCommand):

	#args :
	genome = luigi.Parameter(default="hg19.fa") #the name of reference genome 
	#method : running by defaut when LBwaIndex is called by luigi.
	#Call bwa_command from BwaIndex class
	def run(self):
		tmp = self.bwa_index(self.genome)

#Execute bwa aln with luigi
#Dependencies : os, luigi, BwaAln(from NGS_process), LCutadapt, LBwaIndex (if no .fai)
class LBwaAln(luigi.Task,BwaCommand):

	#args :
	sample = luigi.Parameter() #Final name for sai/sam file : use by self.run() to create the .sai name with all sample_list merged (at least one .fastq)
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.run() and LbwaIndex.run()  
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.run and LCutadapt.run()

		
	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.sai" % self.sample)

	#method : running by defaut when LBwaAln is called by luigi.
	#Call bwa_command from BwaAln class
	def run(self):
		tmp = self.bwa_aln(fastq_list = self.sample_list,output = self.output().path,genome = self.genome)


#Execute bwa aln with luigi
#Dependencies : luigi, BwaSam(from NGS_process), LBwaCommand
class LBwaSamse(luigi.Task,BwaCommand):

	#args :
	sample = luigi.Parameter() #Final name for sai/sam file : use by self.run() to create the .sai name with all sample_list merged (at least one .fastq)
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.run() and LbwaIndex.run()  
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.run and LCutadapt.run()

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.sam" % self.sample)

	#method : running by defaut when LBwaSam is called by luigi.
	#Call bwa_command from BwaSam class
	def run(self):
		tmp = self.bwa_samse(fastq_list = self.sample_list,genome = self.genome,file_input = self.input().path,file_output = self.output().path)

#Execute bwa aln with luigi
#Dependencies : luigi, BwaSam(from NGS_process), LBwaCommand
class LBwaSampe(luigi.Task,BwaCommand):

	#args :
	sample = luigi.Parameter() #Final name for sai/sam file : use by self.run() to create the .sai name with all sample_list merged (at least one .fastq)
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.run() and LbwaIndex.run()  
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.run

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.sam" % self.sample)

	#method : running by defaut when LBwaSam is called by luigi.
	#Call bwa_command from BwaSam class
	def run(self):
		tmp = self.bwa_sampe(sai_list=self.sample_list),fastq_list=self.sample_list,genome=self.genome,file_output=self.output().path)
	
#Execute samtools faidx with luigi
#Dependencies : luigi, SamtoolsCommand(from NGS_process)		
class LSamtools_index_ref(luigi.Task,SamtoolsCommand):

	#args :
	genome = luigi.Parameter(default="hg19.fa") #Reference genome (.fa format) : use by self.run()

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.fai" % self.genome)

	#method : running by defaut when LSamtools_index_ref is called by luigi.
	#Call samtools_index_ref from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_index_ref(self.genome)

#Execute samtools view with luigi
#Dependencies : luigi, SamtoolsCommand(from NGS_process), LSamtools_index_ref, LBwaSam
class LSamtools_sam_to_bam(luigi.Task,SamtoolsCommand):

	#args :
	sample = luigi.Parameter() #Final name for bam file : use by self.run() to create the .bam name with self.sample
	genome = luigi.Parameter() #Reference genome (.fa format) : use by requires.run() especially for LSamtools_index_ref
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.requires()

	#method : give the indexed genome file to convert sam to bam
	#if the LSamtools_index_ref is launch, return his output, else return the self.genome+".fai" path
	def index_input(self):
		return self.input()["samtools_index_ref"] if "samtools_index_ref" in self.input() else luigi.LocalTarget("%s.fai" % self.genome)

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.bam" % self.sample)

	#method : running by defaut when LSamtools_sam_to_bam is called by luigi.
	#Call samtools_sam_to_bam from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_sam_to_bam(file_index = self.index_input().path,file_output = self.output().path,file_input = self.input()["bwa_sam"].path)

#Execute samtools sort with luigi
#Dependencies : luigi, SamtoolsCommand(from NGS_process), LSamtools_sam_to_bam
class LSamtools_sort_bam(luigi.Task,SamtoolsCommand):
	#args :
	sample = luigi.Parameter() #Final name for bam file : use by self.run() to create the .sorted.bam
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.requires()
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.requires()

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.sorted.bam" % self.sample)

	#method : running by defaut when LSamtools_sort_bam is called by luigi.
	#Call samtools_sort_bam from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_sort_bam(bam_file = self.input().path,sorted_bam_file = self.output().path)


#Execute samtools sort with luigi
#Dependencies : luigi, SamtoolsCommand(from NGS_process), LSamtools_sort_bam
class LSamtools_remove_PCR_duplicate(luigi.Task,SamtoolsCommand):
	#args :
	sample = luigi.Parameter() #Final name for bam file : use by self.run() to create the .nodups.bam name with all sample_list merged (at least one .fastq)
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.requires()
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.requires()

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.nodups.bam" % self.sample)

	#method : running by defaut when LSamtools_remove_PCR_duplicate is called by luigi.
	#Call samtools_remove_PCR_duplicate from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_remove_PCR_duplicate(bam_file = self.input().path,nodups_bam_file = self.output().path)

#Execute samtools index with luigi
#Dependencies : luigi, SamtoolsCommand(from NGS_process), LSamtools_sort_bam
class LSamtools_index_bam(luigi.Task,SamtoolsCommand):
	#args :
	sample = luigi.Parameter() #Final name for bam file : use by self.run() to create the .nodups.bam name with self.sample
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.requires()
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.requires()

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.nodups.bam" % self.sample)

	#method : running by defaut when LSamtools_remove_PCR_duplicate is called by luigi.
	#Call samtools_index_bam from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_index_bam(bam_file = self.input().path,index_bam_file = self.output().path)

class LSamtools_merge_bam(luigi.Task,SamtoolsCommand):
	#args :
	sample = luigi.Parameter() #Final name for bam file : use by self.run() to create the .merged.bam with self.sample
	sample_list = luigi.Parameter() #All bam file merged into one .merged.bam file
	
	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.bam" % self.sample)

	#method : running by defaut when LSamtools_remove_PCR_duplicate is called by luigi.
	#Call samtools_merge_bam from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_merge_bam(merged_bam_file=self.output().path,bam_list=self.bam_for_sample(self.sample_list))

