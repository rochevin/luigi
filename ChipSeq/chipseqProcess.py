import os
import luigi
from NGS_process import *

#Execute cutadapt with luigi
#Dependencies : os, luigi, Cutadapt(from NGS_process)
class LCutadapt(Cutadapt,luigi.Task):

	#args :
	sample = luigi.Parameter() #Name of .fastq file : use by self.run()
	adapter = luigi.Parameter() #Sequence given to cutadapt who must be cut for each read

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.cutadapt.fastq" % self.sample)

	#method : running by defaut when LCutadapt is called by luigi.
	#Call cutadapt_command from Cutadapt class
	def run(self):
		file_input = self.sample+".fastq"
		file_output = self.sample+".cutadapt_prev.fastq"
		if os.path.isfile(file_output):
			file_input = file_output
			file_output = self.sample+".cutadapt.fastq"

		tmp = self.cutadapt_command(adapter = self.adapter,file_output = file_output,file_input = file_input)

		
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
	threads = luigi.Parameter() #Number of threads : use by self.run()
	adapter = luigi.Parameter() #Adapter sequence use by LCutadapt.run() for cut the given sequence for each read


	#method : self.requires() method will be launch before self.run()
	#Launch LCutadapt for each .fastq file in sample_list 
	#Launch LBwaIndex for the genome ref (.fa) if self.is_indexing() is False
	def requires(self):
		if self.is_indexing(self.genome):
			return {"cutadapt":[LCutadapt(sample,adapter) for sample in self.sample_list for adapter in self.adapter]}
		else:
			return {"cutadapt":[LCutadapt(sample,adapter) for sample in self.sample_list for adapter in self.adapter],"bwa_index":LBwaIndex(self.genome)}
			
	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.sai" % self.sample)

	#method : give the path of all fastq after cutadapt
	def get_input_path(self):
		return [sample_file.path for sample_file in self.input()["cutadapt"]]

	#method : running by defaut when LBwaAln is called by luigi.
	#Call bwa_command from BwaAln class
	def run(self):
		tmp = self.bwa_aln(sample_list = self.get_input_path(),output = self.output().path,genome = self.genome, threads = self.threads)


#Execute bwa aln with luigi
#Dependencies : luigi, BwaSam(from NGS_process), LBwaAln
class LBwaSam(luigi.Task,BwaCommand):

	#args :
	sample = luigi.Parameter() #Final name for sai/sam file : use by self.run() to create the .sai name with all sample_list merged (at least one .fastq)
	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.run() and LbwaIndex.run()  
	sample_list = luigi.Parameter() #All .fastq file use by bwa for alignment : use by self.run and LCutadapt.run()
	threads = luigi.Parameter() #Number of threads : use by self.run()
	adapter = luigi.Parameter() #Adapter sequence use by LCutadapt.run() for cut the given sequence for each read

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.sam" % self.sample)

	#method : self.requires() method will be launch before self.run()
	#Launch LBwaAln before self.run() to get self.sample+".sai" file
	def requires(self):
		return LBwaAln(self.sample,self.genome,self.sample_list,self.threads,self.adapter)

	#method : get the fastq file produce by cutadapt
	def get_fastq_list(self):
		return self.requires().get_input_path()

	#method : running by defaut when LBwaSam is called by luigi.
	#Call bwa_command from BwaSam class
	def run(self):
		tmp = self.bwa_samse(sample_list = self.get_fastq_list(),genome = self.genome,file_input = self.input().path,file_output = self.output().path)
	
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
	threads = luigi.Parameter() #Number of threads : use by self.requires()
	adapter = luigi.Parameter() #Adapter sequence use by LCutadapt.run() for cut the given sequence for each read : use by self.requires()

	#method : self.requires() method will be launch before self.run()
	#Launch LSamtools_index_ref and LBwaSam before self.run() to get self.sample+".sam" and self.genome+".fai" files
	def requires(self):
		if self.is_indexing(self.genome):
			return {"bwa_sam":LBwaSam(self.sample,self.genome,self.sample_list,self.threads,self.adapter)}
		else:
			return {"samtools_index_ref":LSamtools_index_ref(self.genome),"bwa_sam":LBwaSam(self.sample,self.genome,self.sample_list,self.threads,self.adapter)}

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
	threads = luigi.Parameter() #Number of threads : use by self.requires()
	adapter = luigi.Parameter() #Adapter sequence use by LCutadapt.run() for cut the given sequence for each read : use by self.requires()

	#method : self.requires() method will be launch before self.run()
	#Launch LSamtools_sam_to_bam to get self.sample+".bam"
	def requires(self):
		return LSamtools_sam_to_bam(self.sample,self.genome,self.sample_list,self.threads,self.adapter)
	
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
	threads = luigi.Parameter() #Number of threads : use by self.requires()
	adapter = luigi.Parameter() #Adapter sequence use by LCutadapt.run() for cut the given sequence for each read : use by self.requires()

	#method : self.requires() method will be launch before self.run()
	#Launch LSamtools_sort_bam to get self.sample+".sorted.bam"
	def requires(self):
		return LSamtools_sort_bam(self.sample,self.genome,self.sample_list,self.threads,self.adapter)

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
	threads = luigi.Parameter() #Number of threads : use by self.requires()
	adapter = luigi.Parameter() #Adapter sequence use by LCutadapt.run() for cut the given sequence for each read : use by self.requires()

	#method : self.requires() method will be launch before self.run()
	#Launch LSamtools_sort_bam to get self.sample+".nodups.bam"
	def requires(self):
		return LSamtools_remove_PCR_duplicate(self.sample,self.genome,self.sample_list,self.threads,self.adapter)

	#method : give the output file after running the command self.run()
	def output(self):
		return luigi.LocalTarget("%s.nodups.bai" % self.sample)

	#method : running by defaut when LSamtools_remove_PCR_duplicate is called by luigi.
	#Call samtools_remove_PCR_duplicate from SamtoolsCommand class
	def run(self):
		tmp = self.samtools_index_bam(bam_file = self.input().path,index_bam_file = self.output().path)

#Main class : Launch the pipeline
#Execute LBwaSam for each final_sample_name
class ChipSeqProcess(luigi.Task):
	#args : 
	file_requirement = luigi.Parameter() #The file name in .json format with all requirements for the pipeline

	#method : self.requires() method will be launch before self.run()
	#Call LBwaSam before running self.run()
	def requires(self):
		requirements = json_read(self.file_requirement) 
		final_samples_name = requirements["final_sample_name"]
		sample_list = requirements["sample_list"]
		adapter = requirements["adapter"]
		threads = requirements["threads"]
		genome = requirements["genome"]

		return [LSamtools_index_bam(sample,genome,sample_list[sample],threads,adapter) for sample in final_samples_name]
	#method : running by defaut when LBwaSam is called by luigi.
	#Use just for launching pipeline, no need to execute anything
	def run(self):
		print("running ...")

if __name__=='__main__':
	
	luigi.run()

