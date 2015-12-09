import os
import luigi
from NGS_process import *
from Luigi_process import *

#Launch bwa aln with luigi,
#dependencies : LBwaAln
class BwaAln(LBwaAln):

	def requires(self):
		if not self.is_indexing(self.genome):
			LBwaIndex(self.genome)
		else:
			None
			
#Launch bwa samse with luigi,
#dependencies : LBwaSamse
class BwaSamse(LBwaSamse):

	def requires(self):
		return BwaAln(self.sample,self.genome,self.sample_list)

#Launch samtools view with luigi,
#dependencies : LBwaSamse
class SamToBam(LSamtools_sam_to_bam):

	def requires(self):
		if self.is_indexing(self.genome):
			return {"bwa_sam":BwaSamse(self.sample,self.genome,self.sample_list)}
		else:
			return {"samtools_index_ref":LSamtools_index_ref(self.genome),"bwa_sam":BwaSamse(self.sample,self.genome,self.sample_list)}

#Launch samtools sort with luigi,
#dependencies : LSamtools_sort_bam
class SortBam(LSamtools_sort_bam):

	def requires(self):
		return SamToBam(self.sample,self.genome,self.sample_list)


class MergeBam(LSamtools_merge_bam):

	genome = luigi.Parameter() #Reference genome (.fa format) : use by self.requires()

	def requires(self):
		return [SortBam(fastq_file,self.genome,[fastq_file]) for fastq_file in self.sample_list]


class RemovePCRDuplicate(LSamtools_remove_PCR_duplicate):

	def requires(self):
		return MergeBam(self.sample,self.sample_list,self.genome)
	
class IndexBam(LSamtools_index_bam):

	def requires(self):
		return RemovePCRDuplicate(self.sample,self.genome,self.sample_list)


#Main class : Launch the pipeline
#Execute Pipeline for each final_sample_name
class ChipSeqProcess(luigi.Task):
	#args : 
	file_requirement = luigi.Parameter() #The file name in .json format with all requirements for the pipeline

	#method : self.requires() method will be launch before self.run()
	#Call LBwaSam before running self.run()
	def requires(self):
		requirements = json_read(self.file_requirement) 
		final_samples_name = requirements["final_sample_name"]
		sample_list = requirements["sample_list"]
		genome = requirements["genome"]

		return [IndexBam(sample,genome,sample_list[sample]) for sample in final_samples_name]
	#method : running by defaut when LBwaSam is called by luigi.
	#Use just for launching pipeline, no need to execute anything
	def run(self):
		print("running ...")
